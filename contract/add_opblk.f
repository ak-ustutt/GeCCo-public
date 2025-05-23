*----------------------------------------------------------------------*
      subroutine add_opblk(xnorm2,type,fac,mel_in,mel_out,
     &     iblkin,iblkout,orb_info,reset)
*----------------------------------------------------------------------*
*
*     add block from list mel_in to list mel_out
*     occupation of blocks must be identical
*
*     xnorm2: updated squared norm of the block
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'ifc_memman.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 000

      real(8), intent(out) ::
     &     xnorm2
      real(8), intent(in) ::
     &     fac
      type(me_list), intent(in) ::
     &     mel_in, mel_out
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkin, iblkout, type
      logical, intent(in) ::
     &     reset

      logical ::
     &     ok, bufin, bufout
      integer ::
     &     len_op, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idx, idxst, idxnd, lenbat,
     &     idoffin, idoffout,
     &     njoined_in, njoined_out, ijoin, ijoin_out,
     &     idx_in, idx_out

      integer, pointer ::
     &     occ_try(:,:,:)
      type(filinf), pointer ::
     &     ffin, ffout
      type(operator), pointer ::
     &     opin, opout

      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      logical, external ::
c     &     iocc_equal_n, iocc_equal,
     &     irestr_equal
      real(8), external ::
     &     ddot

      ffin  => mel_in%fhand
      ffout => mel_out%fhand
      opin  => mel_in%op
      opout => mel_out%op

      if (ntest.ge.100) then
        write(lulog,*) '=========================='
        write(lulog,*) ' add_opblk messing around'
        write(lulog,*) '=========================='
        write(lulog,*) ' fac = ',fac
        write(lulog,*) ' mel_in: ',trim(mel_in%label)
        write(lulog,*) ' mel_out: ',trim(mel_out%label)
        write(lulog,*) ' ffin:  ',trim(ffin%name),
     &                   ' rec: ',ffin%current_record
        write(lulog,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(lulog,*) ' opin: ',trim(opin%name),
     &       '  block: ',iblkin
        write(lulog,*) ' opout: ',trim(opout%name),
     &       '  block: ',iblkout
      end if

      njoined_in  = opin%njoined
      njoined_out = opout%njoined
      idx_in  = (iblkin-1)*njoined_in+1
      idx_out = (iblkout-1)*njoined_out+1

      if (njoined_in.eq.njoined_out) then
        ok = iocc_equal_n(opin%ihpvca_occ(1:ngastp,1:2,
     &                                    idx_in:idx_in-1+njoined_in),
     &                   .false.,
     &                   opout%ihpvca_occ(1:ngastp,1:2,
     &                                    idx_out:idx_out-1+njoined_in),
     &                   .false.,
     &                   njoined_in)
      else if (njoined_in.eq.2.and.njoined_out.eq.1) then
        allocate(occ_try(ngastp,2,1))
        occ_try(1:ngastp,1:2,1) =
     &      iocc_overlap(opin%ihpvca_occ(1:ngastp,1:2,idx_in),.false.,
     &                   opin%ihpvca_occ(1:ngastp,1:2,idx_in+1),.false.)
        ok = iocc_zero(occ_try)
        if (ok) then
          occ_try(1:ngastp,1:2,1) =
     &        iocc_add(1,opin%ihpvca_occ(1:ngastp,1:2,idx_in),.false.,
     &                 1,opin%ihpvca_occ(1:ngastp,1:2,idx_in+1),.false.)
          ok = iocc_equal(occ_try,.false.,
     &                   opout%ihpvca_occ(1:ngastp,1:2,idx_out),.false.)
        end if
        deallocate(occ_try)
      else if (njoined_out-njoined_in.eq.1) then
        ! needed ?
        ! density exception: insert a zero-occupation
        allocate(occ_try(ngastp,2,njoined_out))
        do ijoin = 1, njoined_out
          if(ijoin.gt.1)
     &         occ_try(1:ngastp,1:2,1:ijoin-1) =
     &         opin%ihpvca_occ(1:ngastp,1:2,idx_in:idx_in-1+ijoin-1)
          occ_try(1:ngastp,1:2,ijoin) = 0
          if (ijoin.lt.njoined_out)
     &         occ_try(1:ngastp,1:2,ijoin+1:njoined_out) =
     &         opin%ihpvca_occ(1:ngastp,1:2,idx_in-1+ijoin:
     &                             idx_in-1+njoined_out-1)
          ok = iocc_equal_n(occ_try,opin%dagger,
     &          opout%ihpvca_occ(1:ngastp,1:2,
     &                           idx_out:idx_out-1+njoined_out),.false.,
     &          njoined_out)
          if (ok) exit
        end do
        deallocate(occ_try)
      else
        ! ok if only differences are zero vertices
        ok = .true.
        ijoin_out = 0
        in_loop: do ijoin = 1, njoined_in
          if (iocc_zero(opin%ihpvca_occ(1:ngastp,1:2,idx_in-1+ijoin)))
     &       cycle
          ok = .false.
          do while(ijoin_out.lt.njoined_out)
            ijoin_out = ijoin_out + 1
            if (iocc_zero(opout%ihpvca_occ(1:ngastp,1:2,
     &                                     idx_out-1+ijoin_out))) cycle
            if (iocc_equal_n(opin%ihpvca_occ(1:ngastp,1:2,
     &                                       idx_in-1+ijoin),.false.,
     &                   opout%ihpvca_occ(1:ngastp,1:2,
     &                                    idx_out-1+ijoin_out),.false.,
     &                   1)) then
              ok = .true.
              cycle in_loop
            else
              exit in_loop
            end if
          end do
        end do in_loop
        if (ok) then
          do while(ijoin_out.lt.njoined_out)
            ijoin_out = ijoin_out + 1
            if (iocc_zero(opout%ihpvca_occ(1:ngastp,1:2,
     &                                     idx_out-1+ijoin_out))) cycle
            ok = .false.
            exit
          end do
        end if
      end if

      if (.not.ok) then
        write(lulog,*) 'dagger: ',opin%dagger,opout%dagger
        call wrt_occ_n(lulog,opin%ihpvca_occ(1,1,idx_in),njoined_in)
        call wrt_occ_n(lulog,opout%ihpvca_occ(1,1,idx_out),njoined_out)
        call quit(1,'add_opblk','occupations do not fit!')
      end if

c      if (.not.irestr_equal(opin%igasca_restr(1,1,1,1,1,iblkin),
c     &                     opin%dagger,
c     &                     opout%igasca_restr(1,1,1,1,1,iblkout),
c     &                     opout%dagger,
c     &                     orb_info%ngas)) then
c        write(lulog,*) 'dagger: ',opin%dagger,opout%dagger
c        call wrt_rstr(lulog,opin%igasca_restr(1,1,1,1,1,iblkin))
c        call wrt_rstr(lulog,opout%igasca_restr(1,1,1,1,1,iblkout))
c        call quit(1,'add_opblk','occupations do not fit!')
c        ! note: we must be able to handle this case in the future
c      end if

      if (opin%dagger.and..not.opout%dagger .or.
     &    .not.opin%dagger.and.opout%dagger)
     &     call quit(1,'add_opblk','cannot (yet) transpose on-the-fly!')
      ! note: there is add_opblk_transp to do this

      len_op = mel_in%len_op_occ(iblkin)
      ! for the moment this must hold:
      if (len_op.ne.mel_out%len_op_occ(iblkout))then
        write(lulog,*)'len_op = ',len_op
        write(lulog,*)'mel len = ',mel_out%len_op_occ(iblkout)
        write(lulog,*)'mel_in: ',trim(mel_in%label),iblkin
        write(lulog,*)'mel_out: ',trim(mel_out%label),iblkout
        call wrt_occ_n(lulog,opin%ihpvca_occ(1,1,idx_in),njoined_in)
        call wrt_occ_n(lulog,opout%ihpvca_occ(1,1,idx_out),njoined_out)
        write(lulog,*)'formal? ',mel_in%op%formal_blk(iblkin)
        write(lulog,*)'formal? ',mel_out%op%formal_blk(iblkout)
        call quit(1,'add_opblk','unexpected error')      
      endif

      ! buffered data available?
      bufin = .false.
      bufout = .false.
      if (ffin%buffered) bufin = ffin%incore(iblkin).ge.0
      if (ffout%buffered) bufout = ffout%incore(iblkout).ge.0

      ifree = mem_setmark('add_opblk')

      if (.not.bufin.or..not.bufout) then

        ! hopefully, both files have same rec-length
        ! does not matter in principle, but it does for efficiency ...
        lblk = max(ffin%reclen,ffout%reclen)

        nblkmax = ifree/lblk/2
        if (nblkmax.le.0) then
          write(lulog,*) 'free memory (words):  ',ifree
          write(lulog,*) 'block length (words): ',lblk,' * 2'
          call
     &       quit(1,'add_opblk','not even 1 record fits into memory?')
        end if

        nblk = min((len_op-1)/lblk+1,nblkmax)

        if (.not.bufin) then
          nbuff = min(len_op,nblk*lblk)
          ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        end if
        if (.not.bufout) then
          nbuff = min(len_op,nblk*lblk)
          ifree = mem_alloc_real(buffer_out,nbuff,'buffer_out')
        end if
      end if

      ioffin  = mel_in%off_op_occ(iblkin)
      ioffout = mel_out%off_op_occ(iblkout)

      if (.not.bufin.or..not.bufout) then

        ! offset on file (if more than one instance of operator ex.)
        idoffin  = ffin %length_of_record*(ffin %current_record-1)
        idoffout = ffout%length_of_record*(ffout%current_record-1)

        xnorm2 = 0d0
        idxst = 1
        do while(idxst.le.len_op)
          idxnd = min(len_op,idxst-1+nbuff)
          lenbat = idxnd-idxst+1
          if (bufin) then
            if (.not.reset)
     &        call get_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
            if (reset) buffer_out(1:lenbat) = 0d0
            call daxpy(lenbat,fac,ffin%buffer(ioffin+idxst),1,
     &                                   buffer_out,1)
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          else if (bufout) then
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            if (reset) ffout%buffer(ioffout+idxst:ioffout+idxnd) = 0d0
            call daxpy(lenbat,fac,buffer_in,1,
     &                               ffout%buffer(ioffout+idxst),1)
          else
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            if (.not.reset)
     &       call get_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
            if (reset) buffer_out(1:lenbat) = 0d0
            do idx = 1, lenbat
              buffer_out(idx) = buffer_out(idx)+fac*buffer_in(idx)
            end do
c does not alway work correctly (range error at end of buffer_out) ???
c            buffer_out(1:lenbat)
c     &         = fac*buffer_in(1:lenbat)
c     &              +buffer_out(1:lenbat)
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          end if
          idxst = idxnd+1
          if (type.eq.1.and..not.bufout)
     &      xnorm2 = xnorm2 + ddot(lenbat,buffer_out,1,buffer_out,1)
          if (type.eq.1.and.bufout)
     &      xnorm2 = xnorm2 + ddot(lenbat,ffout%buffer(ioffout+idxst),1,
     &                                    ffout%buffer(ioffout+idxst),1)
        end do
        if (type.eq.2.and..not.bufout)
     &      xnorm2 = buffer_out(1)
        if (type.eq.2.and.bufout)
     &      xnorm2 = ffout%buffer(ioffout+1)
      else
        if (reset) ffout%buffer(ioffout+1:ioffout+len_op) = 0d0
        call daxpy(len_op,fac,ffin%buffer(ioffin+1),1,
     &                       ffout%buffer(ioffout+1),1)
        if (type.eq.1) then
          xnorm2 = ddot(len_op,ffout%buffer(ioffout+1),1,
     &                       ffout%buffer(ioffout+1),1)
        else
          xnorm2 = ffout%buffer(ioffout+1)
        end if
      end if

      ifree = mem_flushmark()

      return
      end
