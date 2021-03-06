*----------------------------------------------------------------------*
      subroutine scale_opblk(xnorm2,fac,mel_in,mel_out,
     &     iblkin,iblkout,orb_info)
*----------------------------------------------------------------------*
*
*     scale block from list mel_in and put to list mel_out 
*     (may be identical)
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

      integer, parameter ::
     &     ntest = 00

      real(8), intent(out) ::
     &     xnorm2
      real(8), intent(in) ::
     &     fac
      type(me_list), intent(in) ::
     &     mel_in, mel_out
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkin, iblkout

      logical ::
     &     ok, bufin, bufout
      integer ::
     &     len_op, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idxst, idxnd,
     &     idoffin, idoffout,
     &     njoined_in, njoined_out, ijoin,
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
     &     iocc_equal_n, irestr_equal
      real(8), external ::
     &     ddot

      ffin  => mel_in%fhand
      ffout => mel_out%fhand
      opin  => mel_in%op
      opout => mel_out%op

      if (ntest.ge.100) then
        write(lulog,*) '============================'
        write(lulog,*) ' scale_opblk messing around'
        write(lulog,*) '============================'
        write(lulog,*) ' fac = ',fac
        write(lulog,*) ' mel_in: ',trim(mel_in%label)
        write(lulog,*) ' opout: ',trim(mel_out%label)
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
        ok = iocc_equal_n(opin%ihpvca_occ(1,1,idx_in),opin%dagger,
     &                   opout%ihpvca_occ(1,1,idx_out),opout%dagger,
     &                   njoined_in)
      else if (njoined_out-njoined_in.eq.1) then
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
     &                   opout%ihpvca_occ(1,1,idx_out),opout%dagger,
     &                   njoined_out)
          if (ok) exit
        end do
        deallocate(occ_try)
      else
        ok = .false.
      end if

      if (.not.ok) then
        write(lulog,*) 'dagger: ',opin%dagger,opout%dagger
        call wrt_occ_n(lulog,opin%ihpvca_occ(1,1,idx_in),njoined_in)
        call wrt_occ_n(lulog,opout%ihpvca_occ(1,1,idx_out),njoined_out)
        call quit(1,'scale_opblk','occupations do not fit!')
      end if

      if (opin%dagger.and..not.opout%dagger .or.
     &    .not.opin%dagger.and.opout%dagger)
     &     call quit(1,'scale_opblk',
     &     'cannot (yet) transpose on-the-fly!')

      len_op = mel_in%len_op_occ(iblkin)
      ! for the moment this must hold:
      if (len_op.ne.mel_out%len_op_occ(iblkout))
     &     call quit(1,'scale_opblk','unexpected error')      

      ! buffered data available?
      bufin = .false.
      bufout = .false.
      if (ffin%buffered) bufin = ffin%incore(iblkin).gt.0
      if (ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('scale_opblk')

      if (.not.bufin.or..not.bufout) then

        ! hopefully, both files have same rec-length
        ! does not matter in principle, but it does for efficiency ...
        lblk = max(ffin%reclen,ffout%reclen)

        nblkmax = ifree/lblk/2
        if (nblkmax.le.0) then
          write(lulog,*) 'free memory (words):  ',ifree
          write(lulog,*) 'block length (words): ',lblk,' * 2'
          call
     &       quit(1,'scale_opblk','not even 1 record fits into memory?')
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

        idxst = 1
        do while(idxst.le.len_op)
          idxnd = min(len_op,idxst-1+nbuff)
          if (bufin) then
            buffer_out(1:(idxnd-idxst+1)) =
     &           fac*ffin%buffer(ioffin+idxst:ioffin+idxnd)
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          else if (bufout) then
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            ffout%buffer(ioffout+idxst:ioffout+idxnd) =
     &           fac*buffer_in(1:(idxnd-idxst+1))
          else
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            buffer_out(1:(idxnd-idxst+1))
     &           = fac*buffer_in(1:(idxnd-idxst+1))
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          end if
          idxst = idxnd+1
        end do
        xnorm2 = ddot(len_op,buffer_out,1,buffer_out,1)
      else
        ffout%buffer(ioffout+1:ioffout+len_op) =
     &       fac*ffin%buffer(ioffin+1:ioffin+len_op)
        xnorm2 = ddot(len_op,ffout%buffer(ioffout+1),1,
     &                       ffout%buffer(ioffout+1),1)
      end if

      ifree = mem_flushmark()

      return
      end
