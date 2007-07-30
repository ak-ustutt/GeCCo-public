*----------------------------------------------------------------------*
      subroutine add_opblk(fac,ffin,ffout,
     &     opin,opout,iblkin,iblkout,orb_info)
*----------------------------------------------------------------------*
*
*     add block from file ffin to file ffout
*     occupation of blocks must be identical
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(in) ::
     &     fac
      type(filinf), intent(inout) ::
     &     ffin, ffout
      type(operator), intent(in) ::
     &     opin, opout
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkin, iblkout

      logical ::
     &     bufin, bufout
      integer ::
     &     len_op, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idxst, idxnd,
     &     idoffin, idoffout
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      logical, external ::
     &     iocc_equal, irestr_equal

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' add_opblk messing around'
        write(luout,*) '=========================='
        write(luout,*) ' fac = ',fac
        write(luout,*) ' ffin:  ',trim(ffin%name),
     &                   ' rec: ',ffin%current_record
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opin: ',opin%name(1:len_trim(opin%name)),
     &       '  block: ',iblkin
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name)),
     &       '  block: ',iblkout
      end if

      if (.not.iocc_equal(opin%ihpvca_occ(1,1,iblkin),opin%dagger,
     &                 opout%ihpvca_occ(1,1,iblkout),opout%dagger)) then
        write(luout,*) 'dagger: ',opin%dagger,opout%dagger
        call wrt_occ(luout,opin%ihpvca_occ(1,1,iblkin))
        call wrt_occ(luout,opout%ihpvca_occ(1,1,iblkout))
        call quit(1,'add_opblk','occupations do not fit!')
      end if

      if (.not.irestr_equal(opin%igasca_restr(1,1,1,1,iblkin),
     &                     opin%dagger,
     &                     opout%igasca_restr(1,1,1,1,iblkout),
     &                     opout%dagger,
     &                     orb_info%ngas)) then
        write(luout,*) 'dagger: ',opin%dagger,opout%dagger
        call wrt_rstr(luout,opin%igasca_restr(1,1,1,1,iblkin))
        call wrt_rstr(luout,opout%igasca_restr(1,1,1,1,iblkout))
        call quit(1,'add_opblk','occupations do not fit!')
        ! note: we must be able to handle this case in the future
      end if

      len_op = opin%len_op_occ(iblkin)
      ! for the moment this must hold:
      if (len_op.ne.opout%len_op_occ(iblkout))
     &     call quit(1,'add_opblk','unexpected error')

      ! buffered data available?
      bufin = .false.
      bufout = .false.
      if (ffin%buffered) bufin = ffin%incore(iblkin).gt.0
      if (ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('add_opblk')

      if (.not.bufin.or..not.bufout) then

        ! hopefully, both files have same rec-length
        ! does not matter in principle, but it does for efficiency ...
        lblk = max(ffin%reclen,ffout%reclen)

        nblkmax = ifree/lblk/2
        if (nblkmax.le.0) then
          write(luout,*) 'free memory (words):  ',ifree
          write(luout,*) 'block length (words): ',lblk,' * 2'
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

      ioffin  = opin%off_op_occ(iblkin)
      ioffout = opout%off_op_occ(iblkout)

      if (.not.bufin.or..not.bufout) then

        ! offset on file (if more than one instance of operator ex.)
        idoffin  = ffin %length_of_record*(ffin %current_record-1)
        idoffout = ffout%length_of_record*(ffout%current_record-1)

        idxst = 1
        do while(idxst.le.len_op)
          idxnd = min(len_op,idxst-1+nbuff)
          if (bufin) then
            call get_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
            call daxpy(idxnd-idxst+1,fac,ffin%buffer(ioffin+idxst),1,
     &                                   buffer_out,1)
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          else if (bufout) then
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            call daxpy(idxnd-idxst+1,fac,buffer_in,1,
     &                               ffout%buffer(ioffout+idxst),1)
          else
            call get_vec(ffin,buffer_in,idoffin+ioffin+idxst,
     &                                  idoffin+ioffin+idxnd)
            call get_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
            buffer_out(1:(idxnd-idxst+1))
     &         = fac*buffer_in(1:(idxnd-idxst+1))
     &              +buffer_out(1:(idxnd-idxst+1))
            call put_vec(ffout,buffer_out,idoffout+ioffout+idxst,
     &                                    idoffout+ioffout+idxnd)
          end if
          idxst = idxnd+1
        end do
      else
        call daxpy(len_op,fac,ffin%buffer(ioffin+1),1,
     &                       ffout%buffer(ioffout+1),1)
      end if

      ifree = mem_flushmark()

      return
      end
