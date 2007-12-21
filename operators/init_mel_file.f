*----------------------------------------------------------------------*
      subroutine init_mel_file(mel,rec,rec_lo,rec_hi,incore)
*----------------------------------------------------------------------*
*     initialize file containing ME-list mel
*     rec is the current active record (set to 1 if rec<=0 on input)
*     rec_lo, rec_hi are the bounds for rec (for later calls to
*     switch_mel_record(); both set to value of rec if <=0 on input)
*     incore==0, no incore buffers
*     incore==1, incore buffers for rank one operators
*     incore>1 not supported (I suggest to handle the incore stuff
*         in frm_schedX() when we know how much core memory can be
*         spent for this purpose)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'par_globalmarks.h'
      include 'par_opnames_gen.h'
      include 'ifc_memman.h'
      include 'ioparam.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'

      integer, intent(in) ::
     &     incore, rec, rec_lo, rec_hi
      type(me_list), intent(inout), target ::
     &     mel

      integer ::
     &     ifree, nbuff, iocc, nrecords
      type(filinf), pointer ::
     &     fhand
      type(operator), pointer ::
     &     op

      op => mel%op

      if (associated(mel%fhand))
     &     call quit(1,'init_mel_file',
     &     'file for ME-list is already assigned: '//
     &     trim(mel%label))

      allocate(mel%fhand)
      fhand => mel%fhand
      ! standard name
      call file_init(fhand,
     &       trim(mel%label)//'_list.da',ftyp_da_unf,
     &       lblk_da)

      if (rec.gt.0) then
        fhand%current_record = rec
      else if (rec_lo.gt.0) then
        fhand%current_record = rec_lo
      else
        fhand%current_record = 1
      end if

      ! set super-record length: must be multiple of primitive reclen
      fhand%length_of_record =
     &     ((mel%len_op-1)/fhand%reclen+1)*fhand%reclen

      if (rec_lo.gt.0) then
        if (rec_lo.gt.fhand%current_record) then
          write(luout,*) 'rec_lo, rec: ',rec_lo, rec
          call quit(1,'init_mel_file','record bounds inconsistent')
        end if
        fhand%active_records(1) = rec_lo
      else
        fhand%active_records(1) =
     &       fhand%current_record
      end if

      if (rec_hi.gt.0) then
        if (rec_hi.lt.fhand%current_record.or.rec_hi.lt.rec_lo) then
          write(luout,*) 'rec_lo, rec, rec_hi: ',rec_lo, rec, rec_hi
          call quit(1,'init_mel_file','record bounds inconsistent')
        end if
        fhand%active_records(2) = rec_hi
      else
        fhand%active_records(2) =
     &       fhand%current_record
      end if

      nrecords = fhand%active_records(2)
      allocate(fhand%last_mod(nrecords))
      fhand%last_mod(1:nrecords) = -1

      if (incore.eq.1) then

        fhand%buffered = .true.
        ! make sure, we are in the appropriate memory section
        call mem_pushmark()
        ifree = mem_gotomark(me_list_def)

        ifree = mem_alloc_int(fhand%incore,
     &         op%n_occ_cls,trim(mel%label)//'_incore')
        ifree = mem_alloc_int(fhand%idxrec,
     &         op%n_occ_cls,trim(mel%label)//'_idxrec')
        nbuff = 0
        do iocc = 1, op%n_occ_cls
          ! allocate buffer for 1-hamiltonian only
          if (max(op%ica_occ(1,iocc),
     &              op%ica_occ(2,iocc)).gt.1 .or.
     &         iextr.gt.0.and.  ! R12: ignore currently
     &         op%ihpvca_occ(iextr,1,iocc)+
     &         op%ihpvca_occ(iextr,2,iocc).gt.0 ) then
            fhand%incore(iocc) = -1
            fhand%idxrec(iocc) = -1
          else
            fhand%idxrec(iocc) = nbuff
            fhand%incore(iocc) = mel%len_op_occ(iocc)
            nbuff = nbuff + mel%len_op_occ(iocc)
          end if
        end do
        fhand%nbuffer = nbuff
        ifree = mem_alloc_real(fhand%buffer,
     &       nbuff,trim(mel%label)//'_buffer')

        call mem_popmark()

      else if (incore.gt.1) then
        call quit(1,'init_mel_file','incore>1 not supported!')
      end if

      if (iprlvl.ge.10) then
        write(luout,'(3x,7a)')
     &       'assigned file: ',trim(fhand%name),' to list ',
     &       trim(mel%label),' (operator: ',trim(op%name),')'
        write(luout,'(3x,a,i4,a,i4,a,i4,a)') 'record: ',
     &       fhand%current_record,
     &       ' (active: ',fhand%active_records(1),
     &       ' -- ',      fhand%active_records(2),
     &       ' )'
      end if

      return

      end
