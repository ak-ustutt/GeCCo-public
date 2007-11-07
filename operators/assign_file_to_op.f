*----------------------------------------------------------------------*
      subroutine assign_file_to_op(idxop,standard,ffop_in,
     &                             rec,rec_lo,rec_hi,
     &                             incore,op_info)
*----------------------------------------------------------------------*
*     assign a file to operator #idxop
*     if standard==.true. allocate file-handle and generate the 
*     standard-name, else use file-handle ffop (initialized outside)
*     rec is the current active record (set to 1 if rec<=0 on input)
*     rec_lo, rec_hi are the bounds for rec (for later calls to
*     switch_opfile_record(); both set to value of rec if <=0 on input)
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
      include 'mdef_operator_info.h'

      integer, intent(in) ::
     &     idxop, incore, rec, rec_lo, rec_hi
      logical ::
     &     standard
      type(filinf), intent(inout), target ::
     &     ffop_in
      type(operator_info), intent(inout) ::
     &     op_info

      integer ::
     &     ifree, nbuff, iocc
      type(filinf), pointer ::
     &     ffop
      type(operator), pointer ::
     &     op

      op => op_info%op_arr(idxop)%op

      if (associated(op_info%opfil_arr(idxop)%fhand))
     &     call quit(1,'assign_file_to_op',
     &     'operator is already assigned')

      if (standard) then
        allocate(op_info%opfil_arr(idxop)%fhand)
        ffop => op_info%opfil_arr(idxop)%fhand
        ! standard name
        call file_init(ffop,
     &       'op_'//trim(op%name)//'_elements.da',ftyp_da_unf,
     &       lblk_da)
      else
        op_info%opfil_arr(idxop)%fhand => ffop_in
        ffop => ffop_in
      end if

      if (rec.gt.0) then
        ffop%current_record = rec
      else
        ffop%current_record = 1
      end if

      ! set super-record length: must be multiple of primitive reclen
      ffop%length_of_record = ((op%len_op-1)/ffop%reclen+1)*ffop%reclen

      if (rec_lo.gt.0) then
        if (rec_lo.gt.rec) then
          write(luout,*) 'rec_lo, rec: ',rec_lo, rec
          call quit(1,'assign_file_to_op','record bounds inconsistent')
        end if
        ffop%active_records(1) = rec_lo
      else
        ffop%active_records(1) =
     &       ffop%current_record
      end if

      if (rec_hi.gt.0) then
        if (rec_hi.lt.rec.or.rec_hi.lt.rec_lo) then
          write(luout,*) 'rec_lo, rec, rec_hi: ',rec_lo, rec, rec_hi
          call quit(1,'assign_file_to_op','record bounds inconsistent')
        end if
        ffop%active_records(2) = rec_hi
      else
        ffop%active_records(2) =
     &       ffop%current_record
      end if

      if (incore.eq.1) then

        ffop%buffered = .true.
        ! make sure, we are in the appropriate memory section
        call mem_pushmark()
        ifree = mem_gotomark(op_files)

        ifree = mem_alloc_int(ffop%incore,
     &         op%n_occ_cls,trim(op%name)//'_incore')
        ifree = mem_alloc_int(ffop%idxrec,
     &         op%n_occ_cls,trim(op%name)//'_idxrec')
        nbuff = 0
        do iocc = 1, op%n_occ_cls
          ! allocate buffer for 1-hamiltonian only
          if (max(op%ica_occ(1,iocc),
     &              op%ica_occ(2,iocc)).gt.1 .or.
     &         iextr.gt.0.and.  ! R12: ignore currently
     &         op%ihpvca_occ(iextr,1,iocc)+
     &         op%ihpvca_occ(iextr,2,iocc).gt.0 ) then
            ffop%incore(iocc) = -1
            ffop%idxrec(iocc) = -1
          else
            ffop%idxrec(iocc) = nbuff
            ffop%incore(iocc) = op%len_op_occ(iocc)
            nbuff = nbuff + op%len_op_occ(iocc)
          end if
        end do
        ffop%nbuffer = nbuff
        ifree = mem_alloc_real(ffop%buffer,
     &       nbuff,trim(op%name)//'_buffer')

        call mem_popmark()

      else if (incore.gt.1) then
        call quit(1,'assign_file_to_op','incore>1 not supported!')
      end if

      if (iprlvl.ge.10) then
        write(luout,'(3x,4a)')
     &       'assigned file: ',trim(ffop%name),' to operator ',
     &       trim(op%name)
        write(luout,'(3x,a,i4,a,i4,a,i4,a)') 'record: ',
     &       op_info%opfil_arr(idxop)%fhand%current_record,
     &       ' (active: ',op_info%opfil_arr(idxop)%
     &       fhand%active_records(1),
     &       ' -- ',      op_info%opfil_arr(idxop)%
     &       fhand%active_records(2),
     &       ' )'
      end if

      return

      end
