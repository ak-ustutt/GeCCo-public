*----------------------------------------------------------------------*
      subroutine init_op_files(op_info)
*----------------------------------------------------------------------*
*     basic initialization of files containing operator matrix elements
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'par_globalmarks.h'
      include 'par_opnames_gen.h'
      include 'ifc_memman.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout) ::
     &     op_info

      type(operator_array), pointer ::
     &     ops(:)
      type(operator), pointer ::
     &     cur_op
      type(file_list), pointer ::
     &     file_ptr
      type(filinf), pointer ::
     &     cur_file

      integer ::
     &     iop, ifree, nbuff, iocc, nops

      nops = op_info%nops

      if (    .not.associated(op_info%op_arr)
     &    .or..not.associated(op_info%opfil_list))
     &     call quit(1,'init_op_files',
     &     'op_info not initialized correctly!')

      file_ptr => op_info%opfil_list
      ops => op_info%op_arr

      do iop = 1, nops
        allocate(file_ptr%fhand)
        cur_file => file_ptr%fhand
        cur_op   => ops(iop)%op
        ! basic initialization
        call file_init(cur_file,
     &       'op_'//trim(cur_op%name)//'_elements.da',ftyp_da_unf,
     &       lblk_da)
        ! incore options
        ! preliminary version:
        if (trim(cur_op%name).eq.op_ham.or.cur_op%len_op.eq.1) then
c        if (trim(cur_op%name).eq.op_ham) then
          cur_file%buffered = .true.
          ! make sure, we are in the appropriate memory section
          call mem_pushmark()
          ifree = mem_gotomark(op_files)

          ifree = mem_alloc_int(cur_file%incore,
     &         cur_op%n_occ_cls,trim(cur_op%name)//'_incore')
          ifree = mem_alloc_int(cur_file%idxrec,
     &         cur_op%n_occ_cls,trim(cur_op%name)//'_idxrec')
          nbuff = 0
          do iocc = 1, cur_op%n_occ_cls
            ! allocate buffer for 1-hamiltonian only
            if (max(cur_op%ica_occ(1,iocc),
     &              cur_op%ica_occ(2,iocc)).gt.1 .or.
     &          iextr.gt.0.and. ! R12: ignore currently
     &          cur_op%ihpvca_occ(iextr,1,iocc)+
     &          cur_op%ihpvca_occ(iextr,2,iocc).gt.0 ) then
              cur_file%incore(iocc) = -1
              cur_file%idxrec(iocc) = -1
            else
              cur_file%idxrec(iocc) = nbuff
              cur_file%incore(iocc) = cur_op%len_op_occ(iocc)
              nbuff = nbuff + cur_op%len_op_occ(iocc)
            end if
          end do
          cur_file%nbuffer = nbuff
          ifree = mem_alloc_real(cur_file%buffer,
     &         nbuff,trim(cur_op%name)//'_buffer')

          call mem_popmark()

        end if
        ! advance list pointer
        if (iop.lt.nops) then
          allocate(file_ptr%next)
          file_ptr%next%prev => file_ptr
          file_ptr => file_ptr%next
        end if
      end do

      ! update direct access array
      allocate(op_info%opfil_arr(nops))
      call file_list2arr2(op_info%opfil_list,op_info%opfil_arr,nops)

      return
      end
