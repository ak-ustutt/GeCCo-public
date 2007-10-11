*----------------------------------------------------------------------*
      subroutine add_operator(label,op_info)
*----------------------------------------------------------------------*
*     allocate a new slot for an operator structure in op_info
*     the list is extended and the array is updated
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      character ::
     &     label*(*)

      type(operator_list), pointer ::
     &     list_pnt
      type(file_array), pointer ::
     &     new_opfil_arr(:)

      integer ::
     &     iops

c dbg
      type(operator), pointer :: apnt,bpnt
c dbg

      list_pnt => op_info%op_list
      ! advance to end of operator list:
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%op)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        list_pnt%next => null()
      end if
      allocate (list_pnt%op)

      ! assign unique ID
      op_info%id_cnt = op_info%id_cnt+1
      list_pnt%op%id = op_info%id_cnt

      ! add user-supplied label
      if (len_trim(label).gt.len_opname)
     &   call quit(1,'add_operator','name too long: "'
     &     //trim(label)//'"')
      list_pnt%op%name = '        '
      list_pnt%op%name = trim(label)

      ! init all pointers
      list_pnt%op%ihpvca_occ => null()
      list_pnt%op%ica_occ => null()
      list_pnt%op%igasca_restr => null()
      list_pnt%op%len_op_occ => null()
      list_pnt%op%off_op_occ => null()
      list_pnt%op%len_op_gmo => null()
      list_pnt%op%off_op_gmo => null()
      list_pnt%op%len_op_gmox => null()
      list_pnt%op%off_op_gmox => null()
      list_pnt%op%idx_graph => null()

      ! increment counter
      op_info%nops = op_info%nops+1

c dbg
      if (op_info%nops-1.gt.0)
     &     apnt => op_info%op_arr(op_info%nops-1)%op
c dbg
      ! update operator array
      call update_op_arr(op_info)
c dbg
      if (op_info%nops-1.gt.0) then
        bpnt => op_info%op_arr(op_info%nops-1)%op
        print *,'check in add_operator: ',
     &       associated(apnt,bpnt)
        print *,op_info%op_arr(op_info%nops-1)%op%name
        print *,op_info%op_arr(op_info%nops)%op%name
      end if
c dbg

      ! update operator file array
      allocate(new_opfil_arr(op_info%nops))
      
      do iops = 1, op_info%nops-1
        new_opfil_arr(iops)%fhand => op_info%opfil_arr(iops)%fhand
      end do
      new_opfil_arr(op_info%nops)%fhand => null()

      if (associated(op_info%opfil_arr)) deallocate(op_info%opfil_arr)
      op_info%opfil_arr => new_opfil_arr

c dbg
      if (op_info%nops-1.gt.0) then
        print *,op_info%op_arr(op_info%nops-1)%op%name
        print *,op_info%op_arr(op_info%nops)%op%name
      end if
c dbg

      return
      end
