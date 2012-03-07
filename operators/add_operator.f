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

      integer ::
     &     iops

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

      ! add user-supplied label
      if (len_trim(label).gt.len_opname)
     &   call quit(1,'add_operator','name too long: "'
     &     //trim(label)//'"')
      list_pnt%op%name(1:len_opname) = ' '
      list_pnt%op%name = trim(label)
      list_pnt%op%assoc_list(1:2*len_opname) = ' '

      ! init all pointers
      list_pnt%op%ihpvca_occ => null()
      list_pnt%op%ica_occ => null()
      list_pnt%op%igasca_restr => null()

      ! increment counter
      op_info%nops = op_info%nops+1

      ! update operator array
      call update_op_arr(op_info)

      return
      end
