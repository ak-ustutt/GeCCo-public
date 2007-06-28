*----------------------------------------------------------------------*
      subroutine add_operator(op_info)
*----------------------------------------------------------------------*
*     allocate a new slot for an operator structure in op_info
*     the list is extended and the array is updated
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info

      type(operator_list), pointer ::
     &     list_pnt

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

      ! increment counter
      op_info%nops = op_info%nops+1

      ! update operator array
      call update_op_arr(op_info)

      return
      end
