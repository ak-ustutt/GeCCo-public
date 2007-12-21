*----------------------------------------------------------------------*
      subroutine del_operator(name,op_info)
*----------------------------------------------------------------------*
*     delete info for operator with label "name"
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      character(*), intent(in) ::
     &     name

      type(operator_list), pointer ::
     &     list_pnt
      type(operator), pointer ::
     &     op

      list_pnt => op_info%op_list
      ! advance to respective entry on operator list:
      do while (trim(list_pnt%op%name).ne.trim(name)
     &     .and.associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      if (trim(list_pnt%op%name).ne.trim(name)) then
        call quit(1,'del_operator','unknown label: "'//trim(name)//'"')
      end if

      if (associated(list_pnt%prev)) list_pnt%prev%next => list_pnt%next
      if (associated(list_pnt%next)) list_pnt%next%prev => list_pnt%prev
      op => list_pnt%op
      call dealloc_operator(op)

      ! decrement counter
      op_info%nops = op_info%nops-1

      ! update operator array
      call update_op_arr(op_info)

      return
      end
