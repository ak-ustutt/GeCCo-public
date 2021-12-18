*----------------------------------------------------------------------*
      subroutine op_list2arr(op_list,op_arr,nops)
*----------------------------------------------------------------------*
*     set up an array of pointers that point to the entries of the
*     list on op_list
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_operator_list.h'

      type(operator_list), intent(in), target ::
     &     op_list
      integer, intent(in) ::
     &     nops
      type(operator), intent(out) ::
     &     op_arr(nops)

      type(operator_list), pointer ::
     &     current

      integer ::
     &     iop

c      if (.not.associated(op_list))
c     &     call quit(1,'op_list2arr','operator list not initialized')

      current => op_list
      do iop = 1, nops
        if (.not.associated(current%op))
     &       call quit(1,'op_list2arr','unallocated operator on list')
        op_arr(iop) = current%op
        if (iop.lt.nops.and..not.(associated(current%next)))
     &       call quit(1,'op_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
