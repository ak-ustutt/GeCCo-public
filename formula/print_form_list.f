*----------------------------------------------------------------------*
      subroutine print_form_list(luout,form_head,op_info)
*----------------------------------------------------------------------*
*     print formula on linked list to unit luout
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     luout
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx

      type(formula_item), pointer ::
     &     form_ptr

      idx = 0
      form_ptr => form_head
      do
        call print_form_item(luout,idx,form_ptr,op_info)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

