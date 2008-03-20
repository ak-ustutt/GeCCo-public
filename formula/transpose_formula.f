*----------------------------------------------------------------------*
      subroutine transpose_formula(form_head,op_info)
*----------------------------------------------------------------------*
*     transpose all entries on a formula list
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_ptr

      form_ptr => form_head
      do
        if (form_ptr%command.eq.command_add_contribution)
     &       call transpose_contr(form_ptr%contr,op_info)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

