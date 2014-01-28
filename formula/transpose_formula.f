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

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      type(formula_item), pointer ::
     &     form_ptr

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'transpose_formula')
        write(lulog,*) 'formula on entry'
        call print_form_list(lulog,form_head,op_info)
      end if

      form_ptr => form_head
      do
        form_ptr%target = -form_ptr%target
        if (form_ptr%command.eq.command_add_contribution)
     &       call transpose_contr(form_ptr%contr,op_info)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      if (ntest.ge.100) then
        write(lulog,*) 'formula on exit'
        call print_form_list(lulog,form_head,op_info)
      end if

      return
      end

