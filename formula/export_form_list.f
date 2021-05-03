*----------------------------------------------------------------------*
      subroutine export_form_list(lulog,form_head,op_info)
*----------------------------------------------------------------------*
*     export formula on linked list to unit lulog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog
      type(formula_item), intent(inout), target ::
     &     form_head
      type(operator_info), intent(inout) ::
     &     op_info

      integer ::
     &     idx

      type(formula_item), pointer ::
     &     form_ptr

      idx = 0
      form_ptr => form_head
      do
        call export_form_item(lulog,idx,form_ptr,op_info)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

