*----------------------------------------------------------------------*
      subroutine print_form_item(lulog,idx,fl_item,op_info)
*----------------------------------------------------------------------*
*     print formula item to unit lulog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog
      integer, intent(inout) ::
     &     idx
      type(formula_item), intent(in), target ::
     &     fl_item
      type(operator_info), intent(in) ::
     &     op_info

      ! just a convenience wrapper
      call print_form_item2(lulog,'long',idx,fl_item,op_info)
      
      end
