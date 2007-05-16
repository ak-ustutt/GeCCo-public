*----------------------------------------------------------------------*
      subroutine set_formulae(form_info,op_list,nops)
*----------------------------------------------------------------------*
*     driver routine for setting up formula files
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_input.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_filinf.h'
      include 'mdef_formula_info.h'

      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_list), intent(in) ::
     &     op_list
      integer, intent(in) ::
     &     nops

      type(operator), pointer ::
     &     ops(:)

      ! set up pointer array for operators, which is more
      ! convenient than a linked list
      allocate(ops(nops))
      call op_list2arr(op_list,ops,nops)

      if (is_keyword_set('method.CC').gt.0) then
        call set_cc_formula(form_info,ops,nops)
      end if

      deallocate(ops)

      return
      end
