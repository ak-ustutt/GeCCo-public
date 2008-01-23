*----------------------------------------------------------------------*
      subroutine set_formulae(form_info,op_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up formula files
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_input.h'
      include 'mdef_operator_info.h'
c      include 'def_operator_list.h'
c      include 'def_filinf.h'
      include 'mdef_formula_info.h'
      include 'def_orbinf.h'
      include 'cc_routes.h'

      type(formula_info), intent(inout) ::
     &     form_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf) ::
     &     orb_info

      type(operator), pointer ::
     &     ops(:)


      if (do_cc.or.do_mp) then
c        ! OLD
c        ! set up pointer array for operators, which is more
c        ! convenient than a linked list (DIFFERS FROM NEW op_arr)
c        allocate(ops(op_info%nops))
c        call op_list2arr(op_info%op_list,ops,op_info%nops)
c        ! OLD
c
c        call set_cc_formula(form_info,ops,op_info%nops)
c
c        ! OLD
c        deallocate(ops)
c        ! OLD

        call set_cc_formula2(form_info,op_info,orb_info)

      end if

      return
      end
