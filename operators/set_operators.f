*----------------------------------------------------------------------*
      subroutine set_operators(op_list,nops,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up operator info structures
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ifc_input.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'

      type(operator_list), intent(inout) ::
     &     op_list
      integer, intent(out) ::
     &     nops
      type(orbinf) ::
     &     orb_info   
      logical ::
     &     explicit

      explicit=.false.
      nops=0

      if (is_keyword_set('method.CC').gt.0) then
        call set_cc_operators(op_list,nops,orb_info)
      endif
      if(is_keyword_set('method.R12').gt.0)then
        explicit=.true.
        call set_r12_operators(op_list,nops,orb_info)
      end if

      ! further operator definitions may follow

      return

      end
