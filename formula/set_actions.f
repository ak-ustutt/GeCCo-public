*----------------------------------------------------------------------*
      subroutine set_actions(act_list,nactions,
     &                       form_info,op_list,nops)
*----------------------------------------------------------------------*
*     according to the input requests, set up the necessary steps,
*     e.g. import of integrals, solve ground state equations, solve
*     further equations .... etc.
*     
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_input.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_filinf.h'
      include 'def_file_list.h'
      include 'def_action.h'
      include 'def_action_list.h'
      include 'mdef_formula_info.h'
      include 'explicit.h'

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      type(formula_info), intent(in) ::
     &     form_info
      type(operator_list), intent(in) ::
     &     op_list
      integer, intent(in) ::
     &     nops

      type(operator), pointer ::
     &     ops(:)

      integer, parameter ::
     &     ntest = 100
      integer ::
     &     iprint

      iprint = max(ntest,iprlvl)

      ! set up pointer arrays for operators and formulae
      allocate(ops(nops))
      call op_list2arr(op_list,ops,nops)

      ! Set the actions required to perform certain tasks.
      if(.not.explicit)then
        ! Normal CC GS energy calculations.
        if (is_keyword_set('method.CC').gt.0) then
          call add_cc_default_actions(act_list,nactions,
     &         ops,nops,
     &         form_info)
          
        end if
      else
        ! CC-R12 GS energy calculations
        call add_r12_default_actions(act_list,nactions,
     &       ops,nops,
     &       form_info)
      endif  

      deallocate(ops)

      if (iprint.ge.10) then
        write(luout,*) 'at end of set_actions:'
        call print_action_list(act_list)
      end if

      return
      end 
