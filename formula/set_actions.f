*----------------------------------------------------------------------*
      subroutine set_actions(act_list,nactions,
     &                       form_info,op_info)
*----------------------------------------------------------------------*
*     according to the input requests, set up the necessary steps,
*     e.g. import of integrals, solve ground state equations, solve
*     further equations .... etc.
*     
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_input.h'
c      include 'def_filinf.h'
c      include 'def_file_list.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_action.h'
      include 'def_action_list.h'
      include 'explicit.h'

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      type(formula_info), intent(in) ::
     &     form_info
      type(operator_info), intent(in) ::
     &     op_info

      integer, parameter ::
     &     ntest = 100
      integer ::
     &     iprint

      iprint = max(ntest,iprlvl)

      ! Set the actions required to perform certain tasks.
      if(.not.explicit)then
        ! Normal CC GS energy calculations.
        if (is_keyword_set('method.CC').gt.0) then
          call add_cc_default_actions(act_list,nactions,
     &                              op_info,
     &                              form_info)          
        end if
      else
        ! CC-R12 GS energy calculations
        call add_r12_default_actions(act_list,nactions,
     &       op_info,
     &       form_info)
      endif  

      if (iprint.ge.10) then
        write(luout,*) 'at end of set_actions:'
        call print_action_list(act_list)
      end if

      return
      end 
