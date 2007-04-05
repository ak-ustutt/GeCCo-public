*----------------------------------------------------------------------*
      subroutine set_actions(act_list,nactions,
     &                       form_list,nform,op_list,nops)
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

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      type(file_list), intent(in) ::
     &     form_list
      integer, intent(inout) ::
     &     nform
      type(operator_list), intent(in) ::
     &     op_list
      integer, intent(in) ::
     &     nops

      type(operator), pointer ::
     &     ops(:)
      type(filinf), pointer ::
     &     fform(:)

      ! set up pointer arrays for operators and formulae
      allocate(ops(nops))
      call op_list2arr(op_list,ops,nops)
      allocate(fform(nops))
      call file_list2arr(form_list,fform,nform)
      
      deallocate(ops)
      deallocate(fform)

      return
      end 
