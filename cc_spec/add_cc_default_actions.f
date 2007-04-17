*----------------------------------------------------------------------*
      subroutine add_cc_default_actions(act_list,nactions,
     &                                  ops,nops,
     &                                  fform,nform)
*----------------------------------------------------------------------*
*     set all actions to be taken if a CC calculation is intended
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'ifc_formula.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      integer, intent(in) ::
     &     nform, nops
      type(filinf), intent(in) ::
     &     fform(nform)
      type(operator), intent(in) ::
     &     ops(nops)

      integer ::
     &     idxham, idxtop, idxccen, idxccrs
      

      ! preliminary:
      idxham = 1
      idxtop = 2

      idxccen = 2
      idxccrs = 3
      
      ! import Hamiltonian
      call add_action(act_list,nactions,
     &     iaction_import,nop_out=1,
     &     idxopdef_out=(/idxham/),
     &     idxopfile_out=(/(/idxham,1/)/)
     &     )

      ! solve ground-state equations
      call add_action(act_list,nactions,
     &     iaction_solve_nleq,nop_in=1,nop_out=1,
     &     idxopdef_in=(/idxham/),
     &     idxopfile_in=(/(/idxham,1/)/),
     &     idxopdef_out=(/idxtop/),
     &     idxopfile_out=(/(/idxtop,1/)/),
     &     nform=2,
     &     idx_formula=(/idxccen,idxccrs/)
     &     )

      return
      end
