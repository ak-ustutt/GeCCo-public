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
c      include 'ifc_formula.h'
      include 'def_operator.h'
c      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
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
     &     idxham, idxtop, idxdia, idxccen, idxccrs, idum
  
      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist

      idxham = idx_oplist(op_ham,ops,nops)
      idxtop = idx_oplist(op_top,ops,nops)
      idxdia = idx_oplist(op_dia1,ops,nops)

      ! preliminary:
      idxccen = 2
      idxccrs = 3
      
      ! import Hamiltonian
      call add_action(act_list,nactions,
     &     iaction_import,0,1,
     &     idum,(/idxham/),
     &     idum,(/(/idxham,1/)/),
     &     0,idum
     &     )

      ! set up diagonal preconditioner
      call add_action(act_list,nactions,
     &     iaction_setup_prc,2,1,
     &     (/idxtop,idxham/),(/idxdia/),
     &     (/(/idxdia,1/)/),(/(/idxtop,1/),(/idxham,1/)/),
     &     0,idum
     &     )

      ! solve ground-state equations
      call add_action(act_list,nactions,
     &     iaction_solve_nleq,1,1,
     &     (/idxham/),(/idxtop/),
     &     (/(/idxham,1/)/),(/(/idxtop,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )

      return
      end
