*----------------------------------------------------------------------*
      subroutine add_r12_default_actions(act_list,nactions,
     &                                  ops,nops,
     &                                  form_info)
*----------------------------------------------------------------------*
*     set all actions to be taken if a CC-R12 calculation is intended
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
c      include 'ifc_formula.h'
      include 'def_operator.h'
c      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'
      include 'mdef_formula_info.h'
      include 'explicit.h'

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      integer, intent(in) ::
     &     nops
      type(formula_info), intent(in) ::
     &     form_info
      type(operator), intent(in) ::
     &     ops(nops)

      integer ::
     &     idxham, idxsop, idxdia, idxccen, idxccrs, idxomg, idxhhat,
     &     idum, isim, idxr12
  
      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist, idx_formlist

      idxham = idx_oplist(op_ham,ops,nops)
      idxsop = idx_oplist(op_sop,ops,nops)
      idxomg = idx_oplist(op_omgr12,ops,nops)
      idxdia = idx_oplist(op_diar12,ops,nops)
      idxr12 = idx_oplist(op_rint,ops,nops)

      ! preliminary:
      idxccen = idx_formlist(label_ccen0,form_info)
      idxccrs = idx_formlist(label_ccrs0,form_info)
      
      ! import Hamiltonian
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxham/),
     &     idum,(/(/idxham,1/)/),
     &     0,idum
     &     )

      ! import R12 operator integrals
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxr12/),
     &     idum,(/(/idxr12,1/)/),
     &     0,idum
     &     )

      ! set up diagonal preconditioner
      call add_action(act_list,nactions,
     &     iaction_setup_prc,2,1,0,
     &     (/idxsop,idxham/),(/idxdia/),
     &     (/(/idxsop,1/),(/idxham,1/)/),(/(/idxdia,1/)/),
     &     0,idum
     &     )

      call get_argument_value('calculate.routes','simtraf',ival=isim)

      if (isim.eq.0) then
        ! solve ground-state equations
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,2,1,
     &     (/idxdia,idxham/),(/idxsop,idxomg/),
     &     (/(/idxdia,1/),(/idxham,1/)/),(/(/idxsop,1/),(/idxomg,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )
      else
        idxhhat = idx_oplist(op_hhat,ops,nops)
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,3,1,
     &     (/idxdia,idxham/),(/idxsop,idxomg,idxhhat/),
     &     (/(/idxdia,1/),(/idxham,1/)/),(/(/idxsop,1/),(/idxomg,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )
      end if

      return
      end
