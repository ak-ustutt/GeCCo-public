*----------------------------------------------------------------------*
      subroutine add_r12_default_actions(act_list,nactions,
     &                                  op_info,
     &                                  form_info)
*----------------------------------------------------------------------*
*     set all actions to be taken if a CC-R12 calculation is intended
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'mdef_operator_info.h'
      include 'mdef_formula_info.h'
      include 'def_action.h'
      include 'def_action_list.h'

      type(action_list), intent(inout) ::
     &     act_list
      integer, intent(out) ::
     &     nactions
      type(formula_info), intent(in) ::
     &     form_info
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idxham, idxsop, idxdia, idxccen, idxccrs, idxomg, idxhhat,
     &     idum, isim, idxr12
  
      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2, idx_formlist

      idxham = idx_oplist2(op_ham,op_info)
      idxsop = idx_oplist2(op_sop,op_info)
      idxomg = idx_oplist2(op_omgr12,op_info)
      idxdia = idx_oplist2(op_diar12,op_info)
      idxr12 = idx_oplist2(op_rint,op_info)

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
        idxhhat = idx_oplist2(op_hhat,op_info)
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,3,1,
     &     (/idxdia,idxham/),(/idxsop,idxomg,idxhhat/),
     &     (/(/idxdia,1/),(/idxham,1/)/),(/(/idxsop,1/),(/idxomg,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )
      end if

      return
      end
