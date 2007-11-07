*----------------------------------------------------------------------*
      subroutine add_r12_default_actions(act_list,nactions,
     &                                  op_info,
     &                                  form_info)
*----------------------------------------------------------------------*
*     Set all actions to be taken if a CC-R12 calculation is intended.
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
     &     idxham, idxc12, idxdia, idxccen, idxccrs, idxomg, idxhhat,
     &     idum, isim, idxr12, idxrba, idxsop, idxdel, idxttr,
     &     idxvint, idx_v_form, idxbint, idx_b_form, idx_b_symm,
     &     idx_b_symm_form
  
      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2, idx_formlist

      idxham = idx_oplist2(op_ham,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)
      idxomg = idx_oplist2(op_omgr12,op_info)
      idxdia = idx_oplist2(op_diar12,op_info)
      idxr12 = idx_oplist2(op_rint,op_info)
      idxrba = idx_oplist2(op_rinba,op_info)
      idxttr = idx_oplist2(op_ttr,op_info)
      idxsop = idx_oplist2(op_sop,op_info)
c      idxdel = idx_oplist2(op_del_inter,op_info)

      ! preliminary:
      idxccen = idx_formlist(label_ccen0,form_info)
      idxccrs = idx_formlist(label_ccrs0,form_info)
      
      ! import Hamiltonian
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxham/),
     &     0,idum
     &     )

      ! import R12 operator integrals.
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxr12/),
     &     0,idum
     &     )

      ! import the adjoints of the R12 operator integrals.
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxrba/),
     &     idum,(/(/idxrba,1/)/),
     &     0,idum
     &     )

      ! Import commutator integrals, (pq|[T1+T2,r12]|rs).
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxttr/),
     &     idum,(/(/idxttr,1/)/),
     &     0,idum
     &     )
      
      ! Evaluate the V-intermediate.
      idxvint = idx_oplist2(op_v_inter,op_info)
      idx_v_form = idx_formlist(label_r12_vint,form_info)
      call add_action(act_list,nactions,
     &     iaction_evaluate,0,1,0,
     &     idum,(/idxvint/),
     &     idum,(/(/idxvint,1/)/),
     &     1,(/idx_v_form/)
     &     )

      ! Evaluate the B-intermediate.
      idxbint = idx_oplist2(op_b_inter,op_info)
      idx_b_form = idx_formlist(label_r12_bint,form_info)
      call add_action(act_list,nactions,
     &     iaction_evaluate,0,1,0,
     &     idum,(/idxbint/),
     &     idum,(/(/idxbint,1/)/),
     &     1,(/idx_b_form/)
     &     )

      ! Symmetrise the B-matrix.
      idx_b_symm = idx_oplist2(op_b_symm,op_info)
      idx_b_symm_form = idx_formlist(label_r12_bsymm,form_info)
      call add_action(act_list,nactions,
     &     iaction_symmetrise,1,1,0,
     &     (/idxbint/),(/idx_b_symm/),
     &     (/(/idxbint,1/)/),(/(/idx_b_symm,1/)/),
     &     1,(/idx_b_symm_form/)
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
     &     2,(/idxccen,idxccrs/)
     &     )
      else
        idxhhat = idx_oplist2(op_hhat,op_info)
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,3,1,
     &     (/idxdia,idxham/),(/idxsop,idxomg,idxhhat/),
     &     2,(/idxccen,idxccrs/)
     &     )
      end if

      return
      end
