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
      include 'explicit.h'

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
     &     idum, isim, idxtop, idxr12, idxrba, idxttr,
     &     idxvint, idx_v_form, idxvbint, idx_vb_form, idxbint,
     &     idx_b_form, idx_b_symm, idx_b_symm_form, idxr12dia,
     &     idx_b_inv, idx_b_inv_form, idxecc, idxomg12, idxccrs12,
     &     idxr12sq, idx_xint, idx_x_form

      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2, idx_formlist

      idxham = idx_oplist2(op_ham,op_info)
      idxc12 = idx_oplist2(op_c12,op_info)
      idxomg = idx_oplist2(op_omg,op_info)
      idxomg12 = idx_oplist2(op_omgr12,op_info)
      idxdia = idx_oplist2(op_dia1,op_info)
      idxr12dia = idx_oplist2(op_diar12,op_info)
      idxtop = idx_oplist2(op_top,op_info)
      idxr12 = idx_oplist2(op_rint,op_info)
      idxrba = idx_oplist2(op_rinba,op_info)
      idxr12sq = idx_oplist2(op_f2,op_info)
      idxttr = idx_oplist2(op_ttr,op_info)
      idxecc = idx_oplist2(op_ccen,op_info)

      ! preliminary:
      idxccen = idx_formlist(label_ccen0,form_info)
      idxccrs = idx_formlist(label_ccrs0,form_info)
      idxccrs12 = idx_formlist(label_ccrs12,form_info)
      
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
     &     0,idum
     &     )

      if(trim(r12_apprx).ne.'A')then
        ! Import the R12**2 operator integrals.
        call add_action(act_list,nactions,
     &       iaction_import,0,1,0,
     &       idum,(/idxr12sq/),
     &       0,idum
     &       )
      endif

      ! Import commutator integrals, (pq|[T1+T2,r12]|rs).
      call add_action(act_list,nactions,
     &     iaction_import,0,1,0,
     &     idum,(/idxttr/),
     &     0,idum
     &     )
      
      ! Evaluate the V-intermediate.
      idxvint = idx_oplist2(op_v_inter,op_info)
      idx_v_form = idx_formlist(label_r12_vint,form_info)
      call add_action(act_list,nactions,
     &     iaction_evaluate,0,1,0,
     &     idum,(/idxvint/),
     &     1,(/idx_v_form/)
     &     )

      ! Evaluate the V+-intermediate.
      idxvbint = idx_oplist2(op_vbar_inter,op_info)
      idx_vb_form = idx_formlist(label_r12_vbint,form_info)
      call add_action(act_list,nactions,
     &     iaction_evaluate,0,1,0,
     &     idum,(/idxvbint/),
     &     1,(/idx_vb_form/)
     &     )

      ! Evaluate the X-intermediate, if necessary.
      if(trim(r12_apprx).ne.'A')then
        idx_xint = idx_oplist2(op_x_inter,op_info)
        idx_x_form = idx_formlist(label_r12_xint,form_info)
        call add_action(act_list,nactions,
     &       iaction_evaluate,0,1,0,
     &       idum,(/idx_xint/),
     &       1,(/idx_x_form/)
     &       )
      endif

      ! Evaluate the B-intermediate.
      idxbint = idx_oplist2(op_b_inter,op_info)
      idx_b_form = idx_formlist(label_r12_bint,form_info)
      call add_action(act_list,nactions,
     &     iaction_evaluate,0,1,0,
     &     idum,(/idxbint/),
     &     1,(/idx_b_form/)
     &     )

      ! Symmetrise the B-matrix.
      idx_b_symm = idx_oplist2(op_b_symm,op_info)
      idx_b_symm_form = idx_formlist(label_r12_bsymm,form_info)
      call add_action(act_list,nactions,
     &     iaction_symmetrise,1,1,0,
     &     (/idxbint/),(/idx_b_symm/),
     &     1,(/idx_b_symm_form/)
     &     )

      if(mp2)then
        ! Invert the B-matrix.
        idx_b_inv = idx_oplist2(op_b_inv,op_info)
        idx_b_inv_form = idx_formlist(label_r12_binv,form_info)
        call add_action(act_list,nactions,
     &       iaction_invert,1,1,0,
     &       (/idx_b_symm/),(/idx_b_inv/),
     &       1,(/idx_b_inv_form/)
     &       )

c        ! Multiply the necessary intermediates by the inverse of B.
c        call add_action(act_list,nactions,
c     &       iaction_multiply,2,1,0,
c     &       (/idx_b_inv,idxvbint/),(/idxvbint/),
c     &       1,(/idx_vb_form/)
c     &       )

c        call add_action(act_list,nactions,
c     &       iaction_multiply,2,1,0,
c     &       (/idx_b_inv,idx_b_symm/),(/idx_b_symm/),
c     &       1,(/idx_b_symm_form/)
c     &       )

c        if(trim(r12_apprx).ne.'A')then
c          call add_action(act_list,nactions,
c     &         iaction_multiply,2,1,0,
c     &         (/idx_b_inv,idx_xint/),(/idx_xint/),
c     &         1,(/idx_x_form/)
c     &         )
c        endif

      endif

c dbg
c      ! Evaluate the energy of the 1st iteration of the MP2-R12 equations.
c      call add_action(act_list,nactions,
c     &     iaction_evaluate,0,1,0,
c     &     idum,(/idxecc/),
c     &     1,(/idxccen/),
c     &     )
c dbg

      ! set up standard diagonal preconditioner
      call add_action(act_list,nactions,
     &     iaction_setup_prc,2,1,0,
     &     (/idxtop,idxham/),(/idxdia/),
     &     0,idum
     &     )

c      ! Set up the R12 preconditioner (Diagonal of B).
c      call add_action(act_list,nactions,
c     &     iaction_setup_prc,2,1,0,
c     &     (/idxbint,idxbint/),(/idxr12dia/),
c     &     0,idum
c     &     )

c      call get_argument_value('calculate.routes','simtraf',ival=isim)

c        call add_action(act_list,nactions,
c     &     iaction_solve_nleq,3,5,2,
c     &     (/idxdia,idxr12dia,idxham/),
c     &     (/idxtop,idxc12,idxomg,idxomg12,idxecc/),
c     &     3,(/idxccen,idxccrs,idxccrs12/)
c     &     )

c dbg
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,3,5,2,
     &     (/idxdia,idx_b_inv,idxham/),
     &     (/idxtop,idxc12,idxomg,idxomg12,idxecc/),
     &     3,(/idxccen,idxccrs,idxccrs12/)
     &     )
c dbg

c      if (isim.eq.0) then
c        ! solve ground-state equations
c        call add_action(act_list,nactions,
c     &     iaction_solve_nleq,2,2,1,
c     &     (/idxdia,idxham/),(/idxsop,idxomg/),
c     &     2,(/idxccen,idxccrs/)
c     &     )
c      else
c        idxhhat = idx_oplist2(op_hhat,op_info)
c        call add_action(act_list,nactions,
c     &     iaction_solve_nleq,2,3,1,
c     &     (/idxdia,idxham/),(/idxsop,idxomg,idxhhat/),
c     &     2,(/idxccen,idxccrs/)
c     &     )
c      end if

      return
      end
