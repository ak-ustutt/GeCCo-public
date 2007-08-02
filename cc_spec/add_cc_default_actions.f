*----------------------------------------------------------------------*
      subroutine add_cc_default_actions(act_list,nactions,
     &                                  op_info,
     &                                  form_info)
*----------------------------------------------------------------------*
*     set all actions to be taken if a CC calculation is intended
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'cc_routes.h'
      include 'ifc_input.h'
c      include 'ifc_formula.h'
c      include 'ifc_operators.h'
      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
c      include 'def_filinf.h'
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
     &     idxham, idxtop, idxdia, idxccen, idxccrs, idxomg, idxhhat,
     &     idxrhs, idxlhtr, idxtbar, idxecc, idxtbara, idxeta,
     &     idum, isim
  
      ! explicit interface does not work with ifort
      integer, external ::
     &     idx_oplist2, idx_formlist

      idxecc = idx_oplist2(op_ccen,op_info)
      idxham = idx_oplist2(op_ham,op_info)
      idxtop = idx_oplist2(op_top,op_info)
      idxomg = idx_oplist2(op_omg,op_info)
      idxdia = idx_oplist2(op_dia1,op_info)

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

      ! set up diagonal preconditioner
      call add_action(act_list,nactions,
     &     iaction_setup_prc,2,1,0,
     &     (/idxtop,idxham/),(/idxdia/),
     &     (/(/idxtop,1/),(/idxham,1/)/),(/(/idxdia,1/)/),
     &     0,idum
     &     )

c      call get_argument_value('calculate.routes','simtraf',ival=isim)

      if (ccsimtrf.eq.0) then
        ! solve ground-state equations
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,3,1,
     &     (/idxdia,idxham/),(/idxtop,idxomg,idxecc/),
     &     (/(/idxdia,1/),(/idxham,1/)/),(/(/idxtop,1/),(/idxomg,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )
      else
        idxhhat = idx_oplist2(op_hhat,op_info)
        call add_action(act_list,nactions,
     &     iaction_solve_nleq,2,4,1,
     &     (/idxdia,idxham/),(/idxtop,idxomg,idxecc,idxhhat/),
     &     (/(/idxdia,1/),(/idxham,1/)/),(/(/idxtop,1/),(/idxomg,1/)/),
     &     2,(/idxccen,idxccrs/)
     &     )
      end if

      if (solve_tbar.and.ccsimtrf.eq.0) then
        idxtbar = idx_oplist2(op_tbar,op_info)
        idxtbara = idx_oplist2(op_tbar_a,op_info)
        idxeta =  idx_oplist2(op_eta,op_info)
        idxlhtr = idx_formlist(label_cctbar_a,form_info)
        idxrhs = idx_formlist(label_cceta,form_info)
        call add_action(act_list,nactions,
     &       iaction_solve_leq,4,3,1,
     &       (/idxdia,idxham,idxtop,idxomg/),
     &                                   (/idxtbar,idxtbara,idxeta/),
     &       (/(/idxdia,1/),(/idxham,1/),(/idxtop,1/),(/idxomg,1/)/),
     &       (/(/idxtbar,1/)/),
     &       2,(/idxrhs,idxlhtr/))
      else if (solve_tbar.and.ccsimtrf.eq.1) then
        idxhhat = idx_oplist2(op_hhat,op_info)
        idxtbar = idx_oplist2(op_tbar,op_info)
        idxtbara = idx_oplist2(op_tbar_a,op_info)
        idxeta =  idx_oplist2(op_eta,op_info)
        idxlhtr = idx_formlist(label_cctbar_a,form_info)
        idxrhs = idx_formlist(label_cceta,form_info)
        call add_action(act_list,nactions,
     &       iaction_solve_leq,5,3,1,
     &       (/idxdia,idxham,idxhhat,idxtop,idxomg/),
     &                               (/idxtbar,idxtbara,idxeta/),
     &       (/(/idxdia,1/),(/idxham,1/),(/idxhhat,1/),
     &                      (/idxtop,1/),(/idxomg,1/)/),
     &       (/(/idxtbar,1/)/),
     &       2,(/idxrhs,idxlhtr/))

      end if

      return
      end
