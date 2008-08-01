*----------------------------------------------------------------------*
      subroutine set_cc_gsrsp_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in CC ground state response calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting targets for CC properties ...'

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Entries needed for left-hand equation:
      ! RHS vector eta
      call add_target(op_eta,ttype_op,.false.,tgt_info)
      call set_dependency(op_eta,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_eta,ttype_op,CLONE_OP,
     &              op_eta,1,1,
     &              parameters,1,tgt_info)
      
      ! MV-product Tbar.A
      call add_target(op_tbar_a,ttype_op,.false.,tgt_info)
      call set_dependency(op_tbar_a,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_tbar_a,ttype_op,CLONE_OP,
     &              op_tbar_a,1,1,
     &              parameters,1,tgt_info)
      
      ! one-particle density
      call add_target(op_1dens,ttype_op,.false.,tgt_info)
      call dens_parameters(-1,parameters,
     &                     1,1,1)
      call set_rule(op_1dens,ttype_op,DEF_DENSITY,
     &              op_1dens,1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae:
*----------------------------------------------------------------------*
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cctbar_a
      labels(2) = form_cclg0
      labels(3) = op_tbar_a
      labels(4) = op_top
      labels(5) = ' '
      call add_target(form_cctbar_a,ttype_frm,.false.,tgt_info)
      call add_target(form_cceta,ttype_frm,.false.,tgt_info)
      call set_joined_targets(form_cctbar_a,form_cceta,tgt_info)
      call set_dependency(form_cctbar_a,form_cclg0,tgt_info)
      call set_dependency(form_cctbar_a,op_tbar_a,tgt_info)
      call set_dependency(form_cctbar_a,op_eta,tgt_info)
      call set_rule(form_cctbar_a,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cctbar_a,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cctbar_a
      labels(2) = form_cceta
      labels(3) = form_cctbar_a
      labels(4) = op_tbar_a
      labels(5) = op_eta
      labels(6) = op_tbar
      call set_rule(form_cctbar_a,ttype_frm,LEQ_SPLIT,
     &              labels,6,2,
     &              title_cctbar_a,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc1dens
      labels(2) = form_cclg0
      labels(3) = op_1dens
      labels(4) = op_ham
      labels(5) = ' '
      call add_target(form_cc1dens,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc1dens,form_cclg0,tgt_info)
      call set_dependency(form_cc1dens,op_tbar_a,tgt_info)
      call set_dependency(form_cc1dens,op_eta,tgt_info)
      call set_rule(form_cc1dens,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc1dens,1,tgt_info)

*----------------------------------------------------------------------*
*     Optimized Formulae:
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC ground state left-hand eq.
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cclft0
      labels(2) = form_cceta
      labels(3) = form_cctbar_a
      ncat = 2
      nint = 0
      call add_target(fopt_cclft0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cclft0,form_cctbar_a,tgt_info)
      call set_dependency(fopt_cclft0,form_cceta,tgt_info)
      call set_dependency(fopt_cclft0,mel_tbar_adef,tgt_info)
      call set_dependency(fopt_cclft0,mel_etadef,tgt_info)
      call set_dependency(fopt_cclft0,mel_tbardef,tgt_info)
      call set_dependency(fopt_cclft0,mel_topdef,tgt_info)
      call set_dependency(fopt_cclft0,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_cclft0,form_cchhat,tgt_info)
        call set_dependency(fopt_cclft0,mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cclft0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! CC ground state density
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc1dens
      labels(2) = form_cc1dens
      call add_target(fopt_cc1dens,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cc1dens,form_cc1dens,tgt_info)
      call set_dependency(fopt_cc1dens,mel_tbardef,tgt_info)
      call set_dependency(fopt_cc1dens,mel_topdef,tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule(fopt_cc1dens,ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists:
*----------------------------------------------------------------------*
      ! ETA-list definition
      call add_target(mel_etadef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_etadef,op_eta,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_eta
      labels(2) = op_eta
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_etadef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! TBAR.A-list definition
      call add_target(mel_tbar_adef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_tbar_adef,op_tbar_a,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_tbar_a
      labels(2) = op_tbar_a
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(mel_tbar_adef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! 1DEN definition
      call add_target(meldef_1dens,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_1dens,op_1dens,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_1dens
      labels(2) = op_1dens
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule(meldef_1dens,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      icnt = is_keyword_set('calculate.CC_solve_tbar')
      needed = icnt.gt.0

      ! solve TBAR (LAMBDA) equations
      call add_target(solve_cc_lhwf,ttype_gen,needed,tgt_info)
      call set_dependency(solve_cc_lhwf,mel_dia1,tgt_info)
      call set_dependency(solve_cc_lhwf,fopt_cclft0,tgt_info)
      call set_dependency(solve_cc_lhwf,solve_cc_gs,tgt_info)
      call solve_parameters(-1,parameters,2, 1,1,'DIA')
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_tbar
      labels(2) = mel_dia1
      labels(3) = op_tbar_a
      labels(4) = op_eta
      labels(5) = fopt_cclft0
      call set_rule(solve_cc_lhwf,ttype_opme,SOLVELEQ,
     &     labels,5,1,
     &     parameters,2,tgt_info)

      call add_target(eval_1dens,ttype_gen,.false.,tgt_info)
      call set_dependency(eval_1dens,meldef_1dens,tgt_info)
      call set_dependency(eval_1dens,fopt_cc1dens,tgt_info)
      call set_dependency(eval_1dens,solve_cc_gs,tgt_info)
      call set_dependency(eval_1dens,solve_cc_lhwf,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc1dens
      call set_rule(eval_1dens,ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      icnt = is_keyword_set('calculate.properties')
      needed = icnt.gt.0

      call add_target(eval_props,ttype_gen,needed,tgt_info)
      call set_dependency(eval_props,eval_1dens,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_1dens
      call evalprop_parameters(-1,parameters,1,1,'DALTON')
      call set_rule(eval_props,ttype_opme,EVALPROP,
     &     labels,1,0,
     &     parameters,1,tgt_info)

      return
      end
