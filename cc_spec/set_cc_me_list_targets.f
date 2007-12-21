*----------------------------------------------------------------------*
      subroutine set_cc_me_list_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set all formula targets needed in CC calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'
      include 'def_orbinf.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     icnt, isym, ms, msc,
     &     sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting list targets for CC ...'

      ! Hamilton list:
      call add_target(mel_ham,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ham,op_ham,tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ham
      labels(2) = op_ham
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ham,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ham
      call import_parameters(-1,parameters,'DALTON') ! preliminary
      call set_rule(mel_ham,ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! L0/E0:
      call add_target(mel_cclg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_cclg0,op_cclg,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_cclg0
      labels(2) = op_cclg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_cclg0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_ccen0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccen0def,op_ccen,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ccen0
      labels(2) = op_ccen
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_ccen0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! T-list definition
      call add_target(mel_topdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_topdef,op_top,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = op_top
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_topdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! OMG-list definition
      call add_target(mel_omgdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_omgdef,op_omg,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_omg
      labels(2) = op_omg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_omgdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! TBAR-list definition
      call add_target(mel_tbardef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_tbardef,op_tbar,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_tbar
      labels(2) = op_tbar
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_tbardef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ETA-list definition
      call add_target(mel_etadef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_etadef,op_eta,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_eta
      labels(2) = op_eta
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
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
     &     0,0,1,0,0)
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
     &     0,0,1,0,0)
      call set_rule(meldef_1dens,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! HHat definition
      call add_target(mel_hhatdef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_hhatdef,op_hhat,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_hhat
      labels(2) = op_hhat
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_hhatdef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! DIA-list for all symmetries:
      do isym = 1, orb_info%nsym
        call me_list_label(me_label,mel_dia,isym,0,0,0,.false.)
        call add_target(me_label,ttype_opme,.false.,tgt_info)
        call set_dependency(me_label,mel_ham,tgt_info)
        call set_dependency(me_label,op_dia,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = me_label
        labels(2) = op_dia
        call me_list_parameters(-1,parameters,
     &       0,0,isym,0,0)
        call set_rule(me_label,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        labels(1) = me_label
        labels(2) = mel_ham
        call set_rule(me_label,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)
      end do
      
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      icnt = is_keyword_set('calculate.excitation')
      if (icnt.gt.0) then
        call get_argument_value('calculate.excitation','sym',
     &       iarr=sym_arr)
        call add_target(meldef_rex,ttype_opme,.false.,tgt_info)
        call add_target(meldef_a_rex,ttype_opme,.false.,tgt_info)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ! RE0
          call me_list_label(me_label,mel_rex,isym,0,0,0,.false.)
          call set_dependency(meldef_rex,op_r,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_r
          call me_list_parameters(-1,parameters,
     &         0,0,isym,0,0)
          call set_rule(meldef_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! A.RE0
          call me_list_label(me_label,mel_a_rex,isym,0,0,0,.false.)
          call set_dependency(meldef_a_rex,op_a_r,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_a_r
          call me_list_parameters(-1,parameters,
     &         0,0,isym,0,0)
          call set_rule(meldef_a_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        end do
      end if

      icnt = is_keyword_set('calculate.ionization')
      if (icnt.gt.0) then
        call get_argument_value('calculate.ionization','sym',
     &       iarr=sym_arr)
        call add_target(meldef_rip,ttype_opme,.false.,tgt_info)
        call add_target(meldef_a_rip,ttype_opme,.false.,tgt_info)
        call set_dependency(meldef_rip,op_rip,tgt_info)
        call set_dependency(meldef_a_rip,op_a_rip,tgt_info)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ms = +1
          msc = 0
          ! RI0
          call me_list_label(me_label,mel_rip,isym,0,ms,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_rip
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,ms)
          call set_rule(meldef_rip,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! A.RI0
          call me_list_label(me_label,mel_a_rip,isym,0,ms,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_a_rip
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,ms)
          call set_rule(meldef_a_rip,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        end do
      end if

      ! equation solve targets follow:

      call add_target(solve_cc_gs,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_cc_gs,mel_dia1,tgt_info)
      call set_dependency(solve_cc_gs,label_ccrs0opt,tgt_info)
      call solve_parameters(-1,parameters,1,1)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = mel_omg
      labels(3) = mel_dia1
      labels(4) = mel_ccen0
      labels(5) = label_ccrs0opt
      call set_rule(solve_cc_gs,ttype_opme,SOLVENLEQ,
     &     labels,5,2,
     &     parameters,1,tgt_info)

      icnt = is_keyword_set('calculate.CC_solve_tbar')
      needed = icnt.gt.0

      call add_target(solve_cc_lhwf,ttype_gen,needed,tgt_info)
      call set_dependency(solve_cc_lhwf,mel_dia1,tgt_info)
      call set_dependency(solve_cc_lhwf,label_cclft0opt,tgt_info)
      call set_dependency(solve_cc_lhwf,solve_cc_gs,tgt_info)
      call solve_parameters(-1,parameters,1,1)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_tbar
      labels(2) = mel_dia1
      labels(3) = op_tbar_a
      labels(4) = op_eta
      labels(5) = label_cclft0opt
      call set_rule(solve_cc_lhwf,ttype_opme,SOLVELEQ,
     &     labels,5,1,
     &     parameters,1,tgt_info)

      call add_target(eval_1dens,ttype_gen,.false.,tgt_info)
      call set_dependency(eval_1dens,meldef_1dens,tgt_info)
      call set_dependency(eval_1dens,label_cc1dens_opt,tgt_info)
      call set_dependency(eval_1dens,solve_cc_gs,tgt_info)
      call set_dependency(eval_1dens,solve_cc_lhwf,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc1dens_opt
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

      icnt = is_keyword_set('calculate.excitation')
      if (icnt.gt.0) then
        call get_argument_value('calculate.excitation','sym',
     &       iarr=sym_arr)
        call add_target(solve_cc_rhex,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_cc_rhex,solve_cc_gs,tgt_info)
        call set_dependency(solve_cc_rhex,label_cc_a_r_opt,tgt_info)
        call set_dependency(solve_cc_rhex,meldef_rex,tgt_info)
        call set_dependency(solve_cc_rhex,meldef_a_rex,tgt_info)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          
          call me_list_label(me_label,mel_rex,isym,0,0,0,.false.)
          call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)
          call set_dependency(solve_cc_rhex,dia_label,tgt_info)
          call solve_parameters(-1,parameters,1,sym_arr(isym))
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = dia_label
          labels(3) = op_a_r
          labels(4) = label_cc_a_r_opt
          call set_rule(solve_cc_rhex,ttype_opme,SOLVEEVP,
     &         labels,4,1,
     &         parameters,1,tgt_info)
        end do
          
      end if

      icnt = is_keyword_set('calculate.ionization')
      if (icnt.gt.0) then
        call get_argument_value('calculate.ionization','sym',
     &       iarr=sym_arr)
        call add_target(solve_cc_rhip,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_cc_rhip,solve_cc_gs,tgt_info)
        call set_dependency(solve_cc_rhip,label_cc_a_rip_opt,tgt_info)
        call set_dependency(solve_cc_rhip,op_dia_ip,tgt_info)
        call set_dependency(solve_cc_rhip,meldef_rip,tgt_info)
        call set_dependency(solve_cc_rhip,meldef_a_rip,tgt_info)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ms = +1
          msc = 0
          call me_list_label(me_label,mel_rip,isym,0,ms,msc,.false.)
          call me_list_label(dia_label,mel_dia_ip,isym,0,ms,msc,.false.)
          ! a) make diagonal
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = dia_label
          labels(2) = op_dia_ip
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,ms)
          call set_rule(solve_cc_rhip,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = dia_label
          labels(2) = mel_ham
          call set_rule(solve_cc_rhip,ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)
          ! b) solve the eigenvalue equation
          call solve_parameters(-1,parameters,1,sym_arr(isym))
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = dia_label
          labels(3) = op_a_rip
          labels(4) = label_cc_a_rip_opt
          call set_rule(solve_cc_rhip,ttype_opme,SOLVEEVP,
     &         labels,4,1,
     &         parameters,1,tgt_info)
          ! c) remove the diagonal list
          labels(1) = dia_label
          call set_rule(solve_cc_rhip,ttype_opme,DELETE_ME_LIST,
     &         labels,1,1,
     &         parameters,0,tgt_info)
        end do
          
      end if

      return
      end

