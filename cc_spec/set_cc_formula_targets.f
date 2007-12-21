*----------------------------------------------------------------------*
      subroutine set_cc_formula_targets(tgt_info)
*----------------------------------------------------------------------*
*     set all formula targets needed in CC calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info

      integer ::
     &     isim, ncat, nint, icnt
      character(len_target_name) ::
     &     labels(10)
      character(len_command_par) ::
     &     parameters

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting formula targets for CC ...'

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cclg0
      labels(2) = op_cclg
      labels(3) = op_ham
      labels(4) = op_tbar
      labels(5) = op_top
      call add_target(label_cclg0,ttype_frm,.true.,tgt_info)
      call set_dependency(label_cclg0,op_cclg,tgt_info)
      call set_dependency(label_cclg0,op_ham,tgt_info)
      call set_dependency(label_cclg0,op_tbar,tgt_info)
      call set_dependency(label_cclg0,op_top,tgt_info)
      call set_rule(label_cclg0,ttype_frm,DEF_CC_LAGRANGIAN,
     &              labels,5,1,
     &              title_cclg0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cchhat
      labels(2) = op_hhat
      labels(3) = op_ham
      labels(4) = op_top
      call add_target(label_cchhat,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cchhat,op_hhat,tgt_info)
      call set_dependency(label_cchhat,op_ham,tgt_info)
      call set_dependency(label_cchhat,op_top,tgt_info)
      call set_rule(label_cchhat,ttype_frm,DEF_HHAT,
     &              labels,4,1,
     &              title_cchhat,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_ccen0
      labels(2) = label_cclg0
      labels(3) = op_ccen
      labels(4) = op_tbar
      call add_target(label_ccen0,ttype_frm,.false.,tgt_info)
      call set_dependency(label_ccen0,label_cclg0,tgt_info)
      call set_dependency(label_ccen0,op_ccen,tgt_info)
      call set_rule(label_ccen0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ccen0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_ccrs0
      labels(2) = label_cclg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(label_ccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(label_ccrs0,label_cclg0,tgt_info)
      call set_dependency(label_ccrs0,op_omg,tgt_info)
      call set_rule(label_ccrs0,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccrs0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cctbar_a
      labels(2) = label_cclg0
      labels(3) = op_tbar_a
      labels(4) = op_top
      labels(5) = ' '
      call add_target(label_cctbar_a,ttype_frm,.false.,tgt_info)
      call add_target(label_cceta,ttype_frm,.false.,tgt_info)
      call set_joined_targets(label_cctbar_a,label_cceta,tgt_info)
      call set_dependency(label_cctbar_a,label_cclg0,tgt_info)
      call set_dependency(label_cctbar_a,op_tbar_a,tgt_info)
      call set_dependency(label_cctbar_a,op_eta,tgt_info)
      call set_rule(label_cctbar_a,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cctbar_a,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cctbar_a
      labels(2) = label_cceta
      labels(3) = label_cctbar_a
      labels(4) = op_tbar_a
      labels(5) = op_eta
      labels(6) = op_tbar
      call set_rule(label_cctbar_a,ttype_frm,LEQ_SPLIT,
     &              labels,6,2,
     &              title_cctbar_a,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc1dens
      labels(2) = label_cclg0
      labels(3) = op_1dens
      labels(4) = op_ham
      labels(5) = ' '
      call add_target(label_cc1dens,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc1dens,label_cclg0,tgt_info)
      call set_dependency(label_cc1dens,op_tbar_a,tgt_info)
      call set_dependency(label_cc1dens,op_eta,tgt_info)
      call set_rule(label_cc1dens,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc1dens,1,tgt_info)

      icnt = is_keyword_set('calculate.excitation')
      if (icnt.gt.0) then
      ! right Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc_a_r
      labels(2) = label_ccrs0
      labels(3) = op_a_r
      labels(4) = op_top
      labels(5) = op_r
      call add_target(label_cc_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc_a_r,label_ccrs0,tgt_info)
      call set_dependency(label_cc_a_r,op_a_r,tgt_info)
      call set_dependency(label_cc_a_r,op_r,tgt_info)
      call set_rule(label_cc_a_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_a_r,1,tgt_info)

      ! left Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc_l_a
      labels(2) = label_cctbar_a
      labels(3) = op_l_a
      labels(4) = op_tbar
      labels(5) = op_l
      call add_target(label_cc_l_a,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc_l_a,label_cctbar_a,tgt_info)
      call set_dependency(label_cc_l_a,op_l_a,tgt_info)
      call set_dependency(label_cc_l_a,op_l,tgt_info)
      call set_rule(label_cc_l_a,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_l_a,1,tgt_info)
      end if

      icnt = is_keyword_set('calculate.ionization')
      if (icnt.gt.0) then
      ! right Jacobian transform with IP operator
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc_a_rip
      labels(2) = label_ccrs0
      labels(3) = op_a_rip
      labels(4) = op_top
      labels(5) = op_rip
      call add_target(label_cc_a_rip,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc_a_rip,label_ccrs0,tgt_info)
      call set_dependency(label_cc_a_rip,op_a_rip,tgt_info)
      call set_dependency(label_cc_a_rip,op_rip,tgt_info)
      call set_rule(label_cc_a_rip,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_a_rip,1,tgt_info)
      end if

      ! optimized formulae:

      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC ground state:
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_ccrs0opt
      labels(2) = label_ccen0
      labels(3) = label_ccrs0
      ncat = 2
      nint = 0
      call add_target(label_ccrs0opt,ttype_frm,.false.,tgt_info)
      call set_dependency(label_ccrs0opt,label_ccen0,tgt_info)
      call set_dependency(label_ccrs0opt,label_ccrs0,tgt_info)
      call set_dependency(label_ccrs0opt,mel_omgdef,tgt_info)
      call set_dependency(label_ccrs0opt,mel_topdef,tgt_info)
      call set_dependency(label_ccrs0opt,mel_ham,tgt_info)
      call set_dependency(label_ccrs0opt,mel_ccen0def,tgt_info)      
      if (isim.eq.1) then
        nint = 1
        call set_dependency(label_ccrs0opt,label_cchhat,tgt_info)
        call set_dependency(label_ccrs0opt,mel_hhatdef,tgt_info)
        labels(4) = label_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(label_ccrs0opt,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! CC ground state left-hand eq.
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cclft0opt
      labels(2) = label_cceta
      labels(3) = label_cctbar_a
      ncat = 2
      nint = 0
      call add_target(label_cclft0opt,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cclft0opt,label_cctbar_a,tgt_info)
      call set_dependency(label_cclft0opt,label_cceta,tgt_info)
      call set_dependency(label_cclft0opt,mel_tbar_adef,tgt_info)
      call set_dependency(label_cclft0opt,mel_etadef,tgt_info)
      call set_dependency(label_cclft0opt,mel_tbardef,tgt_info)
      call set_dependency(label_cclft0opt,mel_topdef,tgt_info)
      call set_dependency(label_cclft0opt,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(label_cclft0opt,label_cchhat,tgt_info)
        call set_dependency(label_cclft0opt,mel_hhatdef,tgt_info)
        labels(4) = label_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(label_cclft0opt,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! CC ground state density
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc1dens_opt
      labels(2) = label_cc1dens
      call add_target(label_cc1dens_opt,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc1dens_opt,label_cc1dens,tgt_info)
      call set_dependency(label_cc1dens_opt,mel_tbardef,tgt_info)
      call set_dependency(label_cc1dens_opt,mel_topdef,tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule(label_cc1dens_opt,ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      icnt = is_keyword_set('calculate.excitation')
      if (icnt.gt.0) then
      ! CC right-hand Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc_a_r_opt
      labels(2) = label_cc_a_r
      ncat = 1
      nint = 0
      call add_target(label_cc_a_r_opt,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc_a_r_opt,label_cc_a_r,tgt_info)
      call set_dependency(label_cc_a_r_opt,meldef_a_rex,tgt_info)
      call set_dependency(label_cc_a_r_opt,meldef_rex,tgt_info)
      call set_dependency(label_cc_a_r_opt,mel_topdef,tgt_info)
      call set_dependency(label_cc_a_r_opt,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(label_cc_a_r_opt,label_cchhat,tgt_info)
        call set_dependency(label_cc_a_r_opt,mel_hhatdef,tgt_info)
        labels(3) = label_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(label_cc_a_r_opt,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      end if

      icnt = is_keyword_set('calculate.ionization')
      if (icnt.gt.0) then
      ! CC right-hand Jacobian transform with IP operator
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = label_cc_a_rip_opt
      labels(2) = label_cc_a_rip
      ncat = 1
      nint = 0
      call add_target(label_cc_a_rip_opt,ttype_frm,.false.,tgt_info)
      call set_dependency(label_cc_a_rip_opt,label_cc_a_rip,tgt_info)
      call set_dependency(label_cc_a_rip_opt,meldef_a_rip,tgt_info)
      call set_dependency(label_cc_a_rip_opt,meldef_rip,tgt_info)
      call set_dependency(label_cc_a_rip_opt,mel_topdef,tgt_info)
      call set_dependency(label_cc_a_rip_opt,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(label_cc_a_rip_opt,label_cchhat,tgt_info)
        call set_dependency(label_cc_a_rip_opt,mel_hhatdef,tgt_info)
        labels(3) = label_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(label_cc_a_rip_opt,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      end if

      return
      end

