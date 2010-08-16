*----------------------------------------------------------------------*
      subroutine set_ccr12f_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC-R12 calculations with
*     fixed geminals
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'opdim.h'

      include 'ifc_input.h'
      include 'ifc_targets.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     min_rank, max_rank, ansatz, nop, fix_new,
     &     isim, ncat, nint, icnt, ndef, extend, r12op, RGRc,
     &     isym, ms, msc, msc_s, sym_arr(8), nlabel, ncnt,
     &     ninproj, navoid, nconnect, nreplace, trunc_type, vring_mode,
     &     connect(20), idx_sv(20), iblkmin(20),
     &     iblkmax(20),
     &     occ_def(ngastp,2,20)
      logical ::
     &     needed, r12fix, set_tp, set_tpp, screen, xsp_opt1,
     &     fixed_gem, pf12_trunc, f12x_trunc, pert, use_CS
      character(8) ::
     &     approx, f12x_mode
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)

      character ::
     &     op_bprc*4 = 'BPRC',
     &     op_xprc*4 = 'XPRC',
     &     me_bprc*8 = 'BPRCLIST',
     &     me_xprc*8 = 'XPRCLIST',
     &     medef_bprc*12 = 'DEF-BPRCLIST',
     &     medef_xprc*12 = 'DEF-XPRCLIST'
      character, parameter ::
     &     form_r_t*8 = 'FORM_R_T',
     &     op_r_t*6   = 'OP_R_T',
     &     me_r_t*6   = 'ME_R_T',
     &     medef_r_t*10  = 'DEF-ME_R_T'
      
      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for CC-R12 ...'

      msc = +1  ! assuming closed shell

      approx    = '        '
      f12x_mode = '        '
      ! read keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','fix_new',ival=fix_new)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)
      call get_argument_value('method.R12','screen',lval=screen)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      call get_argument_value('method.R12','vring',ival=vring_mode)
      call get_argument_value('method.R12','use_CS',lval=use_CS)
      call get_argument_value('method.R12','pert',lval=pert)
      call get_argument_value('method.R12','f12x',str=f12x_mode)
      call get_argument_value('method.R12','xsp1',lval=xsp_opt1)
      call get_argument_value('method.R12','RGRcouple',ival=RGRc)
      f12x_trunc = len_trim(f12x_mode).gt.0
      if (f12x_trunc) trunc_type = 0
      if (f12x_trunc) screen = .true.
      pf12_trunc = trunc_type.eq.0
c dbg
      print *,'trunc_Type= ',trunc_type
      print *,'pf12_trunc = ',pf12_trunc
c dbg

      call get_argument_value('calculate.routes','simtraf',ival=isim)

      set_tp = extend.gt.0.or.r12op.eq.1.or.r12op.ge.3
      set_tpp = r12op.eq.2.or.r12op.ge.3

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! pertubational approach: define conventional L and E
      ! Lagrange functional 
      call add_target(op_cclg,ttype_op,.false.,tgt_info)
      call set_rule(op_cclg,ttype_op,DEF_SCALAR,
     &              op_cclg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy 
      call add_target(op_ccen,ttype_op,.false.,tgt_info)
      call set_rule(op_ccen,ttype_op,DEF_SCALAR,
     &              op_ccen,1,1,
     &              parameters,0,tgt_info)


      ! Lagrange functional
      call add_target(op_ccr12lg,ttype_op,.false.,tgt_info)
      call set_rule(op_ccr12lg,ttype_op,DEF_SCALAR,
     &              op_ccr12lg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy
      call add_target(op_ccr12en,ttype_op,.false.,tgt_info)
      call set_rule(op_ccr12en,ttype_op,DEF_SCALAR,
     &              op_ccr12en,1,1,
     &              parameters,0,tgt_info)

      ! prototype formula for metric:
      call add_target(op_ccr12s0,ttype_op,.false.,tgt_info)
      call set_rule(op_ccr12s0,ttype_op,DEF_SCALAR,
     &              op_ccr12s0,1,1,
     &              parameters,0,tgt_info)
      
      if(set_tp)then
        call add_target(op_omgcex,ttype_op,.false.,tgt_info)
        call set_dependency(op_omgcex,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.false.)
        call set_rule(op_omgcex,ttype_op,CLONE_OP,
     &                op_omgcex,1,1,
     &                parameters,1,tgt_info)
      end if

      if(set_tpp)then
        call add_target(op_omgcexx,ttype_op,.false.,tgt_info)
        call set_dependency(op_omgcexx,op_cexx,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cexx,.false.)
        call set_rule(op_omgcexx,ttype_op,CLONE_OP,
     &                op_omgcexx,1,1,
     &                parameters,1,tgt_info)
      end if

      call add_target(op_bprc,ttype_op,.false.,tgt_info)
      occ_def=0
      ndef = 1
      if (extend.eq.1.or.r12op.eq.1) then
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
      else if (extend.eq.2) then
        ndef = 1
      else if (extend.eq.3.or.extend.eq.4) then
        ndef = 2
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
      else if (extend.gt.4.or.r12op.ge.2) then
        ndef = 3
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
        occ_def(IPART,1,3) = 2
        occ_def(IPART,2,3) = 2
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule(op_bprc,ttype_op,DEF_OP_FROM_OCC,
     &              op_bprc,1,1,
     &              parameters,2,tgt_info)

      call add_target(op_xprc,ttype_op,.false.,tgt_info)
      occ_def=0
      ndef = 1
      if (extend.eq.1.or.r12op.eq.1) then
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
      else if (extend.eq.2) then
        ndef = 1
      else if (extend.eq.3.or.extend.eq.4) then
        ndef = 2
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
      else if (extend.gt.4.or.r12op.ge.2) then
        ndef = 3
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
        occ_def(IPART,1,3) = 2
        occ_def(IPART,2,3) = 2
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule(op_xprc,ttype_op,DEF_OP_FROM_OCC,
     &              op_xprc,1,1,
     &              parameters,2,tgt_info)

      ! T vector times Metric
      call add_target(op_s_t,ttype_op,.false.,tgt_info)
      call set_dependency(op_s_t,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_s_t,ttype_op,CLONE_OP,
     &              op_s_t,1,1,
     &              parameters,1,tgt_info)
      
      ! T' vector times Metric
      if (r12op.gt.0) then
        call add_target(op_s_c,ttype_op,.false.,tgt_info)
        call set_dependency(op_s_c,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
        call set_rule(op_s_c,ttype_op,CLONE_OP,
     &              op_s_c,1,1,
     &              parameters,1,tgt_info)

        call add_target(op_evs,ttype_op,.false.,tgt_info)
        call set_dependency(op_evs,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
        call set_rule(op_evs,ttype_op,CLONE_OP,
     &              op_evs,1,1,
     &              parameters,1,tgt_info)

        call add_target(op_s_evs,ttype_op,.false.,tgt_info)
        call set_dependency(op_s_evs,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
        call set_rule(op_s_evs,ttype_op,CLONE_OP,
     &              op_s_evs,1,1,
     &              parameters,1,tgt_info)
      end if

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      ! perturbative: set conv. CC Lagrangian
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cclg0
      labels(2) = op_cclg
      labels(3) = op_ham
      labels(4) = op_tbar
      labels(5) = op_top
      call add_target(form_cclg0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cclg0,op_cclg,tgt_info)
      call set_dependency(form_cclg0,op_ham,tgt_info)
      call set_dependency(form_cclg0,op_tbar,tgt_info)
      call set_dependency(form_cclg0,op_top,tgt_info)
      call set_rule(form_cclg0,ttype_frm,DEF_CC_LAGRANGIAN,
     &              labels,5,1,
     &              title_cclg0,1,tgt_info)
      call set_rule(form_cclg0,ttype_frm,TEX_FORMULA,
     &              labels,5,1,
     &              'cc_lag.tex',1,tgt_info)
  
      ! perturbative: define conv. E
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccen0
      labels(2) = form_cclg0
      labels(3) = op_ccen
      labels(4) = op_tbar
      call add_target(form_ccen0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccen0,form_cclg0,tgt_info)
      call set_dependency(form_ccen0,op_ccen,tgt_info)
      call set_rule(form_ccen0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ccen0,1,tgt_info)

      ! perturbative: define conv. residual
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccrs0
      labels(2) = form_cclg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_ccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccrs0,form_cclg0,tgt_info)
      call set_dependency(form_ccrs0,op_omg,tgt_info)
      call set_rule(form_ccrs0,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccrs0,1,tgt_info)

      ! now R12 starts
      call add_target(form_ccr12lg0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Lagrangian (in 'complete' basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0
      labels(2) = op_ccr12lg
      labels(3) = op_ham
      labels(4) = op_r12
      labels(5) = op_r12
      labels(6) = op_tbar
      labels(7) = op_top
      nlabel = 7
      if (set_tp) then
        call set_dependency(form_ccr12lg0,op_cex,tgt_info)
        call set_dependency(form_ccr12lg0,op_cexbar,tgt_info)
        if (.not.xsp_opt1) then
          labels(nlabel+1) = op_cexbar
          labels(nlabel+2) = op_cex
        else
          labels(nlabel+1) = 'T12FBAR'
          labels(nlabel+2) = 'T12FML'
        end if
        nlabel = nlabel+2
      end if
      if (set_tpp) then
        call set_dependency(form_ccr12lg0,op_cexx,tgt_info)
        call set_dependency(form_ccr12lg0,op_cexxbar,tgt_info)
        labels(nlabel+1) = op_cexxbar
        labels(nlabel+2) = op_cexx
        nlabel = nlabel+2
      end if
      call set_dependency(form_ccr12lg0,op_ccr12lg,tgt_info)
      call set_dependency(form_ccr12lg0,op_ham,tgt_info)
      call set_dependency(form_ccr12lg0,op_r12,tgt_info)
c      call set_dependency(form_ccr12lg0,op_rba,tgt_info)
      call set_dependency(form_ccr12lg0,op_tbar,tgt_info)
      call set_dependency(form_ccr12lg0,op_top,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,ansatz,'---')
      call set_rule(form_ccr12lg0,ttype_frm,DEF_CCR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)
      ! (a-I) apply additional truncation schemes
      if (f12x_trunc) then
        call set_rule2(form_ccr12lg0,SELECT_SPECIAL,tgt_info)
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'LABEL_RES',
     &       1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'LABEL_IN',
     &       1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'OPERATORS',
     &       4,tgt_info,
     &       val_label=(/op_tbar,op_ham,op_top,op_r12/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'TYPE',
     &       1,tgt_info,
     &       val_str='f12x')
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'MODE',
     &       1,tgt_info,
     &       val_str=trim(f12x_mode))
      end if
      ! special XSPopt with tildeT_1 only
      if (xsp_opt1) then
        ! process formula
        call set_rule2(form_ccr12lg0,SELECT_SPECIAL,tgt_info)
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'LABEL_RES',
     &       1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'LABEL_IN',
     &       1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'OPERATORS',
     &       6,tgt_info,
     &       val_label=(/op_tbar,'T12FBAR',op_cexbar,
     &                   op_top, 'T12FML', op_cex/))
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'TYPE',
     &       1,tgt_info,
     &       val_str='OPT1')
        call set_arg(form_ccr12lg0,SELECT_SPECIAL,'MODE',
     &       1,tgt_info,
     &       val_str='unused')
        ! replace remaining formal T12BAR and T12 by their actual def's
      end if
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0 ! output formula (itself)
      labels(2) = form_ccr12lg0 ! input formula
      labels(3) = form_r12_vint    ! the intermediates to be factored
      labels(4) = form_r12_vint//'^+'
      labels(5) = form_r12_bint
      labels(6) = form_r12_bhint
      labels(7) = form_r12_xint
      nint = 5
      call set_dependency(form_ccr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bhint,tgt_info)
      if (ansatz.ne.1) then
        labels(8) = form_r12_cint
        labels(9) = trim(form_r12_cint)//'^+'
c        labels(10) = form_r12_xpint
        call set_dependency(form_ccr12lg0,form_r12_cint,tgt_info)
c        call set_dependency(form_ccr12lg0,form_r12_xpint,tgt_info)
        nint = 7
      end if
      if (vring_mode.gt.0) then
        call set_dependency(form_ccr12lg0,'Vring_formal',tgt_info)
        labels(2+nint+1) = 'Vring_formal'
        labels(2+nint+2) = 'Vring_formal^+'
        if (vring_mode.eq.2) then
          labels(2+nint+3) = 'Vring2_formal^+'
          nint = nint+3
        else
          nint = nint+2
        end if
      end if
      if (use_CS) then
        call set_dependency(form_ccr12lg0,'C1_formal',tgt_info)
        labels(2+nint+1) = 'C1_formal'
        nint = nint+1
      end if 
      if (.not.pf12_trunc) then
        call set_dependency(form_ccr12lg0,form_r12_xhint,tgt_info)
        call set_dependency(form_ccr12lg0,form_r12_pint,tgt_info)
        call set_dependency(form_ccr12lg0,form_r12_zint,tgt_info)
        labels(2+nint+1) = form_r12_xhint
        labels(2+nint+2) = form_r12_pint 
        labels(2+nint+3) = form_r12_zint 
        nint = nint+3
        if (.not.screen) then
          call set_dependency(form_ccr12lg0,'Vpx_formal',tgt_info)
          labels(2+nint+1) = 'Vpx_formal'
          labels(2+nint+2) = 'Vpx_formal^+'
          nint = nint+2
        end if
      end if
      if (r12op.gt.0.and.max_rank.gt.2) then
        if (.not.screen) then
          call set_dependency(form_ccr12lg0,'Vpx_formal',tgt_info)
          call set_dependency(form_ccr12lg0,'Z2INT_R12',tgt_info)
          call set_dependency(form_ccr12lg0,form_r12_xhint,tgt_info)
          labels(2+nint+1) = 'Vpx_formal'
          labels(2+nint+2) = 'Vpx_formal^+'
          labels(2+nint+3)  = 'Z2INT_R12'
          labels(2+nint+4)  = 'Z2INT_R12^+'
          labels(2+nint+5) = form_r12_xhint
          nint = nint+5
        else
          call set_dependency(form_ccr12lg0,'Z2INT_R12',tgt_info)
          call set_dependency(form_ccr12lg0,form_r12_xhint,tgt_info)
          labels(2+nint+1)  = 'Z2INT_R12'
          labels(2+nint+2)  = 'Z2INT_R12^+'
          labels(2+nint+3) = form_r12_xhint
          nint = nint+3
        end if
      end if
 
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,nint,'---')
      call set_rule(form_ccr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)
      ! (c) post-processing: remove terms which do not contribute for
      !     the given R12-approximation
      ! .... to come
      ! fix for A
      if (trim(approx).eq.'A') then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12lg0
        labels(2) = form_ccr12lg0
        labels(3) = op_ccr12lg
        labels(4) = op_x_inter
        call set_rule(form_ccr12lg0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_ccr12lg0,1,tgt_info)
      end if

      ! screen all remaining R12 terms
      ! (argument: would vanish in standard approximation)
      if (screen) then        
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12lg0
        labels(2) = form_ccr12lg0
        labels(3) = op_ccr12lg
        labels(4) = op_r12
        labels(5) = op_r12//'^+'
        call set_rule(form_ccr12lg0,ttype_frm,INVARIANT,
     &              labels,5,1,
     &              title_ccr12lg0,1,tgt_info)
      end if

      ! there remain a few unprocessed R12 contributions
      ! for ansatz > 1
      ! as a first resort we replace r12 by the actual integrals
      if (ansatz.gt.1) then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12lg0
        labels(2) = form_ccr12lg0
        labels(3) = op_r12
        labels(4) = op_rint
        labels(5) = op_r12//'^+'
        labels(6) = op_rint//'^+'
        call set_dependency(form_ccr12lg0,op_rint,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_ccr12lg0,2,'---')
        call set_rule(form_ccr12lg0,ttype_frm,REPLACE,
     &              labels,6,1,
     &              parameters,2,tgt_info)
      end if
      ! perturbative:
      if (pert) then
        ! replace L by T^+
        call set_rule2(form_ccr12lg0,REPLACE,tgt_info)
        call set_arg(form_ccr12lg0,REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/op_tbar,op_top//'^+'/))
        call set_arg(form_ccr12lg0,REPLACE,'TITLE',1,tgt_info,
     &       val_str='CCR12 Lagrangian for pert. eval.')
        ! optional: factor out conv. energy, residual
        ! take only part that does not involve conv. energy, residual
        call set_dependency(form_ccr12lg0,form_ccen0,tgt_info)
        call set_dependency(form_ccr12lg0,form_ccrs0,tgt_info)
        call set_rule2(form_ccr12lg0,FACTOR_OUT,tgt_info)
        call set_arg(form_ccr12lg0,FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        Call set_arg(form_ccr12lg0,FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,FACTOR_OUT,'INTERM',2,tgt_info,
     &       val_label=(/form_ccen0,form_ccrs0/))
        call set_arg(form_ccr12lg0,FACTOR_OUT,'TITLE',1,tgt_info,
     &       val_str='CCR12 Lagrangian for pert. eval.')
        Call set_rule2(form_ccr12lg0,INVARIANT,tgt_info)
        call set_arg(form_ccr12lg0,INVARIANT,'LABEL_RES',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))        
        call set_arg(form_ccr12lg0,INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/form_ccr12lg0/))
        call set_arg(form_ccr12lg0,INVARIANT,'OP_RES',1,tgt_info,
     &       val_label=(/op_ccr12lg/))        
        call set_arg(form_ccr12lg0,INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/op_ccen,op_omg/))  
        call set_arg(form_ccr12lg0,INVARIANT,'TITLE',1,tgt_info,
     &       val_str='CCR12 Lagrangian for pert. eval.')
      end if      
c dbg
      call form_parameters(-2,parameters,2,
     &       'stdout',0,'---')
      call set_rule(form_ccr12lg0,ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)
c dbg
      call set_rule(form_ccr12lg0,ttype_frm,TEX_FORMULA,
     &              labels,5,1,
     &              'ccr12f_lag.tex',1,tgt_info)
      
      call add_target(form_ccr12_s0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Metric (in 'complete' basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12_s0
      labels(2) = op_ccr12s0
      labels(3) = op_r12
      labels(4) = op_tbar
      labels(5) = op_top
      nlabel = 5
      if (set_tp) then
        call set_dependency(form_ccr12_s0,op_cex,tgt_info)
        call set_dependency(form_ccr12_s0,op_cexbar,tgt_info)
        if (.not.xsp_opt1) then
          labels(nlabel+1) = op_cexbar
          labels(nlabel+2) = op_cex
        else
          labels(nlabel+1) = 'T12FBAR'
          labels(nlabel+2) = 'T12FML'
        end if
        nlabel = nlabel+2
      end if
      if (set_tpp) then
        call set_dependency(form_ccr12_s0,op_cexx,tgt_info)
        call set_dependency(form_ccr12_s0,op_cexxbar,tgt_info)
        labels(nlabel+1) = op_cexxbar
        labels(nlabel+2) = op_cexx
        nlabel = nlabel+2
      end if
      call set_dependency(form_ccr12_s0,op_ccr12s0,tgt_info)
      call set_dependency(form_ccr12_s0,op_r12,tgt_info)
      call set_dependency(form_ccr12_s0,op_tbar,tgt_info)
      call set_dependency(form_ccr12_s0,op_top,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,ansatz,'---')
      call set_rule(form_ccr12_s0,ttype_frm,DEF_CCR12_METRIC,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)
      if (xsp_opt1) then
        ! process formula
        call set_rule2(form_ccr12_s0,SELECT_SPECIAL,tgt_info)
        call set_arg(form_ccr12_s0,SELECT_SPECIAL,'LABEL_RES',
     &       1,tgt_info,
     &       val_label=(/form_ccr12_s0/))
        call set_arg(form_ccr12_s0,SELECT_SPECIAL,'LABEL_IN',
     &       1,tgt_info,
     &       val_label=(/form_ccr12_s0/))
        call set_arg(form_ccr12_s0,SELECT_SPECIAL,'OPERATORS',
     &       6,tgt_info,
     &       val_label=(/op_tbar,'T12FBAR',op_cexbar,
     &                   op_top, 'T12FML', op_cex/))
        call set_arg(form_ccr12_s0,SELECT_SPECIAL,'TYPE',
     &       1,tgt_info,
     &       val_str='OPT1')
        call set_arg(form_ccr12_s0,SELECT_SPECIAL,'MODE',
     &       1,tgt_info,
     &       val_str='unused')
        ! replace remaining formal T12BAR and T12 by their actual def's
      end if
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12_s0 ! output formula (itself)
      labels(2) = form_ccr12_s0 ! input formula
      labels(3) = form_r12_xint
      nint = 1
      call set_dependency(form_ccr12_s0,form_r12_xint,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,nint,'---')
      call set_rule(form_ccr12_s0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)
      ! there remain a few unprocessed R12 contributions
      ! for ansatz > 1
      ! as a first resort we replace r12 by the actual integrals
      if (ansatz.gt.1) then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s0
        labels(2) = form_ccr12_s0
        labels(3) = op_r12
        labels(4) = op_rint
        labels(5) = op_r12//'^+'
        labels(6) = op_rint//'^+'
        call set_dependency(form_ccr12_s0,op_rint,tgt_info)
        call form_parameters(-1,
     &       parameters,2,title_ccr12lg0,2,'---')
        call set_rule(form_ccr12_s0,ttype_frm,REPLACE,
     &              labels,6,1,
     &              parameters,2,tgt_info)
      end if
c dbg
      call form_parameters(-2,parameters,2,
     &       'stdout',0,'---')
      call set_rule(form_ccr12_s0,ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)
c dbg


      ! formula for obtaining eigenvalues of R12 metric: S.VEC
      if (.not.r12fix) then
        call add_target(form_ccr12_s_v,ttype_frm,.false.,tgt_info)
        call set_dependency(form_ccr12_s_v,form_ccr12_s0,tgt_info)
        call set_dependency(form_ccr12_s_v,op_s_c,tgt_info)
        call set_dependency(form_ccr12_s_v,op_evs,tgt_info)
        call set_dependency(form_ccr12_s_v,op_s_evs,tgt_info)
        call set_dependency(form_ccr12_s_v,op_cex,tgt_info)
        call set_dependency(form_ccr12_s_v,op_cexbar,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = 'FORM_SCRATCH'
        labels(2) = form_ccr12_s0
        labels(3) = op_s_c
        labels(4) = op_cexbar
        labels(5) = ' '
        call set_rule(form_ccr12_s_v,ttype_frm,DERIVATIVE,
     &               labels,5,1,
     &               'scratch',1,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s_v
        labels(2) = 'FORM_SCRATCH'
        labels(3) = op_s_evs
        labels(4) = op_cex
        labels(5) = op_evs
        call set_rule(form_ccr12_s_v,ttype_frm,DERIVATIVE,
     &               labels,5,1,
     &               'S.EV',1,tgt_info)

      end if
      
      call add_target(form_ccr12en0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccr12en0,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12en0,op_ccr12en,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12en0
      labels(2) = form_ccr12lg0
      labels(3) = op_ccr12en
      labels(4) = op_tbar
      nlabel = 4
      ! test - fix_new: replace L' -> L or T^+ only now
      if (set_tp.and.fix_new.eq.0) then
        labels(nlabel+1) = op_cexbar
        nlabel = nlabel+1
      end if
      if (set_tpp.and.fix_new.eq.0) then
        labels(nlabel+1) = op_cexxbar
        nlabel = nlabel+1
      end if
      call set_rule(form_ccr12en0,ttype_frm,INVARIANT,
     &              labels,nlabel,1,
     &              title_ccr12en0,1,tgt_info)
      if ((set_tp.or.set_tpp).and.fix_new.gt.0) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12en0
        labels(2) = form_ccr12en0
        nlabel = 2
        nreplace = 0
        if (set_tp) then
          labels(nlabel+1) = op_cexbar
          if (fix_new.eq.1) then
            labels(nlabel+2) = op_tbar
          else
            labels(nlabel+2) = op_top//'^+'
          end if
          labels(nlabel+3) = op_cex
          labels(nlabel+4) = op_top
          nlabel = nlabel+4
          nreplace = nreplace+2
        end if
        if (set_tpp) then
          labels(nlabel+1) = op_cexxbar
          if (fix_new.eq.1) then
            labels(nlabel+2) = op_tbar
          else
            labels(nlabel+2) = op_top//'^+'
          end if
          labels(nlabel+3) = op_cexx
          labels(nlabel+4) = op_top
          nlabel = nlabel+4
          nreplace = nreplace+1
        end if
        call form_parameters(-1,
     &       parameters,2,title_ccr12en0,nreplace,'---')
        call set_rule(form_ccr12en0,ttype_frm,REPLACE,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)        
        if (fix_new.gt.1) then
          labels(1) = form_ccr12en0
          labels(2) = form_ccr12en0
          labels(3) = op_ccr12en
          call set_rule(form_ccr12en0,ttype_frm,SUM_HERMIT,
     &              labels,3,1,
     &              title_ccr12en0,1,tgt_info)
        end if
      end if


      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12rs_t
      labels(2) = form_ccr12lg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_ccr12rs_t,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccr12rs_t,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12rs_t,op_omg,tgt_info)
      call set_rule(form_ccr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccr12rs_t,1,tgt_info)
      if ((set_tp.or.set_tpp).and.fix_new.gt.0) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12rs_t
        labels(2) = form_ccr12rs_t
        nlabel = 2
        nreplace = 0
        if (set_tp) then
          labels(nlabel+1) = op_cex
          labels(nlabel+2) = op_top
          nlabel = nlabel+2
          nreplace = nreplace+1
        end if
        if (set_tpp) then
          labels(nlabel+1) = op_cexx
          labels(nlabel+2) = op_top
          nlabel = nlabel+2
          nreplace = nreplace+1
        end if
        call form_parameters(-1,
     &       parameters,2,title_ccr12rs_t,nreplace,'---')
        call set_rule(form_ccr12rs_t,ttype_frm,REPLACE,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)        
      end if

      call add_target(op_r_t,ttype_frm,.false.,tgt_info)
      ndef = 1
      call r12gem_parameters(-1,parameters,
     &     0,2,max_rank,ansatz)
      call set_rule(op_r_t,ttype_op,DEF_R12GEMINAL,
     &              op_r_t,1,1,
     &              parameters,1,tgt_info)

      if (extend.gt.0) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_r_t
        labels(2) = op_r_t
        labels(3) = op_r_t
        labels(4) = op_rint
c     labels(5) = op_top
        labels(5) = op_cex
        labels(6) = op_r_t
        idx_sv(1:4) = (/1,2,3,1/)
        iblkmin(1:4) = (/1,10,1,1/)
        iblkmax(1:4) = (/0,10,0,0/)
        nconnect = 1
        connect(1:2) = (/2,3/)
        navoid = 0
        ninproj = 0
        call add_target(form_r_t,ttype_frm,.false.,tgt_info)
        call set_dependency(form_r_t,op_rint,tgt_info)
c     call set_dependency(form_r_t,op_top,tgt_info)
        call set_dependency(form_r_t,op_cex,tgt_info)
        call set_dependency(form_r_t,op_r_t,tgt_info)
        call expand_parameters(-1,
     &       parameters,3,
     &       'XXX',4,idx_sv,iblkmin,iblkmax,
     &       connect,nconnect,
     &       0,navoid,
     &       0,ninproj)
        call set_rule(form_r_t,ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,6,1,
     &              parameters,3,tgt_info)
      end if

      if(set_tp)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12rs_c
        labels(2) = form_ccr12lg0
        labels(3) = op_omgcex
        labels(4) = op_cexbar
        labels(5) = ' '
        call add_target(form_ccr12rs_c,ttype_frm,.false.,tgt_info)
        call set_dependency(form_ccr12rs_c,form_ccr12lg0,tgt_info)
        call set_dependency(form_ccr12rs_c,op_omgcex,tgt_info)
        call set_rule(form_ccr12rs_c,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_ccr12rs_c,1,tgt_info)
      endif

      if(set_tpp)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_ccr12rs_cpp
        labels(2) = form_ccr12lg0
        labels(3) = op_omgcexx
        labels(4) = op_cexxbar
        labels(5) = ' '
        call add_target(form_ccr12rs_cpp,ttype_frm,.false.,tgt_info)
        call set_dependency(form_ccr12rs_cpp,form_ccr12lg0,tgt_info)
        call set_dependency(form_ccr12rs_cpp,op_omgcexx,tgt_info)
        call set_rule(form_ccr12rs_cpp,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_ccr12rs_cpp,1,tgt_info)
      endif


*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! perturbative: set conv. E + residual opt. formula
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_ccrs0
      labels(2) = form_ccen0
      labels(3) = form_ccrs0
      ncat = 2
      nint = 0
      call add_target(fopt_ccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ccrs0,form_ccen0,tgt_info)
      call set_dependency(fopt_ccrs0,form_ccrs0,tgt_info)
      call set_dependency(fopt_ccrs0,mel_omgdef,tgt_info)
      call set_dependency(fopt_ccrs0,mel_topdef,tgt_info)
      call set_dependency(fopt_ccrs0,mel_ham,tgt_info)
      call set_dependency(fopt_ccrs0,mel_ccen0def,tgt_info)      
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_ccrs0,form_cchhat,tgt_info)
        call set_dependency(fopt_ccrs0,mel_hhatdef,tgt_info)
        labels(4) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_ccrs0,form_cchbar,tgt_info)
        call set_dependency(fopt_ccrs0,meldef_hbar,tgt_info)
        labels(4) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ccrs0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! set opt. R12 pert. expr.
      call add_target2('L_CCR12_PT-OPT',.false.,tgt_info)
      call set_dependency('L_CCR12_PT-OPT',form_ccr12lg0,tgt_info)
      call set_dependency('L_CCR12_PT-OPT',mel_ccr12lg0def,tgt_info)
      call set_dependency('L_CCR12_PT-OPT',mel_topdef,tgt_info)
      call set_dependency('L_CCR12_PT-OPT',mel_ham,tgt_info)
      if (vring_mode.gt.0)
     &       call set_dependency('L_CCR12_PT-OPT','Vring-EVAL',tgt_info)
      if (use_CS)
     &       call set_dependency('L_CCR12_PT-OPT','C1-EVAL',tgt_info)
      if (.not.pf12_trunc) then
        call set_dependency('L_CCR12_PT-OPT',mel_p_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_z_def,tgt_info)      
      end if
      if (max_rank.ge.2) then
        call set_dependency('L_CCR12_PT-OPT',mel_rint,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_v_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_b_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_bh_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_x_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_xh_def,tgt_info)      
        call set_dependency('L_CCR12_PT-OPT',mel_c_def,tgt_info)      
      end if
      if (max_rank.ge.2) then
        call set_dependency('L_CCR12_PT-OPT','DEF-Z2LIST',tgt_info)
        if (.not.screen)
     &       call set_dependency('L_CCR12_PT-OPT','Vpx-INTER',tgt_info)
      end if
      call set_rule2('L_CCR12_PT-OPT',OPTIMIZE,tgt_info)
      call set_arg('L_CCR12_PT-OPT',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &     val_label=(/'L_CCR12_PT-OPT'/))
      call set_arg('L_CCR12_PT-OPT',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &     val_label=(/form_ccr12lg0/))      

      ! CC ground state:
      fixed_gem =  r12fix.or.fix_new.gt.0
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_ccr12_0
      labels(2) = form_ccr12en0
      labels(3) = form_ccr12rs_t
      ncat = 2
      if(set_tp.and..not.fixed_gem)then
        labels(1+ncat+1) = form_ccr12rs_c
        ncat = ncat+1
      endif
      if(set_tpp.and..not.fixed_gem)then
        labels(1+ncat+1) = form_ccr12rs_cpp
        ncat = ncat+1
      endif
      nint = 0
      call add_target(fopt_ccr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12en0,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12rs_t,tgt_info)
      if(set_tp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,form_ccr12rs_c,tgt_info)
      if(set_tpp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,form_ccr12rs_cpp,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_topdef,tgt_info)
      if(set_tp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,mel_cex_def,tgt_info)
      if(set_tpp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,mel_cexx_def,tgt_info)
      if(set_tp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,mel_omgcexdef,tgt_info)
      if(set_tpp.and..not.fixed_gem)
     &     call set_dependency(fopt_ccr12_0,mel_omgcexxdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_ham,tgt_info)
      if (vring_mode.gt.0)
     &       call set_dependency(fopt_ccr12_0,'Vring-EVAL',tgt_info)
      if (use_CS)
     &       call set_dependency(fopt_ccr12_0,'C1-EVAL',tgt_info)
      if (.not.pf12_trunc) then
        call set_dependency(fopt_ccr12_0,mel_p_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_z_def,tgt_info)      
      end if
      if (max_rank.ge.2) then
        call set_dependency(fopt_ccr12_0,mel_rint,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_v_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_b_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_bh_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_x_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_xh_def,tgt_info)      
        call set_dependency(fopt_ccr12_0,mel_c_def,tgt_info)      
      end if
      if (max_rank.ge.2) then
        call set_dependency(fopt_ccr12_0,'DEF-Z2LIST',tgt_info)
        if (.not.screen)
     &       call set_dependency(fopt_ccr12_0,'Vpx-INTER',tgt_info)
      end if
      call set_dependency(fopt_ccr12_0,mel_ccr12en0def,tgt_info)      
      if (isim.eq.1) then
        nint = nint+ 1
        call set_dependency(fopt_ccr12_0,form_cchhat,tgt_info)
        call set_dependency(fopt_ccr12_0,mel_hhatdef,tgt_info)
        labels(ncat+nint+1) = form_cchhat
      else if (isim.eq.2) then
        nint = nint+ 1
        call set_dependency(fopt_ccr12_0,form_cchbar,tgt_info)
        call set_dependency(fopt_ccr12_0,meldef_hbar,tgt_info)
        labels(ncat+nint+1) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ccr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ncnt = is_keyword_set('calculate.check_S')
      if (ncnt.gt.0.and..not.r12fix) then
        call add_target(fopt_ccr12_s_v,ttype_frm,.false.,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = fopt_ccr12_s_v
        labels(2) = form_ccr12_s_v
        ncat = 1
        nint = 0
        call set_dependency(fopt_ccr12_s_v,form_ccr12_s_v,tgt_info)
        call set_dependency(fopt_ccr12_s_v,meldef_s_evs,tgt_info)
        call set_dependency(fopt_ccr12_s_v,meldef_evs,tgt_info)
        call set_dependency(fopt_ccr12_s_v,mel_x_def,tgt_info)
        call opt_parameters(-1,parameters,ncat,nint)
        call set_rule(fopt_ccr12_s_v,ttype_frm,OPTIMIZE,
     &       labels,ncat+nint+1,1,
     &       parameters,1,tgt_info)
      end if
        
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! perturbative: E0 for conv. CC
      call add_target(mel_ccen0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccen0def,op_ccen,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ccen0
      labels(2) = op_ccen
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ccen0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! L0/E0:
      call add_target(mel_ccr12lg0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccr12lg0def,op_ccr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_ccr12lg0
      labels(2) = op_ccr12lg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ccr12lg0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target(mel_ccr12en0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccr12en0def,op_ccr12en,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_ccr12en0
      labels(2) = op_ccr12en
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ccr12en0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target(medef_r_t,ttype_opme,.false.,tgt_info)
      call set_dependency(medef_r_t,op_r_t,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_r_t
      labels(2) = op_r_t
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(medef_r_t,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      if(set_tp)then
        ! OMG-R12-EXT list definition
        call add_target(mel_omgcexdef,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_omgcexdef,op_omgcex,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_omgcex
        labels(2) = op_omgcex
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_omgcexdef,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

        call add_target(mel_cex_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_cex_def,op_cex,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_cex
        labels(2) = op_cex
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_cex_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

        call add_target(mel_cexbar_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_cexbar_def,op_cexbar,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_cexbar
        labels(2) = op_cexbar
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_cexbar_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

      endif

      if(set_tpp)then
        ! OMG-R12-EXT list definition
        call add_target(mel_omgcexxdef,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_omgcexxdef,op_omgcexx,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_omgcexx
        labels(2) = op_omgcexx
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_omgcexxdef,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

        call add_target(mel_cexx_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_cexx_def,op_cexx,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_cexx
        labels(2) = op_cexx
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_cexx_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

        call add_target(mel_cexxbar_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_cexxbar_def,op_cexxbar,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_cexxbar
        labels(2) = op_cexxbar
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_cexxbar_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

      endif

      ncnt = is_keyword_set('calculate.check_S')
      if (ncnt.gt.0.and..not.r12fix) then
        call add_target(meldef_evs,ttype_opme,.false.,tgt_info)
        call add_target(meldef_s_evs,ttype_opme,.false.,tgt_info)
        call set_dependency(meldef_evs,op_evs,tgt_info)
        call set_dependency(meldef_s_evs,op_s_evs,tgt_info)
        do icnt = 1, ncnt
          call get_argument_value('calculate.check_S','sym',
     &         keycount=icnt,
     &         iarr=sym_arr)
          call get_argument_value('calculate.check_S','msc',
     &         keycount=icnt,
     &         ival=msc_s)
          do isym = 1, orb_info%nsym
            if (sym_arr(isym).eq.0) cycle

            call me_list_label(me_label,mel_evs,isym,0,0,msc_s,.false.)
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = op_evs
            call me_list_parameters(-1,parameters,
     &           msc_s,0,isym,0,0,.false.)
            call set_rule(meldef_evs,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)

            call me_list_label(me_label,mel_s_evs,
     &                         isym,0,0,msc_s,.false.)
            labels(1:20)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = op_s_evs
            call me_list_parameters(-1,parameters,
     &           msc_s,0,isym,0,0,.false.)
            call set_rule(meldef_s_evs,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)
          end do
        end do
        
      end if

      call add_target(me_bprc,ttype_opme,.false.,tgt_info)
      call set_dependency(me_bprc,op_bprc,tgt_info)
      call set_dependency(me_bprc,eval_r12_inter,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_bprc
      labels(2) = op_bprc
      call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
      call set_rule(me_bprc,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_bprc
      labels(2) = mel_b_inter
      labels(3) = mel_bh_inter
      call add_parameters(-1,parameters,
     &     2,(/1d0,1d0/),2)
      call set_rule(me_bprc,ttype_opme,ADD,
     &              labels,3,1,
     &              parameters,1,tgt_info)
      
      call add_target(me_xprc,ttype_opme,.false.,tgt_info)
      call set_dependency(me_xprc,op_xprc,tgt_info)
      call set_dependency(me_bprc,eval_r12_inter,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_xprc
      labels(2) = op_xprc
      call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
      call set_rule(me_xprc,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_xprc
      labels(2) = mel_x_inter
      call add_parameters(-1,parameters,
     &     1,(/1d0/),2)
      call set_rule(me_xprc,ttype_opme,ADD,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      needed = is_keyword_set('calculate.skip_E').eq.0

      ! totally symmetric dia for use below:
      fixed_gem =  r12fix.or.fix_new.gt.0
      if (set_tp.and..not.set_tpp.and..not.fixed_gem) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_ccr12_gs,ttype_gen,needed,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
c        call set_dependency(solve_ccr12_gs,mel_diar12,tgt_info)
        call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
        call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
        if (vring_mode.gt.0)
     &       call set_dependency(solve_ccr12_gs,'Vring-EVAL',tgt_info)
        if (use_CS)
     &       call set_dependency(solve_ccr12_gs,'C1-EVAL',tgt_info)
        if (.not.pf12_trunc)
     &       call set_dependency(solve_ccr12_gs,'EVAL_PZ',tgt_info)
        if (max_rank.ge.3.and.RGRc.ne.0)
     &    call set_dependency(solve_ccr12_gs,'Z2INT_R12_EVAL',tgt_info)
        call set_dependency(solve_ccr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_ccr12_gs,me_xprc,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cex
        labels(3) = mel_omg
        labels(4) = mel_omgcex
        labels(5) = mel_dia1    
        labels(6) = mel_dia1 !r12  
        labels(7) = mel_ccr12en0
        labels(8) = fopt_ccr12_0
        labels(9) = me_bprc
        labels(10)= me_xprc
        labels(11) = mel_ham
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
        call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)
        
      else if (.not.set_tp.and.set_tpp.and..not.fixed_gem) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_ccr12_gs,ttype_gen,needed,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
        call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
        if (vring_mode.gt.0)
     &       call set_dependency(solve_ccr12_gs,'Vring-EVAL',tgt_info)
        if (use_CS)
     &       call set_dependency(solve_ccr12_gs,'C1-EVAL',tgt_info)
        if (.not.pf12_trunc)
     &    call set_dependency(solve_ccr12_gs,'EVAL_PZ',tgt_info)
        if (max_rank.ge.3.and.RGRc.ne.0)
     &    call set_dependency(solve_ccr12_gs,'Z2INT_R12_EVAL',tgt_info)
        call set_dependency(solve_ccr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_ccr12_gs,me_xprc,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cexx
        labels(3) = mel_omg
        labels(4) = mel_omgcexx
        labels(5) = mel_dia1    ! dummy
        labels(6) = mel_dia1    ! dummy
        labels(7) = mel_ccr12en0
        labels(8) = fopt_ccr12_0
        labels(9) = me_bprc
        labels(10)= me_xprc
        labels(11) = mel_ham
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
        call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)

      else if (set_tp.and.set_tpp.and..not.fixed_gem) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_ccr12_gs,ttype_gen,needed,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
        call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
        if (vring_mode.gt.0)
     &       call set_dependency(solve_ccr12_gs,'Vring-EVAL',tgt_info)
        if (use_CS)
     &       call set_dependency(solve_ccr12_gs,'C1-EVAL',tgt_info)
        if (.not.pf12_trunc)
     &    call set_dependency(solve_ccr12_gs,'EVAL_PZ',tgt_info)
        if (max_rank.ge.3.and.RGRc.ne.0)
     &    call set_dependency(solve_ccr12_gs,'Z2INT_R12_EVAL',tgt_info)
        call set_dependency(solve_ccr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_ccr12_gs,me_xprc,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cex
        labels(3) = mel_cexx
        labels(4) = mel_omg
        labels(5) = mel_omgcex
        labels(6) = mel_omgcexx
        labels(7) = mel_dia1    
        labels(8) = mel_dia1    ! dummy
        labels(9) = mel_dia1    ! dummy
        labels(10) = mel_ccr12en0
        labels(11) = fopt_ccr12_0
        labels(12) = me_bprc
        labels(13)= me_xprc
        labels(14) = mel_ham
        call solve_parameters(-1,parameters,2, 3,1,'DIA/BLK/BLK')
        call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,14,6,
     &       parameters,2,tgt_info)
        
      else
        ! for present version of pertubation, we are correct here ...

        ! switch pert/non-pert
        if (.not.pert) then
          call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

          call add_target(solve_ccr12_gs,ttype_gen,needed,tgt_info)
          call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
          call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
          if (max_rank.ge.2) 
     &       call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
          if (vring_mode.gt.0)
     &       call set_dependency(solve_ccr12_gs,'Vring-EVAL',tgt_info)
          if (use_CS)
     &       call set_dependency(solve_ccr12_gs,'C1-EVAL',tgt_info)
          if (.not.pf12_trunc)
     &       call set_dependency(solve_ccr12_gs,'EVAL_PZ',tgt_info)
          if (max_rank.ge.3.and.RGRc.ne.0)
     &     call set_dependency(solve_ccr12_gs,'Z2INT_R12_EVAL',tgt_info)
          call solve_parameters(-1,parameters,2, 1,1,'DIA')
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = mel_top
          labels(2) = mel_omg
          labels(3) = mel_dia1
          labels(4) = mel_ccr12en0
          labels(5) = fopt_ccr12_0
          labels(6) = mel_ham
          call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &         labels,6,2,
     &         parameters,2,tgt_info)
   
        else
        ! perturbative case
          call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

          call add_target(solve_ccr12_gs,ttype_gen,needed,tgt_info)
          call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
          call set_dependency(solve_ccr12_gs,fopt_ccrs0,tgt_info)
          call set_dependency(solve_ccr12_gs,'L_CCR12_PT-OPT',tgt_info)
          if (max_rank.ge.2) 
     &       call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
          if (vring_mode.gt.0)
     &       call set_dependency(solve_ccr12_gs,'Vring-EVAL',tgt_info)
          if (use_CS)
     &       call set_dependency(solve_ccr12_gs,'C1-EVAL',tgt_info)
          if (.not.pf12_trunc)
     &       call set_dependency(solve_ccr12_gs,'EVAL_PZ',tgt_info)
          if (max_rank.ge.3.and.RGRc.ne.0)
     &     call set_dependency(solve_ccr12_gs,'Z2INT_R12_EVAL',tgt_info)
          call solve_parameters(-1,parameters,2, 1,1,'DIA')
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
          labels(1:20)(1:len_target_name) = ' '
          labels(1) = mel_top
          labels(2) = mel_omg
          labels(3) = mel_dia1
          labels(4) = mel_ccen0
          labels(5) = fopt_ccrs0
          labels(6) = mel_ham
          ! solve conv. CC equations
          call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &         labels,6,2,
     &         parameters,2,tgt_info)
          ! now pert. R12 evaluation
          call set_rule2(solve_ccr12_gs,EVAL,tgt_info)
          call set_arg(solve_ccr12_gs,EVAL,'FORM',1,tgt_info,
     &         val_label=(/'L_CCR12_PT-OPT'/))
          call set_rule2(solve_ccr12_gs,PRINT_MEL,tgt_info)
          call set_arg(solve_ccr12_gs,PRINT_MEL,'LIST',1,tgt_info,
     &         val_label=(/mel_ccen0/))
          call set_arg(solve_ccr12_gs,PRINT_MEL,'COMMENT',1,tgt_info,
     &         val_str='>>> CC energy (conv) :')
          call set_arg(solve_ccr12_gs,PRINT_MEL,'FORMAT',1,tgt_info,
     &         val_str='SCAL F20.12')
          call set_rule2(solve_ccr12_gs,PRINT_MEL,tgt_info)
          call set_arg(solve_ccr12_gs,PRINT_MEL,'LIST',1,tgt_info,
     &         val_label=(/mel_ccr12lg0/))
          call set_arg(solve_ccr12_gs,PRINT_MEL,'COMMENT',1,tgt_info,
     &         val_str='>>> F12 corr. (pert) :')
          call set_arg(solve_ccr12_gs,PRINT_MEL,'FORMAT',1,tgt_info,
     &         val_str='SCAL F20.12')
          call set_dependency(solve_ccr12_gs,mel_ccr12en0def,tgt_info)
          call set_rule2(solve_ccr12_gs,ADD,tgt_info)
          call set_arg(solve_ccr12_gs,ADD,'LIST_SUM',1,tgt_info,
     &         val_label=(/mel_ccr12en0/))
          call set_arg(solve_ccr12_gs,ADD,'LISTS',2,tgt_info,
     &         val_label=(/mel_ccr12lg0,mel_ccen0/))
          call set_arg(solve_ccr12_gs,ADD,'FAC',2,tgt_info,
     &         val_rl8=(/1d0,1d0/))
          call set_rule2(solve_ccr12_gs,PRINT_MEL,tgt_info)
          call set_arg(solve_ccr12_gs,PRINT_MEL,'LIST',1,tgt_info,
     &         val_label=(/mel_ccr12en0/))
          call set_arg(solve_ccr12_gs,PRINT_MEL,'COMMENT',1,tgt_info,
     &         val_str='>>> CC-F12(pt) energy:')
          call set_arg(solve_ccr12_gs,PRINT_MEL,'FORMAT',1,tgt_info,
     &         val_str='SCAL F20.12')
          

        end if

      end if

      ncnt = is_keyword_set('calculate.check_S')
      if (ncnt.gt.0.and..not.r12fix) then
        call add_target(check_s,ttype_gen,.true.,tgt_info)
        call set_dependency(check_s,fopt_ccr12_s_v,tgt_info)
        call set_dependency(check_s,eval_r12_inter,tgt_info)

        do icnt = 1, ncnt
          call get_argument_value('calculate.check_S','sym',
     &         keycount=icnt,
     &         iarr=sym_arr)
          call get_argument_value('calculate.check_S','msc',
     &         keycount=icnt,
     &         ival=msc_s)
          do isym = 1, orb_info%nsym
            if (sym_arr(isym).eq.0) cycle

            call me_list_label(me_label,mel_evs,isym,0,0,msc_s,.false.)
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = 'SDIAG'
            labels(2) = op_evs
            call me_list_parameters(-1,parameters,
     &           msc_s,0,isym,0,0,.false.)
            call set_rule(check_s,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)
            labels(1) = 'SDIAG'
            labels(2) = mel_x_inter
            call set_rule(check_s,ttype_opme,PRECONDITIONER,
     &           labels,2,1,
     &           parameters,0,tgt_info)
            call solve_parameters(-1,parameters,2,1,sym_arr(isym),'DIA')
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = 'SDIAG'
            labels(3) = op_s_evs
            labels(4) = op_evs
            labels(5) = fopt_ccr12_s_v
            call set_rule(check_s,ttype_opme,SOLVEEVP,
     &           labels,5,1,
     &           parameters,2,tgt_info)
            call set_rule(check_s,ttype_opme,DELETE_ME_LIST,
     &           'SDIAG',1,1,
     &           parameters,0,tgt_info)

          end do
        end do

      end if


      return
      end
