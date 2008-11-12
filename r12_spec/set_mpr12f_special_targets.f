*----------------------------------------------------------------------*
      subroutine set_mpr12f_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in MP-R12 calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
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
     &     min_rank, max_rank, level,
     &     isim, ncat, nint, icnt, ansatz,
     &     isym, ms, msc, sym_arr(8), nlabel, extend,
     &     occ_def(ngastp,2,20), ndef, r12op
      logical ::
     &     needed,r12fix, set_tp, set_tpp
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
      character(12) ::
     &     approx
      character(6), parameter ::
     &     op_jac = 'JACOBI',
     &     form_jac = 'FRMJAC',
     &     fopt_jac = 'OPTJAC',
     &     mel_jac = 'MELJAC',
     &     mel_jac_def = 'DEFJAC'

      character ::
     &     op_bprc*4 = 'BPRC',
     &     op_xprc*4 = 'XPRC',
     &     me_bprc*8 = 'BPRCLIST',
     &     me_xprc*8 = 'XPRCLIST',
     &     medef_bprc*12 = 'DEF-BPRCLIST',
     &     medef_xprc*12 = 'DEF-XPRCLIST',
c     &     mel_mpr12lg0def*8 = 'L(MPR12)',
     &     fopt_mpr12lg0*13 = 'LAG_MPR12_OPT'

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for MP-R12 ...'

      msc = 0

      ! get keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.MP','level',ival=level)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)

      set_tp = extend.gt.0.or.r12op.eq.1.or.r12op.ge.3
      set_tpp = r12op.eq.2.or.r12op.ge.3

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Lagrange functional
      call add_target(op_mpr12lg,ttype_op,.false.,tgt_info)
      call set_rule(op_mpr12lg,ttype_op,DEF_SCALAR,
     &              op_mpr12lg,1,1,
     &              parameters,0,tgt_info)
      
      ! Energy
      call add_target(op_mpr12en,ttype_op,.false.,tgt_info)
      call set_rule(op_mpr12en,ttype_op,DEF_SCALAR,
     &              op_mpr12en,1,1,
     &              parameters,0,tgt_info)

      ! Residual for special T1' operator.
      if (set_tp) then
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

c        call add_target(op_jac,ttype_op,.false.,tgt_info)
c        occ_def = 0
c        occ_def(IHOLE,1,1) = 1
c        occ_def(IPART,1,1) = 1
c        occ_def(IHOLE,2,1) = 1
c        occ_def(IPART,2,1) = 1
c
c        occ_def(IPART,1,2) = 1
c        occ_def(IPART,2,2) = 1
c        call op_from_occ_parameters(-1,parameters,2,
c    &        occ_def,2,1,2)
c        call set_rule(op_jac,ttype_op,DEF_OP_FROM_OCC,
c    &        op_jac,1,1,
c    &        parameters,2,tgt_info)

      ! diagonal

c        call add_target(op_diar12,ttype_op,.false.,tgt_info)
c        call set_dependency(op_diar12,op_omgcex,tgt_info)
c        call cloneop_parameters(-1,parameters,
c     &       op_omgcex,.false.) ! <- dagger=.false.
c        call set_rule(op_diar12,ttype_op,CLONE_OP,
c     &              op_diar12,1,1,
c     &              parameters,1,tgt_info)
c      endif

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
     &     occ_def,ndef,1,ndef)
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
     &     occ_def,ndef,1,ndef)
      call set_rule(op_xprc,ttype_op,DEF_OP_FROM_OCC,
     &              op_xprc,1,1,
     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      call add_target(form_mpr12lg0,ttype_frm,.false.,tgt_info)
      ! (a) set formal Lagrangian (in 'complete' basis)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12lg0
      labels(2) = op_mpr12lg
      labels(3) = op_ham
      labels(4) = op_r12
      labels(5) = op_r12
      labels(6) = op_tbar
      labels(7) = op_top
      nlabel = 7
      if (set_tp) then
        call set_dependency(form_mpr12lg0,op_cex,tgt_info)
        call set_dependency(form_mpr12lg0,op_cexbar,tgt_info)
        labels(nlabel+1) = op_cexbar
        labels(nlabel+2) = op_cex
        nlabel = nlabel+2
      end if
      if (set_tpp) then
        call set_dependency(form_mpr12lg0,op_cexx,tgt_info)
        call set_dependency(form_mpr12lg0,op_cexxbar,tgt_info)
        labels(nlabel+1) = op_cexxbar
        labels(nlabel+2) = op_cexx
        nlabel = nlabel+2
      end if
      call set_dependency(form_mpr12lg0,op_mpr12lg,tgt_info)
      call set_dependency(form_mpr12lg0,op_ham,tgt_info)
      call set_dependency(form_mpr12lg0,op_r12,tgt_info)
c      call set_dependency(form_mpr12lg0,op_rba,tgt_info)
      call set_dependency(form_mpr12lg0,op_tbar,tgt_info)
      call set_dependency(form_mpr12lg0,op_top,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,level,'---')
      call set_rule(form_mpr12lg0,ttype_frm,DEF_MPR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12lg0 ! output formula (itself)
      labels(2) = form_mpr12lg0 ! input formula
      labels(3) = form_r12_vint ! the intermediates to be factored
      labels(4) = trim(form_r12_vint)//'^+'
      labels(5) = form_r12_bint
      labels(6) = form_r12_bhint
      labels(7) = form_r12_xint
      nint = 5
      call set_dependency(form_mpr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_bint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_bhint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_xint,tgt_info)
      if (ansatz.ne.1) then
        labels(2+nint+1) = form_r12_cint
        labels(2+nint+2) = trim(form_r12_cint)//'^+'
        call set_dependency(form_mpr12lg0,form_r12_cint,tgt_info)
        nint = nint+2
      end if
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,nint,'---')
      call set_rule(form_mpr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)
      ! (c) post-processing: remove terms which do not contribute for
      !     the given R12-approximation
      ! .... to come
      ! fix for A
      if (trim(approx).eq.'A') then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12lg0
        labels(2) = form_mpr12lg0
        labels(3) = op_mpr12lg
        labels(4) = op_x_inter
        call set_rule(form_mpr12lg0,ttype_frm,INVARIANT,
     &              labels,4,1,
     &              title_mpr12lg0,1,tgt_info)
      end if

c test
c      if (extend.eq.6) then
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = form_mpr12lg0
c        labels(2) = form_mpr12lg0
c        labels(3) = op_mpr12lg
c        labels(4) = op_cex
c        labels(5) = op_cexbar
c        call set_rule(form_mpr12lg0,ttype_frm,INVARIANT,
c     &              labels,5,1,
c     &              title_mpr12lg0,1,tgt_info)
c      end if
c test

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12en0
      labels(2) = form_mpr12lg0
      labels(3) = op_mpr12en
      labels(4) = op_tbar
      nlabel = 4
      if (set_tp) then
        labels(nlabel+1) = op_cexbar
        nlabel = nlabel+1
      end if
      if (set_tpp) then
        labels(nlabel+1) = op_cexxbar
        nlabel = nlabel+1
      end if
      call add_target(form_mpr12en0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_mpr12en0,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12en0,op_mpr12en,tgt_info)
      call set_rule(form_mpr12en0,ttype_frm,INVARIANT,
     &              labels,nlabel,1,
     &              title_mpr12en0,1,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12rs_t
      labels(2) = form_mpr12lg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_mpr12rs_t,ttype_frm,.false.,tgt_info)
      call set_dependency(form_mpr12rs_t,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12rs_t,op_omg,tgt_info)
      call set_rule(form_mpr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_mpr12rs_t,1,tgt_info)

      if(set_tp)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12rs_cex
        labels(2) = form_mpr12lg0
        labels(3) = op_omgcex
        labels(4) = op_cexbar
        labels(5) = ' '
        call add_target(form_mpr12rs_cex,ttype_frm,.false.,tgt_info)
        call set_dependency(form_mpr12rs_cex,form_mpr12lg0,tgt_info)
        call set_dependency(form_mpr12rs_cex,op_omgcex,tgt_info)
        call set_rule(form_mpr12rs_cex,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_mpr12rs_cex,1,tgt_info)
      endif

      if(set_tpp)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12rs_cexx
        labels(2) = form_mpr12lg0
        labels(3) = op_omgcexx
        labels(4) = op_cexxbar
        labels(5) = ' '
        call add_target(form_mpr12rs_cexx,ttype_frm,.false.,tgt_info)
        call set_dependency(form_mpr12rs_cexx,form_mpr12lg0,tgt_info)
        call set_dependency(form_mpr12rs_cexx,op_omgcexx,tgt_info)
        call set_rule(form_mpr12rs_cexx,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_mpr12rs_cexx,1,tgt_info)
      endif

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! MP ground state:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_mpr12_0
      labels(2) = form_mpr12en0
      labels(3) = form_mpr12rs_t
      ncat = 2
      if(set_tp.and..not.r12fix)then
        labels(1+ncat+1) = form_mpr12rs_cex
        ncat = ncat+1
      endif
      if(set_tpp.and..not.r12fix)then
        labels(1+ncat+1) = form_mpr12rs_cexx
        ncat = ncat+1
      endif
      nint = 0
      call add_target(fopt_mpr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12en0,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12rs_t,tgt_info)
      if (set_tp.and..not.r12fix)
     &     call set_dependency(fopt_mpr12_0,form_mpr12rs_cex,tgt_info)
      if (set_tpp.and..not.r12fix)
     &     call set_dependency(fopt_mpr12_0,form_mpr12rs_cexx,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_topdef,tgt_info)
      if(set_tp.and..not.r12fix)then
        call set_dependency(fopt_mpr12_0,mel_omgcexdef,tgt_info)
        call set_dependency(fopt_mpr12_0,mel_cex_def,tgt_info)
      endif
      if(set_tpp.and..not.r12fix)then
        call set_dependency(fopt_mpr12_0,mel_omgcexxdef,tgt_info)
        call set_dependency(fopt_mpr12_0,mel_cexx_def,tgt_info)
      endif
      call set_dependency(fopt_mpr12_0,mel_ham,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_mpr12en0def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_v_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_b_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_bh_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_c_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_x_def,tgt_info)
      call set_dependency(fopt_mpr12_0,eval_r12_inter,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_mpr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      if(extend.gt.0)then
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = fopt_jac
c        labels(2) = form_jac
c        ncat = 1
c        nint = 0
c        call add_target(fopt_jac,ttype_frm,.true.,tgt_info)
c        call set_dependency(fopt_jac,form_jac,tgt_info)
c        call set_dependency(fopt_jac,mel_jac_def,tgt_info)
c        call set_dependency(fopt_jac,mel_omgcexdef,tgt_info)
c        call set_dependency(fopt_jac,mel_cex_def,tgt_info)
c        call set_dependency(fopt_jac,mel_ham,tgt_info)
c        call set_dependency(fopt_jac,mel_mpr12en0def,tgt_info)      
c        call set_dependency(fopt_jac,mel_v_def,tgt_info)      
c        call set_dependency(fopt_jac,mel_b_def,tgt_info)      
c        call set_dependency(fopt_jac,mel_bh_def,tgt_info)      
c        call set_dependency(fopt_jac,mel_c_def,tgt_info)      
c        call set_dependency(fopt_jac,mel_x_def,tgt_info)
c        call set_dependency(fopt_jac,eval_r12_inter,tgt_info)
c        call opt_parameters(-1,parameters,ncat,nint)
c        call set_rule(fopt_jac,ttype_frm,OPTIMIZE,
c     &       labels,ncat+nint+1,1,
c     &       parameters,1,tgt_info)
c
c      ! lagrangian for testing:
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = fopt_mpr12lg0
c      labels(2) = form_mpr12lg0
c      ncat = 1
c      nint = 0
c      call add_target(fopt_mpr12lg0,ttype_frm,.false.,tgt_info)
c      call set_dependency(fopt_mpr12lg0,form_mpr12lg0,tgt_info)
c      call set_dependency(fopt_mpr12lg0,mel_topdef,tgt_info)
c      if(extend.gt.0)then
c        call set_dependency(fopt_mpr12lg0,mel_cex_def,tgt_info)
c      endif
c      call set_dependency(fopt_mpr12lg0,mel_ham,tgt_info)
c      call set_dependency(fopt_mpr12lg0,mel_mpr12lg0def,tgt_info)      
c      call set_dependency(fopt_mpr12lg0,mel_v_def,tgt_info)      
c      call set_dependency(fopt_mpr12lg0,mel_b_def,tgt_info)      
c      call set_dependency(fopt_mpr12lg0,mel_bh_def,tgt_info)      
c      call set_dependency(fopt_mpr12lg0,mel_c_def,tgt_info)      
c      call set_dependency(fopt_mpr12lg0,mel_x_def,tgt_info)
c      call set_dependency(fopt_mpr12lg0,eval_r12_inter,tgt_info)
c      call opt_parameters(-1,parameters,ncat,nint)
c      call set_rule(fopt_mpr12lg0,ttype_frm,OPTIMIZE,
c     &              labels,ncat+nint+1,1,
c     &              parameters,1,tgt_info)
c      endif
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! L0/E0:
      call add_target(mel_mpr12lg0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12lg0def,op_mpr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12lg0
      labels(2) = op_mpr12lg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_mpr12lg0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_mpr12en0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12en0def,op_mpr12en,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12en0
      labels(2) = op_mpr12en
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_mpr12en0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      if (set_tp) then
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

        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

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

c        call add_target(mel_jac_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_jac_def,op_jac,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_jac
c        labels(2) = op_jac
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0,.false.)
c        call set_rule(mel_jac_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)

        call add_target('DIATEST',ttype_opme,.false.,tgt_info)
        call set_dependency('DIATEST',op_diar12,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = 'DIATEST'
        labels(2) = op_diar12
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0,.false.)
        call set_rule('DIATEST',ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        call scale_parameters(-1,parameters,1,1,1,0.02d0,12)
        labels(1) = 'DIATEST'
        labels(2) = mel_dia1
        call set_rule('DIATEST',ttype_opme,SCALE,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        
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

c      if(.not.r12fix)then
c        ! OMG-R12 list definition
c        call add_target(mel_omgr12def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_omgr12def,op_omgr12,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_omgr12
c        labels(2) = op_omgr12
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0,.false.)
c        call set_rule(mel_omgr12def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

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
      ! totally symmetric dia for use below:
      if (set_tp.and..not.set_tpp.and..not.r12fix) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_diar12,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call set_dependency(solve_mpr12_gs,eval_r12_inter,tgt_info)
        call set_dependency(solve_mpr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_mpr12_gs,me_xprc,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cex
        labels(3) = mel_omg
        labels(4) = mel_omgcex
        labels(5) = mel_dia1    
        labels(6) = mel_dia1  !mel_diar12  
        labels(7) = mel_mpr12en0
        labels(8) = fopt_mpr12_0
        labels(9) = me_bprc
        labels(10)= me_xprc
        labels(11) = mel_ham
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)
        
      else if (.not.set_tp.and.set_tpp.and..not.r12fix) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call set_dependency(solve_mpr12_gs,eval_r12_inter,tgt_info)
        call set_dependency(solve_mpr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_mpr12_gs,me_xprc,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cexx
        labels(3) = mel_omg
        labels(4) = mel_omgcexx
        labels(5) = mel_dia1    ! dummy
        labels(6) = mel_dia1    ! dummy
        labels(7) = mel_mpr12en0
        labels(8) = fopt_mpr12_0
        labels(9) = me_bprc
        labels(10)= me_xprc
        labels(11) = mel_ham
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)

      else if (set_tp.and.set_tpp.and..not.r12fix) then
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call set_dependency(solve_mpr12_gs,eval_r12_inter,tgt_info)
        call set_dependency(solve_mpr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_mpr12_gs,me_xprc,tgt_info)
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
        labels(10) = mel_mpr12en0
        labels(11) = fopt_mpr12_0
        labels(12) = me_bprc
        labels(13)= me_xprc
        labels(14) = mel_ham
        call solve_parameters(-1,parameters,2, 3,1,'DIA/BLK/BLK')
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,14,6,
     &       parameters,2,tgt_info)
        
      else
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call set_dependency(solve_mpr12_gs,eval_r12_inter,tgt_info)
        call solve_parameters(-1,parameters,2, 1,1,'DIA')
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_omg
        labels(3) = mel_dia1
        labels(4) = mel_mpr12en0
        labels(5) = fopt_mpr12_0
        labels(6) = mel_ham
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,6,2,
     &       parameters,2,tgt_info)
      end if


c      if(set_tp.gt.0.and.set_tpp.and..not.r12fix)then
c        ! totally symmetric dia for use below:
c        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
c        
c        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
c
cc        labels(1:10)(1:len_target_name) = ' '
cc        labels(1) = fopt_jac
cc        call set_rule(solve_mpr12_gs,ttype_opme,EVAL,
cc     &       labels,1,0,
cc     &       parameters,0,tgt_info)
c
c        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
cc        call set_dependency(solve_mpr12_gs,'DIATEST',tgt_info)
cc        call set_dependency(solve_mpr12_gs,mel_b_inv,tgt_info)
cc        call set_dependency(solve_mpr12_gs,mel_b_dia,tgt_info)
cc        call set_dependency(solve_mpr12_gs,mel_x_inv,tgt_info)
c        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
c        call set_dependency(solve_mpr12_gs,me_bprc,tgt_info)
c        call set_dependency(solve_mpr12_gs,me_xprc,tgt_info)
cc        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_top
c        labels(2) = mel_cex
c        labels(3) = mel_omg
c        labels(4) = mel_omgcex
c        labels(5) = mel_dia1
c        labels(6) = mel_dia1
c        labels(7) = mel_mpr12en0
c        labels(8) = fopt_mpr12_0
c        labels(9) = me_bprc
c        labels(10) = me_xprc
c        labels(11) = mel_ham
c        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
c     &       labels,11,4,
c     &       parameters,2,tgt_info)
cc testing
cc        labels(1) = fopt_mpr12lg0
cc        call set_dependency(solve_mpr12_gs,fopt_mpr12lg0,tgt_info)
cc        call set_rule(solve_mpr12_gs,ttype_opme,EVAL,
cc     &       labels,1,1,
cc     &       parameters,2,tgt_info)
c
c      else
c        ! totally symmetric dia for use below:
c        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
c        
c        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
c        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
cc        if(.not.extend)then
c          call solve_parameters(-1,parameters,2, 1,1,'DIA/BLK')
cc        else
cc          call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
cc        endif
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_top
c        labels(2) = mel_omg
c        labels(3) = mel_dia1
c        labels(4) = mel_mpr12en0
c        labels(5) = fopt_mpr12_0
c        labels(6) = mel_ham
c        nlabel = 6
cc        if(extend.gt.0)then
cc          labels(7) = mel_cex
cc          labels(8) = mel_cexbar
cc          labels(9) = mel_omgcex
cc          nlabel = 9
cc        endif
c        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
c     &       labels,nlabel,2,
c     &       parameters,2,tgt_info)
c
c      endif

      return
      end
