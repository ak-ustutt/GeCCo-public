*----------------------------------------------------------------------*
      subroutine set_mpr12f_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in MP-R12 calculations
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'opdim.h'
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
     &     isym, ms, msc, sym_arr(8), nlabel, mode,
     &     occ_def(ngastp,2,20), ndef
      logical ::
     &     needed,r12fix
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
      character(12) ::
     &     approx

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

      ! get keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.MP','level',ival=level)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=mode)

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

      ! No residual if R12-amplitudes are fixed.
      if(.not.r12fix)then
        ! residual
        call add_target(op_omgr12,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .false.,min_rank,max_rank,0,2)
        call set_rule(op_omgr12,ttype_op,DEF_R12INTERM,
     &                op_omgr12,1,1,
     &                parameters,1,tgt_info)
      endif

      ! Residual for special T1' operator.
      if(mode.gt.0)then
        call add_target(op_omgcex,ttype_op,.false.,tgt_info)
        call set_dependency(op_omgcex,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.false.)
        call set_rule(op_omgcex,ttype_op,CLONE_OP,
     &                op_omgcex,1,1,
     &                parameters,1,tgt_info)

c      ! diagonal

        call add_target(op_diar12,ttype_op,.false.,tgt_info)
        call set_dependency(op_diar12,op_omgr12,tgt_info)
        call cloneop_parameters(-1,parameters,
     &       op_omgcex,.false.) ! <- dagger=.false.
        call set_rule(op_diar12,ttype_op,CLONE_OP,
     &              op_diar12,1,1,
     &              parameters,1,tgt_info)
      endif

      call add_target(op_bprc,ttype_op,.false.,tgt_info)
      occ_def=0
      ndef = 1
      if (mode.eq.1) then
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
      else if (mode.eq.2) then
        ndef = 1
      else if (mode.eq.3) then
        ndef = 2
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,ndef)
      call set_rule(op_bprc,ttype_op,DEF_OP_FROM_OCC,
     &              op_bprc,1,1,
     &              parameters,2,tgt_info)

      call add_target(op_xprc,ttype_op,.false.,tgt_info)
      occ_def=0
      ndef = 1
      if (mode.eq.1) then
        ndef = 1
        occ_def(IPART,1,1) = 1
        occ_def(IPART,2,1) = 1
      else if (mode.eq.2) then
        ndef = 1
      else if (mode.eq.3) then
        ndef = 2
        occ_def(IPART,1,2) = 1
        occ_def(IPART,2,2) = 1
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
      if(mode.gt.0)then
        labels(8) = op_cexbar
        labels(9) = op_cex
        nlabel = 9
      endif
      call set_dependency(form_mpr12lg0,op_mpr12lg,tgt_info)
      call set_dependency(form_mpr12lg0,op_ham,tgt_info)
      call set_dependency(form_mpr12lg0,op_r12,tgt_info)
c      call set_dependency(form_mpr12lg0,op_rba,tgt_info)
      call set_dependency(form_mpr12lg0,op_tbar,tgt_info)
      call set_dependency(form_mpr12lg0,op_top,tgt_info)

      if(mode.gt.0)then
        call set_dependency(form_mpr12lg0,op_cexbar,tgt_info)
        call set_dependency(form_mpr12lg0,op_cex,tgt_info)
      endif
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,level,'---')
      call set_rule(form_mpr12lg0,ttype_frm,DEF_MPR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)

c dbg
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

c dbg

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

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12en0
      labels(2) = form_mpr12lg0
      labels(3) = op_mpr12en
      labels(4) = op_tbar
      nlabel = 4
      if(mode.gt.0)then
        labels(5) = op_cexbar
        nlabel = 5
      endif
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

      if(mode.gt.0)then
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

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! MP ground state:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_mpr12_0
      labels(2) = form_mpr12en0
      labels(3) = form_mpr12rs_t
      ncat = 2
      if(mode.gt.0)then
        labels(4) = form_mpr12rs_cex
        ncat = 3
      endif
      nint = 0
      call add_target(fopt_mpr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12en0,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12rs_t,tgt_info)
      if(mode.gt.0)
     &     call set_dependency(fopt_mpr12_0,form_mpr12rs_cex,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_topdef,tgt_info)
      if(mode.gt.0)then
        call set_dependency(fopt_mpr12_0,mel_omgcexdef,tgt_info)
        call set_dependency(fopt_mpr12_0,mel_cex_def,tgt_info)
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

      ! lagrangian for testing:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_mpr12lg0
      labels(2) = form_mpr12lg0
      ncat = 1
      nint = 0
      call add_target(fopt_mpr12lg0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_mpr12lg0,form_mpr12lg0,tgt_info)
      call set_dependency(fopt_mpr12lg0,mel_topdef,tgt_info)
      if(mode.gt.0)then
        call set_dependency(fopt_mpr12lg0,mel_cex_def,tgt_info)
      endif
      call set_dependency(fopt_mpr12lg0,mel_ham,tgt_info)
      call set_dependency(fopt_mpr12lg0,mel_mpr12lg0def,tgt_info)      
      call set_dependency(fopt_mpr12lg0,mel_v_def,tgt_info)      
      call set_dependency(fopt_mpr12lg0,mel_b_def,tgt_info)      
      call set_dependency(fopt_mpr12lg0,mel_bh_def,tgt_info)      
      call set_dependency(fopt_mpr12lg0,mel_c_def,tgt_info)      
      call set_dependency(fopt_mpr12lg0,mel_x_def,tgt_info)
      call set_dependency(fopt_mpr12lg0,eval_r12_inter,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_mpr12lg0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

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
     &     0,0,1,0,0)
      call set_rule(mel_mpr12lg0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target(mel_mpr12en0def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12en0def,op_mpr12en,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12en0
      labels(2) = op_mpr12en
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(mel_mpr12en0def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      if(mode.gt.0)then
        call add_target(mel_cex_def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_cex_def,op_cex,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_cex
        labels(2) = op_cex
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
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
     &       0,0,1,0,0)
        call set_rule(mel_cexbar_def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

        call add_target('DIATEST',ttype_opme,.false.,tgt_info)
        call set_dependency('DIATEST',op_diar12,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = 'DIATEST'
        labels(2) = op_diar12
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
        call set_rule('DIATEST',ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
        call scale_parameters(-1,parameters,1,1,0.02d0,12)
        labels(1) = 'DIATEST'
        labels(2) = mel_dia1
        call set_rule('DIATEST',ttype_opme,SCALE,
     &       labels,2,1,
     &       parameters,1,tgt_info)
        
      endif

c      if(.not.r12fix)then
c        ! OMG-R12 list definition
c        call add_target(mel_omgr12def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_omgr12def,op_omgr12,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_omgr12
c        labels(2) = op_omgr12
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0)
c        call set_rule(mel_omgr12def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

      if(mode.gt.0)then
        ! OMG-R12-EXT list definition
        call add_target(mel_omgcexdef,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_omgcexdef,op_omgcex,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_omgcex
        labels(2) = op_omgcex
        call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
        call set_rule(mel_omgcexdef,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      endif

      call add_target(me_bprc,ttype_opme,.false.,tgt_info)
      call set_dependency(me_bprc,op_bprc,tgt_info)
      call set_dependency(me_bprc,eval_r12_inter,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = me_bprc
      labels(2) = op_bprc
      call me_list_parameters(-1,parameters,
     &       0,0,1,0,0)
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
     &       0,0,1,0,0)
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

      if(mode.gt.0)then
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
        
        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
c        call set_dependency(solve_mpr12_gs,'DIATEST',tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_b_inv,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_b_dia,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_x_inv,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
        call set_dependency(solve_mpr12_gs,me_bprc,tgt_info)
        call set_dependency(solve_mpr12_gs,me_xprc,tgt_info)
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_cex
        labels(3) = mel_omg
        labels(4) = mel_omgcex
        labels(5) = mel_dia1
        labels(6) = mel_dia1
        labels(7) = mel_mpr12en0
        labels(8) = fopt_mpr12_0
        labels(9) = me_bprc
        labels(10) = me_xprc
        labels(11) = mel_ham
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)
c testing
c        labels(1) = fopt_mpr12lg0
c        call set_dependency(solve_mpr12_gs,fopt_mpr12lg0,tgt_info)
c        call set_rule(solve_mpr12_gs,ttype_opme,EVAL,
c     &       labels,1,1,
c     &       parameters,2,tgt_info)

      else
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
        
        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
c        if(.not.mode)then
          call solve_parameters(-1,parameters,2, 1,1,'DIA/BLK')
c        else
c          call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
c        endif
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_omg
        labels(3) = mel_dia1
        labels(4) = mel_mpr12en0
        labels(5) = fopt_mpr12_0
        labels(6) = mel_ham
        nlabel = 6
c        if(mode.gt.0)then
c          labels(7) = mel_cex
c          labels(8) = mel_cexbar
c          labels(9) = mel_omgcex
c          nlabel = 9
c        endif
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,nlabel,2,
     &       parameters,2,tgt_info)

      endif

      return
      end
