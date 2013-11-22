*----------------------------------------------------------------------*
      subroutine set_mpr12_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in MP-R12 calculations
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
     &     min_rank, max_rank, level,
     &     isim, ncat, nint, icnt, ansatz,
     &     isym, ms, msc, sym_arr(8), nlabel, mode
      logical ::
     &     needed,r12fix,extend
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
      character(12) ::
     &     approx

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting special targets for MP-R12 ...'

      msc = +1    ! assuming closed shell

      ! get keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.MP','level',ival=level)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      approx(1:12) = ' '
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=mode)
      extend = mode.gt.0

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
      if(extend)then
        call add_target(op_omgcex,ttype_op,.false.,tgt_info)
        call set_dependency(op_omgcex,op_cex,tgt_info)
        call cloneop_parameters(-1,parameters,
     &                          op_cex,.false.)
        call set_rule(op_omgcex,ttype_op,CLONE_OP,
     &                op_omgcex,1,1,
     &                parameters,1,tgt_info)
      endif

c      ! diagonal
c      call add_target(op_diar12,ttype_op,.false.,tgt_info)
c      call set_dependency(op_diar12,op_omgr12,tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &                        op_omgr12,.false.) ! <- dagger=.false.
c      call set_rule(op_diar12,ttype_op,CLONE_OP,
c     &              op_diar12,1,1,
c     &              parameters,1,tgt_info)

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
c      nlabel = 7
c      if(.not.r12fix)then
      labels(8) = op_cba
      labels(9) = op_c12
      nlabel = 9
c      endif
      if(extend)then
        labels(10) = op_cexbar
        labels(11) = op_cex
        nlabel = 11
      endif
      call set_dependency(form_mpr12lg0,op_mpr12lg,tgt_info)
      call set_dependency(form_mpr12lg0,op_ham,tgt_info)
      call set_dependency(form_mpr12lg0,op_r12,tgt_info)
c      call set_dependency(form_mpr12lg0,op_rba,tgt_info)
      call set_dependency(form_mpr12lg0,op_tbar,tgt_info)
      call set_dependency(form_mpr12lg0,op_top,tgt_info)

c      if(.not.r12fix)then
      call set_dependency(form_mpr12lg0,op_cba,tgt_info)
      call set_dependency(form_mpr12lg0,op_c12,tgt_info)
c      endif
      if(extend)then
        call set_dependency(form_mpr12lg0,op_cexbar,tgt_info)
        call set_dependency(form_mpr12lg0,op_cex,tgt_info)
      endif
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,level,'---')
      call set_rule(form_mpr12lg0,ttype_frm,DEF_MPR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)

c dbg
      if(.not.extend)then
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = form_mpr12lg0 ! output formula (itself)
      labels(2) = form_mpr12lg0 ! input formula
c      if(.not.r12fix)then
      labels(3) = form_r12_vint ! the intermediates to be factored
      labels(4) = trim(form_r12_vint)//'^+'
      labels(5) = form_r12_xint
      labels(6) = form_r12_bint
      nint = 4
c      else
c        labels(3) = form_r12_v0int ! the intermediates to be factored
c        labels(4) = trim(form_r12_v0int)//'^+'
c        labels(5) = form_r12_x1int
c        labels(6) = form_r12_b0int
cc dbg
cc        labels(6) = form_r12_xint
cc dbg
c        nint = 4
c      endif
c      if(.not.r12fix)then
      call set_dependency(form_mpr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_mpr12lg0,form_r12_bint,tgt_info)
c      else
c        call set_dependency(form_mpr12lg0,form_r12_v0int,tgt_info)
c        call set_dependency(form_mpr12lg0,form_r12_x1int,tgt_info)
c        call set_dependency(form_mpr12lg0,form_r12_b0int,tgt_info)
c      endif
      if (ansatz.ne.1) then
c        if(.not.r12fix)then
          labels(7) = form_r12_cint
          labels(8) = trim(form_r12_cint)//'^+'
          call set_dependency(form_mpr12lg0,form_r12_cint,tgt_info)
          nint = 6
c        else
c          labels(7) = form_r12_cint
c          labels(8) = trim(form_r12_cint)//'^+'
c          call set_dependency(form_mpr12lg0,form_r12_cint,tgt_info)
c          nint = 6
c        endif
      end if
      call form_parameters(-1,
     &     parameters,2,title_mpr12lg0,nint,'---')
      call set_rule(form_mpr12lg0,ttype_frm,FACTOR_OUT,
     &              labels,nint+2,1,
     &              parameters,2,tgt_info)

      endif
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
      if(.not.r12fix)then
        labels(5) = op_cba
        nlabel = 5
      else
        if(extend)then
          labels(5) = op_cexbar
          nlabel = 5
        endif
      endif
      call add_target(form_mpr12en0,ttype_frm,.true.,tgt_info)
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
      call add_target(form_mpr12rs_t,ttype_frm,.true.,tgt_info)
      call set_dependency(form_mpr12rs_t,form_mpr12lg0,tgt_info)
      call set_dependency(form_mpr12rs_t,op_omg,tgt_info)
      call set_rule(form_mpr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_mpr12rs_t,1,tgt_info)

      if(.not.r12fix)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12rs_c
        labels(2) = form_mpr12lg0
        labels(3) = op_omgr12
        labels(4) = op_cba
        labels(5) = ' '
        call add_target(form_mpr12rs_c,ttype_frm,.true.,tgt_info)
        call set_dependency(form_mpr12rs_c,form_mpr12lg0,tgt_info)
        call set_dependency(form_mpr12rs_c,op_omgr12,tgt_info)
        call set_rule(form_mpr12rs_c,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_mpr12rs_c,1,tgt_info)
      endif

      if(extend)then
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = form_mpr12rs_cex
        labels(2) = form_mpr12lg0
        labels(3) = op_omgcex
        labels(4) = op_cexbar
        labels(5) = ' '
        call add_target(form_mpr12rs_cex,ttype_frm,.true.,tgt_info)
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
      if(.not.r12fix)then
        labels(4) = form_mpr12rs_c
        ncat = 3
      endif
c      if(extend)then
c        labels(4) = form_mpr12rs_cex
c        ncat = 3
c      endif
      nint = 0
      call add_target(fopt_mpr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12en0,tgt_info)
      call set_dependency(fopt_mpr12_0,form_mpr12rs_t,tgt_info)
      if(.not.r12fix)
     &     call set_dependency(fopt_mpr12_0,form_mpr12rs_c,tgt_info)
c      if(extend)
c     &     call set_dependency(fopt_mpr12_0,form_mpr12rs_cex,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_topdef,tgt_info)
      if(.not.r12fix)
     &     call set_dependency(fopt_mpr12_0,mel_omgr12def,tgt_info)
c      if(extend)then
c        call set_dependency(fopt_mpr12_0,mel_omgcexdef,tgt_info)
c        call set_dependency(fopt_mpr12_0,mel_cex_def,tgt_info)
c      endif
      call set_dependency(fopt_mpr12_0,mel_c12def,tgt_info)
      if(r12fix)then
        call set_dependency(fopt_mpr12_0,mel_cbardef,tgt_info)
      endif
      call set_dependency(fopt_mpr12_0,mel_ham,tgt_info)
      call set_dependency(fopt_mpr12_0,mel_mpr12en0def,tgt_info)      
c      if(.not.r12fix)then
      call set_dependency(fopt_mpr12_0,mel_v_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_b_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_c_def,tgt_info)      
      call set_dependency(fopt_mpr12_0,mel_x_def,tgt_info)
c      else
c        call set_dependency(fopt_mpr12_0,mel_v0_def,tgt_info)      
c        call set_dependency(fopt_mpr12_0,mel_b0_def,tgt_info)
cc dbg
cc        call set_dependency(fopt_mpr12_0,mel_x_def,tgt_info)
cc dbg
c        call set_dependency(fopt_mpr12_0,mel_x1_def,tgt_info)
c      endif
      call set_dependency(fopt_mpr12_0,eval_r12_inter,tgt_info)
c dbg
c      if(r12fix)
c     &     call set_dependency(fopt_mpr12_0,'EVALINTS',tgt_info)
c dbg
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_mpr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)


*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! L0/E0:
      call add_target(mel_mpr12lg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_mpr12lg0,op_mpr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_mpr12lg0
      labels(2) = op_mpr12lg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_mpr12lg0,ttype_opme,DEF_ME_LIST,
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

c      if(.not.r12fix)then
      ! C12 list definition
      call add_target(mel_c12def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_c12def,op_c12,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_c12
      labels(2) = op_c12
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_c12def,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      
      if(r12fix)then
        ! Add unity.
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_c12
c     call add_parameters(-1,parameters,1,1d0,1)
        call set_rule(mel_c12def,ttype_opme,UNITY,
     &       labels,1,1,
     &       parameters,1,tgt_info)
      endif

      ! CBAR list definition
      call add_target(mel_cbardef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_cbardef,op_cba,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_cbar
      labels(2) = op_cba
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_cbardef,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      if(r12fix)then
        ! Add unity.
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_cbar
c     call add_parameters(-1,parameters,1,1d0,1)
        call set_rule(mel_cbardef,ttype_opme,UNITY,
     &       labels,1,1,
     &       parameters,1,tgt_info)
      endif

c      if(extend)then
c        call add_target(mel_cex_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_cex_def,op_cex,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_cex
c        labels(2) = op_cex
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0,.false.)
c        call set_rule(mel_cex_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c
c        call add_target(mel_cexbar_def,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_cexbar_def,op_cexbar,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_cexbar
c        labels(2) = op_cexbar
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0,.false.)
c        call set_rule(mel_cexbar_def,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif

      if(.not.r12fix)then
        ! OMG-R12 list definition
        call add_target(mel_omgr12def,ttype_opme,.false.,tgt_info)
        call set_dependency(mel_omgr12def,op_omgr12,tgt_info)
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_omgr12
        labels(2) = op_omgr12
        call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
        call set_rule(mel_omgr12def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)
      endif

c      if(extend)then
c        ! OMG-R12-EXT list definition
c        call add_target(mel_omgcexdef,ttype_opme,.false.,tgt_info)
c        call set_dependency(mel_omgcexdef,op_omgcex,tgt_info)
c        labels(1:20)(1:len_target_name) = ' '
c        labels(1) = mel_omgcex
c        labels(2) = op_omgcex
c        call me_list_parameters(-1,parameters,
c     &       0,0,1,0,0,.false.)
c        call set_rule(mel_omgcexdef,ttype_opme,DEF_ME_LIST,
c     &                labels,2,1,
c     &                parameters,1,tgt_info)
c      endif
        

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      if(.not.r12fix)then
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
        
        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_b_inv,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_b_dia,tgt_info)
c        call set_dependency(solve_mpr12_gs,mel_x_inv,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
        call solve_parameters(-1,parameters,2, 2,1,'DIA/BLK')
c        call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
        labels(1:20)(1:len_target_name) = ' '
        labels(1) = mel_top
        labels(2) = mel_c12
        labels(3) = mel_omg
        labels(4) = mel_omgr12
        labels(5) = mel_dia1
        labels(6) = mel_dia1    ! dummy
c        labels(6) = mel_b_dia
        labels(7) = mel_mpr12en0
        labels(8) = fopt_mpr12_0
        if(trim(approx).eq.'A')then
          labels(9) = mel_b_inv   ! or mel_b_inter
        else
          labels(9) = mel_b_inter
        endif
        labels(10) = mel_x_inter
        labels(11) = mel_ham
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)

      else
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)
        
        call add_target(solve_mpr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_mpr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_mpr12_gs,fopt_mpr12_0,tgt_info)
c        if(.not.extend)then
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
c        if(extend)then
c          labels(7) = mel_cex
c          labels(8) = mel_cexbar
c          labels(9) = mel_omgcex
c          nlabel = 9
c        endif
c        labels(7) = mel_c12
c        labels(8) = mel_cbar
        call set_rule(solve_mpr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,nlabel,2,
     &       parameters,2,tgt_info)

      endif

      return
      end
