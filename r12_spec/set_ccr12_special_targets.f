*----------------------------------------------------------------------*
      subroutine set_ccr12_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in CC-R12 calculations
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
     &     min_rank, max_rank, ansatz,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8), nlabel, trunc_type
      logical ::
     &     needed, r12fix, truncate, pz_eval, screen
      character(8) ::
     &     approx
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)
c dbg
      character(8), parameter ::
     &     op_omgr12_dum = 'OMG_RDUM'
c dbg

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting special targets for CC-R12 ...'
      
      msc = +1  ! assuming closed shell

      approx = '        '
      ! read keyword values
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','pz_eval',lval=pz_eval)
      call get_argument_value('method.R12','screen',lval=screen)
      call get_argument_value('method.R12','trunc',ival=trunc_type)
      truncate = trunc_type.ge.0

      if (is_keyword_set('method.truncate').gt.0) then
        call warn('set_ccr12_special_targets',
     &       '"method truncate" is obsolete, use "method R12 truncate"')
        truncate = is_keyword_set('method.truncate').gt.0
        call get_argument_value('method.truncate','trunc_type',
     &       ival=trunc_type)
      end if

c dbg
      ! for debugging only:
      r12fix = .false.
c dbg

      call get_argument_value('calculate.routes','simtraf',ival=isim)

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
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

      ! QUICK FIX: Set up dummy residual then clone its adjoint.
      ! dummy residual
c      call add_target(op_omgr12_dum,ttype_op,.false.,tgt_info)

      if(.not.r12fix)then
        call add_target(op_omgr12,ttype_op,.false.,tgt_info)
        call xop_parameters(-1,parameters,
     &       .true.,min_rank,max_rank,0,2)
c      call set_rule(op_omgr12_dum,ttype_op,DEF_R12INTERM,
c     &              op_omgr12_dum,1,1,
c     &              parameters,1,tgt_info)
        call set_rule(op_omgr12,ttype_op,DEF_R12INTERM,
     &                op_omgr12,1,1,
     &                parameters,1,tgt_info)
      endif
      
c      ! actual residual
c      call add_target(op_omgr12,ttype_op,.false.,tgt_info)
c      call set_dependency(op_omgr12,op_omgr12_dum,tgt_info)
c      call cloneop_parameters(-1,parameters,op_omgr12_dum,.true.)
c      call set_rule(op_omgr12,ttype_op,CLONE_OP,
c     &              op_omgr12,1,1,
c     &              parameters,1,tgt_info)

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
      labels(8) = op_cba
      labels(9) = op_c12
      nlabel = 9
      call set_dependency(form_ccr12lg0,op_ccr12lg,tgt_info)
      call set_dependency(form_ccr12lg0,op_ham,tgt_info)
      call set_dependency(form_ccr12lg0,op_r12,tgt_info)
c      call set_dependency(form_ccr12lg0,op_rba,tgt_info)
      call set_dependency(form_ccr12lg0,op_tbar,tgt_info)
      call set_dependency(form_ccr12lg0,op_top,tgt_info)
      call set_dependency(form_ccr12lg0,op_cba,tgt_info)
      call set_dependency(form_ccr12lg0,op_c12,tgt_info)
      call form_parameters(-1,
     &     parameters,2,title_ccr12lg0,ansatz,'---')
      call set_rule(form_ccr12lg0,ttype_frm,DEF_CCR12_LAGRANGIAN,
     &              labels,nlabel,1,
     &              parameters,2,tgt_info)
      ! (b) Factor out the R12 intermediates 
      ! (effectively removing all reference to the complete basis)
      labels(1:15)(1:len_target_name) = ' '
      labels(1) = form_ccr12lg0 ! output formula (itself)
      labels(2) = form_ccr12lg0 ! input formula
c      labels(3) = form_r12_z4int
c      labels(4) = form_r12_zint
      nint = 0
      if(.not.truncate.or.(truncate.and.trunc_type.gt.0))then
        labels(nint+3) = form_r12_pint
        call set_dependency(form_ccr12lg0,form_r12_pint,tgt_info)
        nint = nint+1
      endif
      if(.not.truncate.or.(truncate.and.trunc_type.eq.1))then
        labels(nint+3) = form_r12_zint
        call set_dependency(form_ccr12lg0,form_r12_zint,tgt_info)
        nint = nint+1
c dbg
c        labels(nint+3) = form_z_test
c        call set_dependency(form_ccr12lg0,form_z_test,tgt_info)
c        nint = nint+1
c dbg        
      endif
      if (ansatz.ne.1) then
c        if(.not.truncate)then
c          labels(nint+3) = form_r12_p3gint
c          nint = nint + 1
c          call set_dependency(form_ccr12lg0,form_r12_p3gint,tgt_info)
c        endif
        labels(nint+3) = form_r12_cint
        labels(nint+4) = trim(form_r12_cint)//'^+'
        call set_dependency(form_ccr12lg0,form_r12_cint,tgt_info)
        nint = nint + 2
      end if
      labels(nint+3) = form_r12_vint    ! the intermediates to be factored
      labels(nint+4) = form_r12_vint//'^+'
      labels(nint+5) = form_r12_bint
      labels(nint+6) = form_r12_xint
      nint = nint + 4
      call set_dependency(form_ccr12lg0,form_r12_vint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_xint,tgt_info)
      call set_dependency(form_ccr12lg0,form_r12_bint,tgt_info)
c      call set_dependency(form_ccr12lg0,form_r12_z4int,tgt_info)
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
      call set_rule(form_ccr12lg0,ttype_frm,TEX_FORMULA,
     &              labels,5,1,
     &              'ccr12_lag.tex',1,tgt_info)
      
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12en0
      labels(2) = form_ccr12lg0
      labels(3) = op_ccr12en
      labels(4) = op_tbar
      nlabel = 4
      if(.not.r12fix)then
        labels(5) = op_cba
        nlabel = 5
      endif
      call add_target(form_ccr12en0,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12en0,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12en0,op_ccr12en,tgt_info)
      call set_rule(form_ccr12en0,ttype_frm,INVARIANT,
     &              labels,nlabel,1,
     &              title_ccr12en0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12rs_t
      labels(2) = form_ccr12lg0
      labels(3) = op_omg
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_ccr12rs_t,ttype_frm,.true.,tgt_info)
      call set_dependency(form_ccr12rs_t,form_ccr12lg0,tgt_info)
      call set_dependency(form_ccr12rs_t,op_omg,tgt_info)
      call set_rule(form_ccr12rs_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_ccr12rs_t,1,tgt_info)

      if(.not.r12fix)then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12rs_c
        labels(2) = form_ccr12lg0
        labels(3) = op_omgr12
        labels(4) = op_cba
        labels(5) = ' '
        call add_target(form_ccr12rs_c,ttype_frm,.true.,tgt_info)
        call set_dependency(form_ccr12rs_c,form_ccr12lg0,tgt_info)
        call set_dependency(form_ccr12rs_c,op_omgr12,tgt_info)
        call set_rule(form_ccr12rs_c,ttype_frm,DERIVATIVE,
     &                labels,5,1,
     &                title_ccr12rs_c,1,tgt_info)
      endif

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

      ! CC ground state:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = fopt_ccr12_0
      labels(2) = form_ccr12en0
      labels(3) = form_ccr12rs_t
      ncat = 2
      nint = 0
      if(.not.r12fix)then
        labels(4) = form_ccr12rs_c
        ncat = 3
      endif
      call add_target(fopt_ccr12_0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12en0,tgt_info)
      call set_dependency(fopt_ccr12_0,form_ccr12rs_t,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_omgdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_topdef,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_b_def,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_v_def,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_x_def,tgt_info)
      if(.not.truncate.or.(truncate.and.(trunc_type.gt.0)))
     &           call set_dependency(fopt_ccr12_0,mel_p_def,tgt_info)
      if(.not.truncate.or.(truncate.and.(trunc_type.eq.1)))
     &           call set_dependency(fopt_ccr12_0,mel_z_def,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_rint,tgt_info)
      if (ansatz.gt.1)
     &     call set_dependency(fopt_ccr12_0,mel_c_def,tgt_info)
      if(.not.r12fix)then
        call set_dependency(fopt_ccr12_0,form_ccr12rs_c,tgt_info)
        call set_dependency(fopt_ccr12_0,mel_omgr12def,tgt_info)
      endif
      call set_dependency(fopt_ccr12_0,mel_c12def,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_ham,tgt_info)
      call set_dependency(fopt_ccr12_0,mel_ccr12en0def,tgt_info)      
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_ccr12_0,form_cchhat,tgt_info)
        call set_dependency(fopt_ccr12_0,mel_hhatdef,tgt_info)
        labels(1+ncat+1) = form_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_ccr12_0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*
      ! L0/E0:
      call add_target(mel_ccr12lg0,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_ccr12lg0,op_ccr12lg,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_ccr12lg0
      labels(2) = op_ccr12lg
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule(mel_ccr12lg0,ttype_opme,DEF_ME_LIST,
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

      ! CBAR list definition
      call add_target(mel_cbardef,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_cbardef,op_cba,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_cbar
      labels(2) = op_cba
      call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
      call set_rule(mel_cbardef,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

      ! C12  list definition
      call add_target(mel_c12def,ttype_opme,.false.,tgt_info)
      call set_dependency(mel_c12def,op_c12,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = mel_c12
      labels(2) = op_c12
      call me_list_parameters(-1,parameters,
     &       msc,0,1,0,0,.false.)
      call set_rule(mel_c12def,ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                parameters,1,tgt_info)

      if (.not.r12fix) then
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
      
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      if(.not.r12fix)then
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        if(.not.pz_eval)then
          call add_target(solve_ccr12_gs,ttype_gen,.true.,tgt_info)
        else
          call add_target(solve_ccr12_gs,ttype_gen,.false.,tgt_info)
        endif

        call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_b_inv,tgt_info)
c        call set_dependency(solve_ccr12_gs,mel_x_inv,tgt_info)
c        call set_dependency(solve_ccr12_gs,mel_b_dia,tgt_info)
        call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
        call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
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
        labels(7) = mel_ccr12en0
        labels(8) = fopt_ccr12_0
        if(trim(approx).eq.'A')then
          labels(9) = mel_b_inv ! or mel_b_inter
        else
          labels(9) = mel_b_inter
        endif
        labels(10) = mel_x_inter
        labels(11) = mel_ham
        call set_rule(solve_ccr12_gs,ttype_opme,SOLVENLEQ,
     &       labels,11,4,
     &       parameters,2,tgt_info)

      else
        ! totally symmetric dia for use below:
        call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

        call add_target(solve_ccr12_gs,ttype_gen,.true.,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_dia1,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_c12def,tgt_info)
        call set_dependency(solve_ccr12_gs,mel_cbardef,tgt_info)
        call set_dependency(solve_ccr12_gs,fopt_ccr12_0,tgt_info)
        call set_dependency(solve_ccr12_gs,eval_r12_inter,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_c12
        call set_rule(solve_ccr12_gs,ttype_opme,UNITY,
     &       labels,1,1,
     &       parameters,0,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = mel_cbar
        call set_rule(solve_ccr12_gs,ttype_opme,UNITY,
     &       labels,1,1,
     &       parameters,0,tgt_info)
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
     &       labels,6,2,
     &       parameters,2,tgt_info)
      endif

      return
      end
