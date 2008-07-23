*----------------------------------------------------------------------*
      subroutine set_ccr12_exst_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in CC-R12 excited state calculations
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
     &     isim, ncat, nint, icnt, ncnt,
     &     isym, ms, msc, sym_arr(8), r12op
      logical ::
     &     needed, r12fix
      character(len_target_name) ::
     &     me_label, mep_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)


      ! skip this section if not requested
      ncnt = is_keyword_set('calculate.excitation')
      if (ncnt.eq.0) return

      if (iprlvl.gt.0)
     &    write(luout,*) 'setting targets for CC-R12 excited states ...'

      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','r12op',ival=r12op)

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! particle conserving response
      ! right-response vector
      call add_target(op_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_r,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_r,ttype_op,CLONE_OP,
     &              op_r,1,1,
     &              parameters,1,tgt_info)
      
      ! right-response vector times Jacobian
      call add_target(op_a_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_a_r,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_a_r,ttype_op,CLONE_OP,
     &              op_a_r,1,1,
     &              parameters,1,tgt_info)

      ! T vector times Metric
      call add_target(op_s_t,ttype_op,.false.,tgt_info)
      call set_dependency(op_s_t,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_s_t,ttype_op,CLONE_OP,
     &              op_s_t,1,1,
     &              parameters,1,tgt_info)
      
      ! T' vector times Metric
      call add_target(op_s_c,ttype_op,.false.,tgt_info)
      call set_dependency(op_s_c,op_cex,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
      call set_rule(op_s_c,ttype_op,CLONE_OP,
     &              op_s_c,1,1,
     &              parameters,1,tgt_info)

      ! right-response vector times Metric
      call add_target(op_s_r,ttype_op,.false.,tgt_info)
      call set_dependency(op_s_r,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_s_r,ttype_op,CLONE_OP,
     &              op_s_r,1,1,
     &              parameters,1,tgt_info)

      ! if applicable: define here additional R12-part of response:
      ! R' response vector
      call add_target(op_rp,ttype_op,.false.,tgt_info)
      call set_dependency(op_rp,op_cex,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
      call set_rule(op_rp,ttype_op,CLONE_OP,
     &              op_rp,1,1,
     &              parameters,1,tgt_info)
      ! R' times Jacobian
      call add_target(op_a_rp,ttype_op,.false.,tgt_info)
      call set_dependency(op_a_rp,op_cex,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
      call set_rule(op_a_rp,ttype_op,CLONE_OP,
     &              op_a_rp,1,1,
     &              parameters,1,tgt_info)
      ! R' times Metric
      call add_target(op_s_rp,ttype_op,.false.,tgt_info)
      call set_dependency(op_s_rp,op_cex,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_cex,.false.)
      call set_rule(op_s_rp,ttype_op,CLONE_OP,
     &              op_s_rp,1,1,
     &              parameters,1,tgt_info)


      ! left-response vector
      call add_target(op_l,ttype_op,.false.,tgt_info)
      call set_dependency(op_l,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_l,ttype_op,CLONE_OP,
     &              op_l,1,1,
     &              parameters,1,tgt_info)

      ! left-response vector times Jacobian
      call add_target(op_l_a,ttype_op,.false.,tgt_info)
      call set_dependency(op_l_a,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_l_a,ttype_op,CLONE_OP,
     &              op_l_a,1,1,
     &              parameters,1,tgt_info)



*----------------------------------------------------------------------*
*     Formulae:
*----------------------------------------------------------------------*

      ! right Jacobian transform
      ! 2-components: make derivative wrt both comp.s
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_a_r
      labels(2) = form_ccr12rs_t
      labels(3) = op_a_r
      call add_target(form_cc_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc_a_r,form_ccr12rs_t,tgt_info)
      call set_dependency(form_cc_a_r,op_a_r,tgt_info)
      call set_dependency(form_cc_a_r,op_r,tgt_info)
      if (r12fix) then
        labels(4) = op_top
        labels(5) = op_r
        call form_parameters(-1,
     &       parameters,2,
     &       title_cc_a_r,1,'-')
        call set_rule(form_cc_a_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      else
        labels(4) = op_top
        labels(5) = op_cex
        labels(6) = op_r
        labels(7) = op_rp
        call set_dependency(form_cc_a_r,op_rp,tgt_info)
        call form_parameters(-1,
     &       parameters,2,
     &       title_cc_a_r,2,'-')
        call set_rule(form_cc_a_r,ttype_frm,DERIVATIVE,
     &              labels,7,1,
     &              parameters,2,tgt_info)
      end if

      ! 2-components: add here derivative of rs_c component
      if (.not.r12fix) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_cc_a_rp
        labels(2) = form_ccr12rs_c
        labels(3) = op_a_rp
        call add_target(form_cc_a_rp,ttype_frm,.false.,tgt_info)
        call set_dependency(form_cc_a_rp,form_ccr12rs_c,tgt_info)
        call set_dependency(form_cc_a_rp,op_a_rp,tgt_info)
        call set_dependency(form_cc_a_rp,op_rp,tgt_info)
        labels(4) = op_top
        labels(5) = op_cex
        labels(6) = op_r
        labels(7) = op_rp
        call form_parameters(-1,
     &       parameters,2,
     &       title_cc_a_rp,2,'-')
        call set_rule(form_cc_a_rp,ttype_frm,DERIVATIVE,
     &              labels,7,1,
     &              parameters,2,tgt_info)
      end if

      call add_target(form_ccr12_s_r,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ccr12_s_r,form_ccr12_s0,tgt_info)
      call set_dependency(form_ccr12_s_r,op_s_t,tgt_info)
      call set_dependency(form_ccr12_s_r,op_s_r,tgt_info)
      call set_dependency(form_ccr12_s_r,op_r,tgt_info)
      ! we also need the transform with the metric
      ! 2-components: make derivative wrt both comp.s
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ccr12_s_t
      labels(2) = form_ccr12_s0
      labels(3) = op_s_t
      labels(4) = op_tbar
      labels(5) = ' '
      call form_parameters(-1,
     &       parameters,2,
     &       title_ccr12_s_t,1,'-')
      call set_rule(form_ccr12_s_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      if (r12fix) then
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s_r
        labels(2) = form_ccr12_s_t
        labels(3) = op_s_r
        labels(4) = op_top
        labels(5) = op_r
        call form_parameters(-1,
     &       parameters,2,
     &       title_ccr12_s_r,1,'-')
        call set_rule(form_ccr12_s_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      else
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s_r
        labels(2) = form_ccr12_s_t
        labels(3) = op_s_r
        labels(4) = op_top
        labels(5) = op_cex
        labels(6) = op_r
        labels(7) = op_rp
        call form_parameters(-1,
     &       parameters,2,
     &       title_ccr12_s_t,2,'-')
        call set_rule(form_ccr12_s_r,ttype_frm,DERIVATIVE,
     &              labels,7,1,
     &              parameters,2,tgt_info)
      end if

      ! 2-components: add here derivative of rs_c component

      if (.not.r12fix) then
        call add_target(form_ccr12_s_rp,ttype_frm,.false.,tgt_info)
        call set_dependency(form_ccr12_s_rp,form_ccr12_s0,tgt_info)
        call set_dependency(form_ccr12_s_rp,op_s_c,tgt_info)
        call set_dependency(form_ccr12_s_rp,op_s_rp,tgt_info)
        call set_dependency(form_ccr12_s_rp,op_r,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s_c
        labels(2) = form_ccr12_s0
        labels(3) = op_s_c
        labels(4) = op_cexbar
        labels(5) = ' '
        call form_parameters(-1,
     &       parameters,2,
     &       title_ccr12_s_rp,1,'-')
        call set_rule(form_ccr12_s_rp,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = form_ccr12_s_rp
        labels(2) = form_ccr12_s_c
        labels(3) = op_s_rp
        labels(4) = op_top
        labels(5) = op_cex
        labels(6) = op_r
        labels(7) = op_rp
        call form_parameters(-1,
     &       parameters,2,
     &       title_ccr12_s_rp,2,'-')
        call set_rule(form_ccr12_s_rp,ttype_frm,DERIVATIVE,
     &              labels,7,1,
     &              parameters,2,tgt_info)
      end if

c      ! left Jacobian transform
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_cc_l_a
c      labels(2) = form_cctbar_a
c      labels(3) = op_l_a
c      labels(4) = op_tbar
c      labels(5) = op_l
c      call add_target(form_cc_l_a,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_cc_l_a,form_cctbar_a,tgt_info)
c      call set_dependency(form_cc_l_a,op_l_a,tgt_info)
c      call set_dependency(form_cc_l_a,op_l,tgt_info)
c      call set_rule(form_cc_l_a,ttype_frm,DERIVATIVE,
c     &              labels,5,1,
c     &              title_cc_l_a,1,tgt_info)


*----------------------------------------------------------------------*
*     Optimized Formulae:
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC right-hand Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_a_r
      labels(2) = form_cc_a_r
      labels(3) = form_ccr12_s_r
      ncat = 2
      nint = 0
      call add_target(fopt_cc_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cc_a_r,form_cc_a_r,tgt_info)
      call set_dependency(fopt_cc_a_r,form_ccr12_s_r,tgt_info)
      call set_dependency(fopt_cc_a_r,meldef_s_rex,tgt_info)
      call set_dependency(fopt_cc_a_r,meldef_a_rex,tgt_info)
      call set_dependency(fopt_cc_a_r,meldef_rex,tgt_info)
      call set_dependency(fopt_cc_a_r,mel_topdef,tgt_info)
      call set_dependency(fopt_cc_a_r,mel_ham,tgt_info)
      if (.not.r12fix) then
        labels(4) = form_cc_a_rp
        labels(5) = form_ccr12_s_rp
        call set_dependency(fopt_cc_a_r,form_cc_a_rp,tgt_info)
        call set_dependency(fopt_cc_a_r,form_ccr12_s_rp,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_s_rpex,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_a_rpex,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_rpex,tgt_info)
        ncat = 4
      end if
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_cc_a_r,form_cchhat,tgt_info)
        call set_dependency(fopt_cc_a_r,mel_hhatdef,tgt_info)
        labels(1+ncat+1) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_cc_a_r,form_cchbar,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_hbar,tgt_info)
        labels(1+ncat+1) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_a_r,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

c      ! CC left-hand Jacobian transform
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = fopt_cc_l_a
c      labels(2) = form_cc_l_a
c      ncat = 1
c      nint = 0
c      call add_target(fopt_cc_l_a,ttype_frm,.false.,tgt_info)
c      call set_dependency(fopt_cc_l_a,form_cc_l_a,tgt_info)
c      call set_dependency(fopt_cc_l_a,meldef_lex_a,tgt_info)
c      call set_dependency(fopt_cc_l_a,meldef_lex,tgt_info)
c      call set_dependency(fopt_cc_l_a,mel_topdef,tgt_info)
c      call set_dependency(fopt_cc_l_a,mel_ham,tgt_info)
c      if (isim.eq.1) then
c        nint = 1
c        call set_dependency(fopt_cc_l_a,form_cchhat,tgt_info)
c        call set_dependency(fopt_cc_l_a,mel_hhatdef,tgt_info)
c        labels(3) = form_cchhat
c      else if (isim.eq.2) then
c        nint = 1
c        call set_dependency(fopt_cc_a_r,form_cchbar,tgt_info)
c        call set_dependency(fopt_cc_a_r,meldef_hbar,tgt_info)
c        labels(3) = form_cchbar
c      end if
c      call opt_parameters(-1,parameters,ncat,nint)
c      call set_rule(fopt_cc_l_a,ttype_frm,OPTIMIZE,
c     &              labels,ncat+nint+1,1,
c     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists:
*----------------------------------------------------------------------*
      call add_target(meldef_rex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_a_rex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_s_rex,ttype_opme,.false.,tgt_info)
      if (.not.r12fix) then
        call add_target(meldef_rpex,ttype_opme,.false.,tgt_info)
        call add_target(meldef_a_rpex,ttype_opme,.false.,tgt_info)
        call add_target(meldef_s_rpex,ttype_opme,.false.,tgt_info)
      end if
      call add_target(meldef_lex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_lex_a,ttype_opme,.false.,tgt_info)
      do icnt = 1, ncnt 
        call get_argument_value('calculate.excitation','sym',
     &       keycount=icnt,
     &       iarr=sym_arr)
        call get_argument_value('calculate.excitation','msc',
     &       keycount=icnt,
     &       ival=msc)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ! RE0
          call me_list_label(me_label,mel_rex,isym,0,0,msc,.false.)
          call set_dependency(meldef_rex,op_r,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_r
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,0)
          call set_rule(meldef_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! A.RE0
          call me_list_label(me_label,mel_a_rex,isym,0,0,msc,.false.)
          call set_dependency(meldef_a_rex,op_a_r,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_a_r
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,0)
          call set_rule(meldef_a_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! S.RE0
          call me_list_label(me_label,mel_s_rex,isym,0,0,msc,.false.)
          call set_dependency(meldef_s_rex,op_s_r,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_s_r
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,0)
          call set_rule(meldef_s_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          if (.not.r12fix) then
            ! RE0''
            call me_list_label(me_label,mel_rpex,isym,0,0,msc,.false.)
            call set_dependency(meldef_rpex,op_rp,tgt_info)
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = op_rp
            call me_list_parameters(-1,parameters,
     &           msc,0,isym,0,0)
            call set_rule(meldef_rpex,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)
            ! A.RE0
            call me_list_label(me_label,mel_a_rpex,isym,0,0,msc,.false.)
            call set_dependency(meldef_a_rpex,op_a_rp,tgt_info)
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = op_a_rp
            call me_list_parameters(-1,parameters,
     &           msc,0,isym,0,0)
            call set_rule(meldef_a_rpex,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)
            ! S.RE0
            call me_list_label(me_label,mel_s_rpex,isym,0,0,msc,.false.)
            call set_dependency(meldef_s_rpex,op_s_rp,tgt_info)
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = op_s_rp
            call me_list_parameters(-1,parameters,
     &           msc,0,isym,0,0)
            call set_rule(meldef_s_rpex,ttype_opme,DEF_ME_LIST,
     &           labels,2,1,
     &           parameters,1,tgt_info)            
          end if ! .not.r12fix

        end do

        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ! LE0
          call me_list_label(me_label,mel_lex,isym,0,0,msc,.false.)
          call set_dependency(meldef_lex,op_l,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_l
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,0)
          call set_rule(meldef_lex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! LE0.A
          call me_list_label(me_label,mel_lex_a,isym,0,0,msc,.false.)
          call set_dependency(meldef_lex_a,op_l_a,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_l_a
          call me_list_parameters(-1,parameters,
     &         msc,0,isym,0,0)
          call set_rule(meldef_lex_a,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        end do
      end do

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      call add_target(solve_cc_rhex,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_cc_rhex,solve_ccr12_gs,tgt_info)
      call set_dependency(solve_cc_rhex,fopt_cc_a_r,tgt_info)
      call set_dependency(solve_cc_rhex,meldef_rex,tgt_info)
      call set_dependency(solve_cc_rhex,meldef_a_rex,tgt_info)
      call set_dependency(solve_cc_rhex,meldef_s_rex,tgt_info)

c      call add_target(solve_cc_lhex,ttype_gen,.false.,tgt_info)
c      call set_dependency(solve_cc_lhex,solve_ccr12_gs,tgt_info)
c      call set_dependency(solve_cc_lhex,fopt_cc_l_a,tgt_info)
c      call set_dependency(solve_cc_lhex,meldef_lex,tgt_info)
c      call set_dependency(solve_cc_lhex,meldef_lex_a,tgt_info)

      do icnt = 1, ncnt 
        call get_argument_value('calculate.excitation','sym',
     &       keycount=icnt,
     &       iarr=sym_arr)
        call get_argument_value('calculate.excitation','msc',
     &       keycount=icnt,
     &       ival=msc)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          
          if (r12fix) then
            call me_list_label(me_label,mel_rex,isym,0,0,msc,.false.)
            call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)
            call set_dependency(solve_cc_rhex,dia_label,tgt_info)
            call solve_parameters(-1,parameters,2,1,sym_arr(isym),'DIA')
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = dia_label
            labels(3) = op_a_r
            labels(4) = op_s_r
            labels(5) = fopt_cc_a_r
            call set_rule(solve_cc_rhex,ttype_opme,SOLVEEVP,
     &           labels,5,1,
     &           parameters,2,tgt_info)
          else
            call me_list_label(me_label,mel_rex,isym,0,0,msc,.false.)
            call me_list_label(mep_label,mel_rpex,isym,0,0,msc,.false.)
            call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)
            call set_dependency(solve_cc_rhex,dia_label,tgt_info)
c            call set_dependency(solve_cc_rhex,me_bprc,tgt_info)
c            call set_dependency(solve_cc_rhex,me_xprc,tgt_info)
            call set_dependency(solve_cc_rhex,'BPRCLIST',tgt_info)
            call set_dependency(solve_cc_rhex,'XPRCLIST',tgt_info)
            call solve_parameters(-1,parameters,2,2,sym_arr(isym),
     &           'DIA/BLK')
            labels(1:10)(1:len_target_name) = ' '
            labels(1) = me_label
            labels(2) = mep_label
            labels(3) = dia_label
            labels(4) = dia_label
            labels(5) = op_a_r
            labels(6) = op_a_rp
            labels(7) = op_s_r
            labels(8) = op_s_rp
            labels(9) = fopt_cc_a_r
c            labels(10) = me_bprc
c            labels(11)= me_xprc
            labels(10) = 'BPRCLIST'
            labels(11)= 'XPRCLIST'
            labels(12) = mel_ham
            call set_rule(solve_cc_rhex,ttype_opme,SOLVEEVP,
     &           labels,12,2,
     &           parameters,2,tgt_info)
          end if
        end do
      
c        do isym = 1, orb_info%nsym
c          if (sym_arr(isym).eq.0) cycle          
c          call me_list_label(me_label,mel_lex,isym,0,0,msc,.false.)
c          call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)
c          call set_dependency(solve_cc_lhex,dia_label,tgt_info)
c          call solve_parameters(-1,parameters,2,1,sym_arr(isym),'DIA')
c          labels(1:10)(1:len_target_name) = ' '
c          labels(1) = me_label
c          labels(2) = dia_label
c          labels(3) = op_l_a
c          labels(4) = fopt_cc_l_a
c          call set_rule(solve_cc_lhex,ttype_opme,SOLVEEVP,
c     &         labels,4,1,
c     &         parameters,2,tgt_info)
c        end do
      end do
      
      return
      end
