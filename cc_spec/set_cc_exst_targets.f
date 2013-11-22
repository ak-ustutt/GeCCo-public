*----------------------------------------------------------------------*
      subroutine set_cc_exst_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in CC excited state calculations
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
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     me_label_l, me_label_r, me_label_lr,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(4)


      ! skip this section if not requested
      ncnt = is_keyword_set('calculate.excitation')
      if (ncnt.eq.0) return

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting targets for CC excited states ...'

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

      ! result of <L|R>
      call add_target(op_l_s_r,ttype_op,.false.,tgt_info)
      call set_rule(op_l_s_r,ttype_op,DEF_SCALAR,
     &              op_l_s_r,1,1,
     &              parameters,0,tgt_info)

      ! result of <L|A|R>
      call add_target(op_l_a_r,ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r,ttype_op,DEF_SCALAR,
     &              op_l_a_r,1,1,
     &              parameters,0,tgt_info)

      call add_target(op_l_a_r//'_1',ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r//'_1',ttype_op,DEF_SCALAR,
     &              op_l_a_r//'_1',1,1,
     &              parameters,0,tgt_info)
      call add_target(op_l_a_r//'_2',ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r//'_2',ttype_op,DEF_SCALAR,
     &              op_l_a_r//'_2',1,1,
     &              parameters,0,tgt_info)
      call add_target(op_l_a_r//'_3',ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r//'_3',ttype_op,DEF_SCALAR,
     &              op_l_a_r//'_3',1,1,
     &              parameters,0,tgt_info)
      call add_target(op_l_a_r//'_4',ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r//'_4',ttype_op,DEF_SCALAR,
     &              op_l_a_r//'_4',1,1,
     &              parameters,0,tgt_info)
      call add_target(op_l_a_r//'_5',ttype_op,.false.,tgt_info)
      call set_rule(op_l_a_r//'_5',ttype_op,DEF_SCALAR,
     &              op_l_a_r//'_5',1,1,
     &              parameters,0,tgt_info)

*----------------------------------------------------------------------*
*     Formulae:
*----------------------------------------------------------------------*

      ! right Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_a_r
      labels(2) = form_ccrs0
      labels(3) = op_a_r
      labels(4) = op_top
      labels(5) = op_r
      call add_target(form_cc_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc_a_r,form_ccrs0,tgt_info)
      call set_dependency(form_cc_a_r,op_a_r,tgt_info)
      call set_dependency(form_cc_a_r,op_r,tgt_info)
      call set_rule(form_cc_a_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_a_r,1,tgt_info)

      ! left Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_l_a
      labels(2) = form_cctbar_a
      labels(3) = op_l_a
      labels(4) = op_tbar
      labels(5) = op_l
      call add_target(form_cc_l_a,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc_l_a,form_cctbar_a,tgt_info)
      call set_dependency(form_cc_l_a,op_l_a,tgt_info)
      call set_dependency(form_cc_l_a,op_l,tgt_info)
      call set_rule(form_cc_l_a,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_l_a,1,tgt_info)

      ! make <L|R> formula
      call add_target(form_cc_l_s_r,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc_l_s_r,op_r,tgt_info)
      call set_dependency(form_cc_l_s_r,op_l,tgt_info)
      call set_dependency(form_cc_l_s_r,op_l_s_r,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_l_s_r
      labels(2) = op_l_s_r
      labels(3) = op_l_s_r
      labels(4) = op_l
      labels(5) = op_r
      labels(6) = op_l_s_r
      call expand_parameters(-1,parameters,3,
     &     title_cc_l_s_r,4,
     &     (/1,2,3,1/),-1,-1,
     &     0,0,0,0,0,0)
      call set_rule(form_cc_l_s_r,ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,6,1,
     &              parameters,3,tgt_info)

      ! make <L|A|R> formula
      call add_target(form_cc_l_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(form_cc_l_a_r,form_cclg0,tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r,tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r//'_1',tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r//'_2',tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r//'_3',tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r//'_4',tgt_info)
      call set_dependency(form_cc_l_a_r,op_l_a_r//'_5',tgt_info)
      call set_dependency(form_cc_l_a_r,op_r,tgt_info)
      call set_dependency(form_cc_l_a_r,op_l,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_SCRATCH'
      labels(2) = form_cclg0
      labels(3) = op_l_a_r
      labels(4) = op_tbar
      labels(5) = op_l
      call form_parameters(-1,
     &       parameters,2,
     &       'scratch',1,'-')
      call set_rule(form_cc_l_a_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_l_a_r
      labels(2) = 'FORM_SCRATCH'
      labels(3) = op_l_a_r
      labels(4) = op_top
      labels(5) = op_r
      call form_parameters(-1,
     &     parameters,2,
     &     title_cc_l_a_r,1,'-')
      call set_rule(form_cc_l_a_r,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_1'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_1'
      labels(4) = op_l
      labels(5) = op_r
      labels(6) = op_top
      call select_parameters(-1,
     &     parameters,3,
     &     2,0,1,(/1,1/),(/0/),(/2/))
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,6,1,
     &              parameters,3,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_2'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_2'
      labels(4) = op_l
      labels(5) = op_r
      labels(6) = op_top
      call select_parameters(-1,
     &     parameters,3,
     &     3,0,0,(/1,1,2/),0,0)
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,6,1,
     &              parameters,3,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_2'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_2'
      labels(4) = op_l
      labels(5) = op_r
      labels(6) = op_top
      call select_parameters(-1,
     &     parameters,3,
     &     3,0,0,(/1,1,2/),0,0)
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,6,1,
     &              parameters,3,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_3'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_3'
      labels(4) = op_l
      labels(5) = op_r
      call select_parameters(-1,
     &     parameters,3,
     &     2,0,0,(/1,2/),0,0)
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,5,1,
     &              parameters,3,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_4'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_4'
      labels(4) = op_l
      labels(5) = op_r
      call select_parameters(-1,
     &     parameters,3,
     &     2,0,0,(/2,1/),0,0)
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,5,1,
     &              parameters,3,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = trim(form_cc_l_a_r)//'_5'
      labels(2) = form_cc_l_a_r
      labels(3) = op_l_a_r//'_5'
      labels(4) = op_l
      labels(5) = op_r
      call select_parameters(-1,
     &     parameters,3,
     &     2,0,0,(/2,2/),0,0)
      call set_rule(form_cc_l_a_r,ttype_frm,SELECT_TERMS,
     &              labels,5,1,
     &              parameters,3,tgt_info)
      
*----------------------------------------------------------------------*
*     Optimized Formulae:
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC right-hand Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_a_r
      labels(2) = form_cc_a_r
      ncat = 1
      nint = 0
      call add_target(fopt_cc_a_r,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cc_a_r,form_cc_a_r,tgt_info)
      call set_dependency(fopt_cc_a_r,meldef_a_rex,tgt_info)
      call set_dependency(fopt_cc_a_r,meldef_rex,tgt_info)
      call set_dependency(fopt_cc_a_r,mel_topdef,tgt_info)
      call set_dependency(fopt_cc_a_r,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_cc_a_r,form_cchhat,tgt_info)
        call set_dependency(fopt_cc_a_r,mel_hhatdef,tgt_info)
        labels(3) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_cc_a_r,form_cchbar,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_hbar,tgt_info)
        labels(3) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_a_r,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! CC left-hand Jacobian transform
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_l_a
      labels(2) = form_cc_l_a
      ncat = 1
      nint = 0
      call add_target(fopt_cc_l_a,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cc_l_a,form_cc_l_a,tgt_info)
      call set_dependency(fopt_cc_l_a,meldef_lex_a,tgt_info)
      call set_dependency(fopt_cc_l_a,meldef_lex,tgt_info)
      call set_dependency(fopt_cc_l_a,mel_topdef,tgt_info)
      call set_dependency(fopt_cc_l_a,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_cc_l_a,form_cchhat,tgt_info)
        call set_dependency(fopt_cc_l_a,mel_hhatdef,tgt_info)
        labels(3) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_cc_a_r,form_cchbar,tgt_info)
        call set_dependency(fopt_cc_a_r,meldef_hbar,tgt_info)
        labels(3) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_l_a,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! <L|R>
      call add_target(fopt_cc_l_s_r,ttype_frm,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_l_s_r
      labels(2) = form_cc_l_s_r
      ncat = 1
      nint = 0
      call set_dependency(fopt_cc_l_s_r,form_cc_l_s_r,tgt_info)
      call set_dependency(fopt_cc_l_s_r,meldef_rex,tgt_info)
      call set_dependency(fopt_cc_l_s_r,meldef_lex,tgt_info)
      call set_dependency(fopt_cc_l_s_r,meldef_l_s_r,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_l_s_r,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      
      ! <L|A|R>
      call add_target(fopt_cc_l_a_r,ttype_frm,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_l_a_r
      labels(2) = form_cc_l_a_r
      labels(3) = form_cc_l_a_r//'_5'
      labels(4) = form_cc_l_a_r//'_4'
      labels(5) = form_cc_l_a_r//'_3'
      labels(6) = form_cc_l_a_r//'_2'
      labels(7) = form_cc_l_a_r//'_1'
      ncat = 1
      nint = 5
      call set_dependency(fopt_cc_l_a_r,form_cc_l_a_r,tgt_info)
      call set_dependency(fopt_cc_l_a_r,meldef_rex,tgt_info)
      call set_dependency(fopt_cc_l_a_r,meldef_lex,tgt_info)
      call set_dependency(fopt_cc_l_a_r,meldef_l_a_r,tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_l_a_r,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists:
*----------------------------------------------------------------------*
      call add_target(meldef_rex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_a_rex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_lex,ttype_opme,.false.,tgt_info)
      call add_target(meldef_lex_a,ttype_opme,.false.,tgt_info)
      call add_target(meldef_l_s_r,ttype_opme,.false.,tgt_info)
      call add_target(meldef_l_a_r,ttype_opme,.false.,tgt_info)
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
     &         msc,0,isym,0,0,.false.)
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
     &         msc,0,isym,0,0,.false.)
          call set_rule(meldef_a_rex,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
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
     &         msc,0,isym,0,0,.false.)
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
     &         msc,0,isym,0,0,.false.)
          call set_rule(meldef_lex_a,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        end do

        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle
          ! LE0.S.RE0
          call me_list_label(me_label,mel_l_s_r,isym,0,0,msc,.false.)
          call set_dependency(meldef_l_s_r,op_l,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_l_s_r
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_s_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          ! LE0.A.RE0
          call me_list_label(me_label,mel_l_a_r,isym,0,0,msc,.false.)
          call set_dependency(meldef_l_a_r,op_l,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_l_a_r
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          call me_list_label(me_label,'L1HR1',isym,0,0,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(me_label)
          labels(2) = op_l_a_r//'_1'
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          call me_list_label(me_label,'L1HcR1',isym,0,0,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(me_label)
          labels(2) = op_l_a_r//'_2'
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          call me_list_label(me_label,'L1HR2',isym,0,0,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(me_label)
          labels(2) = op_l_a_r//'_3'
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          call me_list_label(me_label,'L2HR1',isym,0,0,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(me_label)
          labels(2) = op_l_a_r//'_4'
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
          call me_list_label(me_label,'L2HR2',isym,0,0,msc,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = trim(me_label)
          labels(2) = op_l_a_r//'_5'
          call me_list_parameters(-1,parameters,
     &         abs(msc),0,1,0,0,.false.)
          call set_rule(meldef_l_a_r,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)

        end do

      end do

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      call add_target(solve_cc_rhex,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_cc_rhex,solve_cc_gs,tgt_info)
      call set_dependency(solve_cc_rhex,fopt_cc_a_r,tgt_info)
      call set_dependency(solve_cc_rhex,meldef_rex,tgt_info)
      call set_dependency(solve_cc_rhex,meldef_a_rex,tgt_info)

      call add_target(solve_cc_lhex,ttype_gen,.false.,tgt_info)
      call set_dependency(solve_cc_lhex,solve_cc_gs,tgt_info)
      call set_dependency(solve_cc_lhex,fopt_cc_l_a,tgt_info)
      call set_dependency(solve_cc_lhex,meldef_lex,tgt_info)
      call set_dependency(solve_cc_lhex,meldef_lex_a,tgt_info)

      needed = is_keyword_set('calculate.excitation.normalize').gt.0
      call add_target(norm_cc_lhex,ttype_gen,needed,tgt_info)
      call set_dependency(norm_cc_lhex,solve_cc_rhex,tgt_info)
      call set_dependency(norm_cc_lhex,solve_cc_lhex,tgt_info)
      call set_dependency(norm_cc_lhex,fopt_cc_l_s_r,tgt_info)      

      needed = is_keyword_set('calculate.excitation.analyze').gt.0
      call add_target(analyze_cc_ex,ttype_gen,needed,tgt_info)
      call set_dependency(analyze_cc_ex,solve_cc_rhex,tgt_info)
      call set_dependency(analyze_cc_ex,solve_cc_lhex,tgt_info)
      call set_dependency(analyze_cc_ex,norm_cc_lhex,tgt_info)      
      call set_dependency(analyze_cc_ex,fopt_cc_l_a_r,tgt_info)      

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_l_a_r
      call form_parameters(-1,
     &     parameters,2,
     &     title_cc_l_a_r,1,'-')
      call set_rule(analyze_cc_ex,ttype_frm,PRINT_FORMULA,
     &         labels,1,1,
     &         parameters,2,tgt_info)


      do icnt = 1, ncnt 
        call get_argument_value('calculate.excitation','sym',
     &       keycount=icnt,
     &       iarr=sym_arr)
        call get_argument_value('calculate.excitation','msc',
     &       keycount=icnt,
     &       ival=msc)
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          
          call me_list_label(me_label,mel_rex,isym,0,0,msc,.false.)
          call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_r
          call set_rule(solve_cc_rhex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call set_dependency(solve_cc_rhex,dia_label,tgt_info)
          call solve_parameters(-1,parameters,2,1,sym_arr(isym),'DIA')
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = dia_label
          labels(3) = op_a_r
          labels(4) = op_r
          labels(5) = fopt_cc_a_r
          call set_rule(solve_cc_rhex,ttype_opme,SOLVEEVP,
     &         labels,5,1,
     &         parameters,2,tgt_info)
        end do
      
        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          
          call me_list_label(me_label,mel_lex,isym,0,0,msc,.false.)
          call me_list_label(dia_label,mel_dia,isym,0,0,0,.false.)

          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = op_l
          call set_rule(solve_cc_lhex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call set_dependency(solve_cc_lhex,dia_label,tgt_info)
          call solve_parameters(-1,parameters,2,1,sym_arr(isym),'DIA')
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label
          labels(2) = dia_label
          labels(3) = op_l_a
          labels(4) = op_l
          labels(5) = fopt_cc_l_a
          call set_rule(solve_cc_lhex,ttype_opme,SOLVEEVP,
     &         labels,5,1,
     &         parameters,2,tgt_info)
        end do

        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          

          call me_list_label(me_label_l,mel_lex,isym,0,0,msc,.false.)
          call me_list_label(me_label_r,mel_rex,isym,0,0,msc,.false.)
          call me_list_label(me_label_lr,mel_l_s_r,isym,0,0,msc,.false.)

          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_l
          labels(2) = op_l
          call set_rule(norm_cc_lhex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_r
          labels(2) = op_r
          call set_rule(norm_cc_lhex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_lr
          labels(2) = op_l_s_r
          call set_rule(norm_cc_lhex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = fopt_cc_l_s_r
          call set_rule(norm_cc_lhex,ttype_opme,EVAL,
     &         labels,1,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_l
          labels(2) = me_label_l
          labels(3) = me_label_lr
          call scale_parameters(-1,parameters,3,0,0,0d0,1)
          call set_rule(norm_cc_lhex,ttype_opme,SCALE,
     &         labels,3,1,
     &         parameters,1,tgt_info)
        end do

        do isym = 1, orb_info%nsym
          if (sym_arr(isym).eq.0) cycle          

          call me_list_label(me_label_l,mel_lex,isym,0,0,msc,.false.)
          call me_list_label(me_label_r,mel_rex,isym,0,0,msc,.false.)
          call me_list_label(me_label_lr,mel_l_s_r,isym,0,0,msc,.false.)

          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_l
          labels(2) = op_l
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_r
          labels(2) = op_r
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          labels(1) = me_label_lr
          labels(2) = op_l_a_r
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          labels(1:10)(1:len_target_name) = ' '
          call me_list_label(me_label_lr,'L1HR1',isym,0,0,msc,.false.)
          labels(1) = me_label_lr
          labels(2) = op_l_a_r//'_1'
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call me_list_label(me_label_lr,'L1HcR1',isym,0,0,msc,.false.)
          labels(1) = me_label_lr
          labels(2) = op_l_a_r//'_2'
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call me_list_label(me_label_lr,'L1HR2',isym,0,0,msc,.false.)
          labels(1) = me_label_lr
          labels(2) = op_l_a_r//'_3'
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call me_list_label(me_label_lr,'L2HR1',isym,0,0,msc,.false.)
          labels(1) = me_label_lr
          labels(2) = op_l_a_r//'_4'
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)
          call me_list_label(me_label_lr,'L2HR2',isym,0,0,msc,.false.)
          labels(1) = me_label_lr
          labels(2) = op_l_a_r//'_5'
          call set_rule(analyze_cc_ex,ttype_opme,ASSIGN_ME2OP,
     &         labels,2,1,
     &         parameters,0,tgt_info)

          labels(1:10)(1:len_target_name) = ' '
          labels(1) = fopt_cc_l_a_r
          call set_rule(analyze_cc_ex,ttype_opme,EVAL,
     &         labels,1,1,
     &         parameters,0,tgt_info)
        end do


      end do
      

      return
      end
