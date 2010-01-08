*----------------------------------------------------------------------*
      subroutine set_cc_ipst_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed in CC ionized state calculations
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

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
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed, setr12
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      ! skip this section if no IP calculation requested
      icnt = is_keyword_set('calculate.ionization')
      if (icnt.eq.0) return

      setr12 = is_keyword_set('method.R12').gt.0

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting targets for CC ionized states ...'

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! particle removing response
      ! right-response vector
      call add_target(op_rip,ttype_op,.false.,tgt_info)
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      call xop_parameters(-1,parameters,
     &                   .false.,min_rank,max_rank,-1,1)
      call set_rule(op_rip,ttype_op,DEF_EXCITATION,
     &              op_rip,1,1,
     &              parameters,1,tgt_info)

      ! the appropriate diagonal:
      call add_target(op_dia_ip,ttype_op,.false.,tgt_info)
      call set_dependency(op_dia_ip,op_rip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rip,.false.)
      call set_rule(op_dia_ip,ttype_op,CLONE_OP,
     &              op_dia_ip,1,1,
     &              parameters,1,tgt_info)
      
      ! right-response vector times Jacobian
      call add_target(op_a_rip,ttype_op,.false.,tgt_info)
      call set_dependency(op_a_rip,op_rip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rip,.false.)
      call set_rule(op_a_rip,ttype_op,CLONE_OP,
     &              op_a_rip,1,1,
     &              parameters,1,tgt_info)

      ! left vector
      call add_target(op_lip,ttype_op,.false.,tgt_info)
      call set_dependency(op_lip,op_rip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_rip,.true.)
      call set_rule(op_lip,ttype_op,CLONE_OP,
     &              op_lip,1,1,
     &              parameters,1,tgt_info)

      ! left vector times Jacobian
      call add_target(op_lip_a,ttype_op,.false.,tgt_info)
      call set_dependency(op_lip_a,op_lip,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_lip,.false.)
      call set_rule(op_lip_a,ttype_op,CLONE_OP,
     &              op_lip_a,1,1,
     &              parameters,1,tgt_info)
      

*----------------------------------------------------------------------*
*     Formulae:
*----------------------------------------------------------------------*

      ! right Jacobian transform with IP operator
      call add_target(form_cc_a_rip,ttype_frm,.false.,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_cc_a_rip
      if (.not.setr12) then
        labels(2) = form_ccrs0
        call set_dependency(form_cc_a_rip,form_ccrs0,tgt_info)
      else
        labels(2) = form_ccr12rs_t
        call set_dependency(form_cc_a_rip,form_ccr12rs_t,tgt_info)
      end if
      labels(3) = op_a_rip
      labels(4) = op_top
      labels(5) = op_rip
      call set_dependency(form_cc_a_rip,op_a_rip,tgt_info)
      call set_dependency(form_cc_a_rip,op_rip,tgt_info)
      call set_rule(form_cc_a_rip,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_cc_a_rip,1,tgt_info)

c      ! left Jacobian transform with IP operator
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = form_cc_lip_a
c      labels(2) = form_cctbar_a
c      labels(3) = op_lip_a
c      labels(4) = op_tbar
c      labels(5) = op_lip
c      call add_target(form_cc_lip_a,ttype_frm,.false.,tgt_info)
c      call set_dependency(form_cc_lip_a,form_ccrs0,tgt_info)
c      call set_dependency(form_cc_lip_a,op_lip_a,tgt_info)
c      call set_dependency(form_cc_lip_a,op_lip,tgt_info)
c      call set_rule(form_cc_lip_a,ttype_frm,DERIVATIVE,
c     &              labels,5,1,
c     &              title_cc_lip_a,1,tgt_info)

*----------------------------------------------------------------------*
*     Optimized Formulae:
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! CC right-hand Jacobian transform with IP operator
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_cc_a_rip
      labels(2) = form_cc_a_rip
      ncat = 1
      nint = 0
      call add_target(fopt_cc_a_rip,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_cc_a_rip,form_cc_a_rip,tgt_info)
      call set_dependency(fopt_cc_a_rip,meldef_a_rip,tgt_info)
      call set_dependency(fopt_cc_a_rip,meldef_rip,tgt_info)
      call set_dependency(fopt_cc_a_rip,mel_topdef,tgt_info)
      call set_dependency(fopt_cc_a_rip,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_cc_a_rip,form_cchhat,tgt_info)
        call set_dependency(fopt_cc_a_rip,mel_hhatdef,tgt_info)
        labels(3) = form_cchhat
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_cc_a_rip,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists:
*----------------------------------------------------------------------*
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
     &         msc,0,isym,0,ms,.false.)
        call set_rule(meldef_rip,ttype_opme,DEF_ME_LIST,
     &         labels,2,1,
     &         parameters,1,tgt_info)
        ! A.RI0
        call me_list_label(me_label,mel_a_rip,isym,0,ms,msc,.false.)
        labels(1:10)(1:len_target_name) = ' '
        labels(1) = me_label
        labels(2) = op_a_rip
        call me_list_parameters(-1,parameters,
     &       msc,0,isym,0,ms,.false.)
        call set_rule(meldef_a_rip,ttype_opme,DEF_ME_LIST,
     &       labels,2,1,
     &       parameters,1,tgt_info)
      end do

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      ! solve the eigenvalue equations
      call add_target2(solve_cc_rhip,.true.,tgt_info)
      if (.not.setr12)
     &     call set_dependency(solve_cc_rhip,solve_cc_gs,tgt_info)
      if (setr12)
     &     call set_dependency(solve_cc_rhip,solve_ccr12_gs,tgt_info)
      call set_dependency(solve_cc_rhip,fopt_cc_a_rip,tgt_info)
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
        call set_rule2(solve_cc_rhip,DEF_ME_LIST,tgt_info)
        call set_arg(solve_cc_rhip,DEF_ME_LIST,'LIST',1,tgt_info,
     &       val_label=(/dia_label/))
        call set_arg(solve_cc_rhip,DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &       val_label=(/op_dia_ip/))
        call set_arg(solve_cc_rhip,DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &       val_int=(/0/))
        call set_arg(solve_cc_rhip,DEF_ME_LIST,'IRREP',1,tgt_info,
     &       val_int=(/isym/))
        call set_arg(solve_cc_rhip,DEF_ME_LIST,'MS',1,tgt_info,
     &       val_int=(/ms/))

        call set_rule2(solve_cc_rhip,PRECONDITIONER,tgt_info)
        call set_arg(solve_cc_rhip,PRECONDITIONER,'LIST_PRC',1,tgt_info,
     &       val_label=(/dia_label/))
        call set_arg(solve_cc_rhip,PRECONDITIONER,'LIST_INP',1,tgt_info,
     &       val_label=(/mel_ham/))
        call set_arg(solve_cc_rhip,PRECONDITIONER,'MODE',1,tgt_info,
     &       val_str='dia-F')

        ! b) solve the eigenvalue equation
        call set_rule2(solve_cc_rhip,SOLVEEVP,tgt_info)
        call set_arg(solve_cc_rhip,SOLVEEVP,'LIST_OPT',1,tgt_info,
     &       val_label=(/me_label/))
        call set_arg(solve_cc_rhip,SOLVEEVP,'MODE',1,tgt_info,
     &       val_str='DIA')
        call set_arg(solve_cc_rhip,SOLVEEVP,'N_ROOTS',1,tgt_info,
     &       val_int=(/sym_arr(isym)/))
        call set_arg(solve_cc_rhip,SOLVEEVP,'LIST_PRC',1,tgt_info,
     &       val_label=(/dia_label/))
        call set_arg(solve_cc_rhip,SOLVEEVP,'OP_MVP',1,tgt_info,
     &       val_label=(/op_a_rip/))
        call set_arg(solve_cc_rhip,SOLVEEVP,'OP_SVP',1,tgt_info,
     &       val_label=(/op_rip/))
        call set_arg(solve_cc_rhip,SOLVEEVP,'FORM',1,tgt_info,
     &       val_label=(/fopt_cc_a_rip/))

        ! c) remove the diagonal list
        call set_rule2(solve_cc_rhip,DELETE_ME_LIST,tgt_info)
        call set_arg(solve_cc_rhip,DELETE_ME_LIST,'LIST',1,tgt_info,
     &       val_label=(/dia_label/))
      end do


      return
      end
