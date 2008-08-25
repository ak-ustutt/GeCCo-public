*----------------------------------------------------------------------*
      subroutine set_ecc_special_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets needed specifically in ECC calculations
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
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting special targets for ECC ...'

      msc = 1 ! presently assuming closed shell

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Lagrange functional 
      call add_target(op_ecclg,ttype_op,.false.,tgt_info)
      call set_rule(op_ecclg,ttype_op,DEF_SCALAR,
     &              op_ecclg,1,1,
     &              parameters,0,tgt_info)

      ! T Residual
      call add_target(op_omg_t,ttype_op,.false.,tgt_info)
      call set_dependency(op_omg_t,op_top,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_top,.false.)
      call set_rule(op_omg_t,ttype_op,CLONE_OP,
     &              op_omg_t,1,1,
     &              parameters,1,tgt_info)
      
      ! TBAR Residual
      call add_target(op_omg_l,ttype_op,.false.,tgt_info)
      call set_dependency(op_omg_l,op_tbar,tgt_info)
      call cloneop_parameters(-1,parameters,
     &                        op_tbar,.false.)
      call set_rule(op_omg_l,ttype_op,CLONE_OP,
     &              op_omg_l,1,1,
     &              parameters,1,tgt_info)
      
*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_ecclg0
      labels(2) = op_ecclg
      labels(3) = op_ham
      labels(4) = op_tbar
      labels(5) = op_top
      call add_target(form_ecclg0,ttype_frm,.false.,tgt_info)
      call set_dependency(form_ecclg0,op_ecclg,tgt_info)
      call set_dependency(form_ecclg0,op_ham,tgt_info)
      call set_dependency(form_ecclg0,op_tbar,tgt_info)
      call set_dependency(form_ecclg0,op_top,tgt_info)
      call set_rule(form_ecclg0,ttype_frm,DEF_ECC_LAGRANGIAN,
     &              labels,5,1,
     &              title_ecclg0,1,tgt_info)

      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_eccrs0_t
      labels(2) = form_ecclg0
      labels(3) = op_omg_t
      labels(4) = op_tbar
      labels(5) = ' '
      call add_target(form_eccrs0_t,ttype_frm,.false.,tgt_info)
      call set_dependency(form_eccrs0_t,form_ecclg0,tgt_info)
      call set_dependency(form_eccrs0_t,op_omg_t,tgt_info)
      call set_rule(form_eccrs0_t,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_eccrs0_t,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = form_eccrs0_l
      labels(2) = form_ecclg0
      labels(3) = op_omg_l
      labels(4) = op_top
      labels(5) = ' '
      call add_target(form_eccrs0_l,ttype_frm,.false.,tgt_info)
      call set_dependency(form_eccrs0_l,form_ecclg0,tgt_info)
      call set_dependency(form_eccrs0_l,op_omg_l,tgt_info)
      call set_rule(form_eccrs0_l,ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              title_eccrs0_l,1,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*
      call get_argument_value('calculate.routes','simtraf',ival=isim)

      ! ECC ground state:
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = fopt_eccrs0
      labels(2) = form_ecclg0
      labels(3) = form_eccrs0_l
      labels(4) = form_eccrs0_t
      ncat = 3
      nint = 0
      call add_target(fopt_eccrs0,ttype_frm,.false.,tgt_info)
      call set_dependency(fopt_eccrs0,form_eccrs0_l,tgt_info)
      call set_dependency(fopt_eccrs0,form_eccrs0_t,tgt_info)
      call set_dependency(fopt_eccrs0,mel_omgdef,tgt_info)
      call set_dependency(fopt_eccrs0,mel_topdef,tgt_info)
      call set_dependency(fopt_eccrs0,mel_tbardef,tgt_info)
      call set_dependency(fopt_eccrs0,meldef_ecclg0,tgt_info)
      call set_dependency(fopt_eccrs0,meldef_omg_t,tgt_info)
      call set_dependency(fopt_eccrs0,meldef_omg_l,tgt_info)
      call set_dependency(fopt_eccrs0,mel_ham,tgt_info)
      if (isim.eq.1) then
        nint = 1
        call set_dependency(fopt_eccrs0,form_cchhat,tgt_info)
        call set_dependency(fopt_eccrs0,mel_hhatdef,tgt_info)
        labels(1+ncat+1) = form_cchhat
      else if (isim.eq.2) then
        nint = 1
        call set_dependency(fopt_eccrs0,form_cchbar,tgt_info)
        call set_dependency(fopt_eccrs0,meldef_hbar,tgt_info)
        labels(1+ncat+1) = form_cchbar
      end if
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule(fopt_eccrs0,ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! L0/E0: 
      call add_target(meldef_ecclg0,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_ecclg0,op_ecclg,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_ecclg0
      labels(2) = op_ecclg
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0)
      call set_rule(meldef_ecclg0,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! OMG_T-list definition
      call add_target(meldef_omg_t,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_omg_t,op_omg_t,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_omg_t
      labels(2) = op_omg_t
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0)
      call set_rule(meldef_omg_t,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! OMG_L-list definition
      call add_target(meldef_omg_l,ttype_opme,.false.,tgt_info)
      call set_dependency(meldef_omg_l,op_omg_l,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_omg_l
      labels(2) = op_omg_l
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0)
      call set_rule(meldef_omg_l,ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c      ! Hbar definition
c      call add_target(meldef_hbar,ttype_opme,.false.,tgt_info)
c      call set_dependency(meldef_hbar,op_hbar,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = mel_hbar
c      labels(2) = op_hbar
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0)
c      call set_rule(meldef_hbar,ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)

      
*----------------------------------------------------------------------*
*     "phony" targets 
*----------------------------------------------------------------------*
      ! totally symmetric dia for use below:
      call me_list_label(mel_dia1,mel_dia,1,0,0,0,.false.)

      ! solve GS equations
      call add_target(solve_cc_gs,ttype_gen,.true.,tgt_info)
      call set_dependency(solve_cc_gs,mel_dia1,tgt_info)
      call set_dependency(solve_cc_gs,fopt_eccrs0,tgt_info)
      call solve_parameters(-1,parameters,2, 2,1,'DIA/DIA')
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = mel_top
      labels(2) = mel_tbar
      labels(3) = mel_omg_t
      labels(4) = mel_omg_l
      labels(5) = mel_dia1
      labels(6) = mel_dia1
      labels(7) = mel_ecclg0
      labels(8) = fopt_eccrs0
      call set_rule(solve_cc_gs,ttype_opme,SOLVENLEQ,
     &     labels,8,4,
     &     parameters,2,tgt_info)

      return
      end
