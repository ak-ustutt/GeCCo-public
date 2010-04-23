*----------------------------------------------------------------------*
      subroutine set_unc_mrci_targets(tgt_info,orb_info,calc)
*----------------------------------------------------------------------*
*     set targets for CASSCF or uncontracted CI coefficient optimization
*
*     matthias, fall 2009
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

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      logical, intent(in) ::
     &     calc

      integer ::
     &     ndef, occ_def(ngastp,2,60),
     &     isym, msc, ip, ih, 
     &     cminh, cmaxh, cminp, cmaxp, cmaxexc
      character(len_target_name) ::
     &     dia_label, labels(20)
      character(len_command_par) ::
     &     parameters(3)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting targets for multiref. wave function'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('calculate.multiref','cminh',
     &     ival=cminh)
      call get_argument_value('calculate.multiref','cmaxh',
     &     ival=cmaxh)
      call get_argument_value('calculate.multiref','cminp',
     &     ival=cminp)
      call get_argument_value('calculate.multiref','cmaxp',
     &     ival=cmaxp)
      call get_argument_value('calculate.multiref','cmaxexc',
     &     ival=cmaxexc)
      if (cmaxh.lt.0) cmaxh = cmaxexc
      if (cmaxp.lt.0) cmaxp = cmaxexc

      if (ntest.ge.100) then
        write(luout,*) 'cminh   = ',cminh
        write(luout,*) 'cmaxh   = ',cmaxh
        write(luout,*) 'cminp   = ',cminp
        write(luout,*) 'cmaxp   = ',cmaxp
        write(luout,*) 'cmaxexc = ',cmaxexc
        write(luout,*) 'nactel  = ',orb_info%nactel
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define scalar reference energy
      call add_target('E_REF',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('E_REF',ttype_op,DEF_HAMILTONIAN,'E_REF',
     &              1,1,parameters,1,tgt_info)

      ! define active electron creation operator C0
      call add_target('C0',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = cminp, cmaxp
        do ih = cminh, cmaxh
          if (orb_info%nactel+ih-ip.lt.0) cycle
          ndef = ndef + 1
          occ_def(IHOLE,2,ndef) = ih
          occ_def(IPART,1,ndef) = ip
          occ_def(IVALE,1,ndef) = orb_info%nactel + ih - ip
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('C0',ttype_op,DEF_OP_FROM_OCC,
     &              'C0',1,1,
     &              parameters,2,tgt_info)
 
      ! define product Jacobian times C0
      call add_target('A_C0',ttype_op,.false.,tgt_info)
      call set_dependency('A_C0','C0',tgt_info)
      call cloneop_parameters(-1,parameters,'C0',.false.)
      call set_rule('A_C0',ttype_op,CLONE_OP,'A_C0',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner for reference
      call add_target(op_dia//'_C0',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_C0','C0',tgt_info)
      call cloneop_parameters(-1,parameters,'C0',.false.)
      call set_rule(op_dia//'_C0',ttype_op,CLONE_OP,op_dia//'_C0',1,1,
     &              parameters,1,tgt_info)

      ! define Fock operator wrt reference function
      call add_target('FREF',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,1,1,.false.)
      call set_rule('FREF',ttype_op,DEF_HAMILTONIAN,'FREF',
     &              1,1,parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! reference energy expression
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_REF'
      labels(2) = 'E_REF'
      labels(3) = 'C0^+'
      labels(4) = op_ham
      labels(5) = 'C0'
      call add_target('F_REF',ttype_frm,.false.,tgt_info)
      call set_dependency('F_REF','E_REF',tgt_info)
      call set_dependency('F_REF',op_ham,tgt_info)
      call set_dependency('F_REF','C0',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'reference energy expression',3,
     &     (/2,3,4/),
     &     (/-1,-1,-1/),
     &     (/-1,-1,-1/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('F_REF',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,5,1,
     &              parameters,3,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_REF',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! hamiltonian times casscf coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'F_A_C0'
      labels(2) = 'F_REF'
      labels(3) = 'A_C0'
      labels(4) = 'C0^+'
      labels(5) = ' '
      call add_target('F_A_C0',ttype_frm,.false.,tgt_info)
      call set_dependency('F_A_C0','F_REF',tgt_info)
      call set_dependency('F_A_C0','A_C0',tgt_info)
      call set_dependency('F_A_C0','C0',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'CI matrix times right hand vector',1,'---')
      call set_rule('F_A_C0',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_A_C0',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! Fock operator wrt reference function
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_FREF'
      labels(2) = 'FREF'
      labels(3) = 'FREF'
      labels(4) = 'C0^+'
      labels(5) = op_ham
      labels(6) = 'C0'
      labels(7) = 'FREF'
      call add_target('F_FREF',ttype_frm,.false.,tgt_info)
      call set_dependency('F_FREF','FREF',tgt_info)
      call set_dependency('F_FREF',op_ham,tgt_info)
      call set_dependency('F_FREF','C0',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'reference energy expression',5,
     &     (/1,2,3,4,1/),
     &     (/-1,-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1,-1/),
     &     0,0,
     &     (/1,4,2,5/),2,
     &     0,0)
      call set_rule('F_FREF',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,7,1,
     &              parameters,3,tgt_info)
      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
      call set_rule('F_FREF',ttype_frm,PRINT_FORMULA,
     &                labels,2,1,parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! hamiltonian times casscf coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_A_C0'
      labels(2) = 'F_A_C0'
      call add_target('FOPT_A_C0',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_A_C0','F_A_C0',tgt_info)
      call set_dependency('FOPT_A_C0',mel_ham,tgt_info)
      call set_dependency('FOPT_A_C0','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_A_C0','DEF_ME_A_C0',tgt_info)
      call set_dependency('FOPT_A_C0','DEF_ME_E_REF',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_A_C0',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! reference energy expression
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_REF'
      labels(2) = 'F_REF'
      call add_target('FOPT_REF',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_REF','F_REF',tgt_info)
      call set_dependency('FOPT_REF',mel_ham,tgt_info)
      call set_dependency('FOPT_REF','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_REF','DEF_ME_E_REF',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_REF',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Fock operator wrt reference function
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_FREF'
      labels(2) = 'F_FREF'
      call add_target('FOPT_FREF',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_FREF','F_FREF',tgt_info)
      call set_dependency('FOPT_FREF',mel_ham,tgt_info)
      call set_dependency('FOPT_FREF','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_FREF','DEF_ME_FREF',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_FREF',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME_C0
      call add_target('DEF_ME_C0',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_C0','C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_C0'
      labels(2) = 'C0'
      call me_list_parameters(-1,parameters,
     &     msc,0,orb_info%lsym,
     &     0,0,.false.)
      call set_rule('DEF_ME_C0',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_A_C0
      call add_target('DEF_ME_A_C0',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_A_C0','A_C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_A_C0'
      labels(2) = 'A_C0'
      call me_list_parameters(-1,parameters,
     &     msc,0,orb_info%lsym,
     &     0,0,.false.)
      call set_rule('DEF_ME_A_C0',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_E_REF
      call add_target('DEF_ME_E_REF',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_E_REF','E_REF',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_E_REF'
      labels(2) = 'E_REF'
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DEF_ME_E_REF',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner for reference
      call me_list_label(dia_label,mel_dia,orb_info%lsym,
     &     0,0,0,.false.)
      call add_target(trim(dia_label)//'C0',ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label)//'C0',mel_ham,tgt_info)
      call set_dependency(trim(dia_label)//'C0',
     &                    op_dia//'_'//'C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C0'
      labels(2) = op_dia//'_'//'C0'
      call me_list_parameters(-1,parameters,
     &     0,0,orb_info%lsym,
     &     0,0,.false.)
      call set_rule(trim(dia_label)//'C0',ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      labels(1) = trim(dia_label)//'C0'
      labels(2) = mel_ham
      call set_rule(trim(dia_label)//'C0',ttype_opme,
     &              PRECONDITIONER,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! ME_FREF
      call add_target('DEF_ME_FREF',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_FREF','FREF',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_FREF'
      labels(2) = 'FREF'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_FREF',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! SOLVE reference eigenvalue equation
      call add_target('SOLVE_REF',ttype_gen,.false.,tgt_info)
      call set_dependency('SOLVE_REF','FOPT_A_C0',tgt_info)
      call me_list_label(dia_label,mel_dia,orb_info%lsym,
     &     0,0,0,.false.)
      call set_dependency('SOLVE_REF',trim(dia_label)//'C0',tgt_info)
      call solve_parameters(-1,parameters,2,1,1,'DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_C0'
      labels(2) = trim(dia_label)//'C0'
      labels(3) = 'A_C0'
      labels(4) = 'C0'
      labels(5) = 'FOPT_A_C0'
      call set_rule('SOLVE_REF',ttype_opme,SOLVEEVP,
     &     labels,5,1,
     &     parameters,2,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'CI coefficients :',0,'LIST')
      call set_rule('SOLVE_REF',ttype_opme,PRINT_MEL,
     &     'ME_C0',1,0,
     &     parameters,2,tgt_info)

      ! Evaluate reference energy
      call add_target('EVAL_E_REF',ttype_gen,calc,tgt_info)
      call set_dependency('EVAL_E_REF','SOLVE_REF',tgt_info)
      call set_dependency('EVAL_E_REF','FOPT_REF',tgt_info)
      call set_rule('EVAL_E_REF',ttype_opme,EVAL,
     &     'FOPT_REF',1,0,
     &     parameters,0,tgt_info)

      ! Evaluate Fock operator wrt reference function
      call add_target('EVAL_FREF',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_FREF','FOPT_FREF',tgt_info)
      call set_dependency('EVAL_FREF','SOLVE_REF',tgt_info)
      call set_rule('EVAL_FREF',ttype_opme,EVAL,
     &     'FOPT_FREF',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'effective Fock operator:',0,'LIST')
c      call set_rule('EVAL_FREF',ttype_opme,PRINT_MEL,
c     &     'ME_FREF',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      return
      end
