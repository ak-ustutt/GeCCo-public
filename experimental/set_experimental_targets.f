*----------------------------------------------------------------------*
      subroutine set_experimental_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets for experiments with GeCCo
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
     &     labels(20)
      character(len_command_par) ::
     &     parameters(2)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting experimental targets ...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! e.g. 
c      call add_target('MY_OP',ttype_op,.false.,tgt_info)
c      call set_rule('MY_OP',ttype_op,DEF_SCALAR,
c     &              'MY_OP',1,1,
c     &              parameters,0,tgt_info)
      ! cf. set_xxxx_targets.f routines in e.g. cc_special for
      ! further examples

      call add_target('H_0',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,1,2,1,.false.)
      call set_rule('H_0',ttype_op,DEF_HAMILTONIAN,'H_0',
     &              1,1,parameters,1,tgt_info)
      call ord_parameters(-1,parameters,0)
      call set_rule('H_0',ttype_op,SET_ORDER,'H_0',
     &              1,1,parameters,1,tgt_info)
c      call add_target('PHI',ttype_op,.false.,tgt_info)
c      call hop_parameters(-1,parameters,2,2,1,.false.)
c      call set_rule('PHI',ttype_op,DEF_HAMILTONIAN,'PHI',
c     &              1,1,parameters,1,tgt_info)
      call add_target('LRESP',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('LRESP',ttype_op,DEF_HAMILTONIAN,'LRESP',
     &              1,1,parameters,1,tgt_info)
      call add_target('T2',ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,.false.,2,2,0,1)
      call set_rule('T2',ttype_op,DEF_EXCITATION,'T2',
     &              1,1,parameters,1,tgt_info)
c      call add_target('O2',ttype_op,.false.,tgt_info)
c      call set_dependency('O2','T2',tgt_info)
c      call cloneop_parameters(-1,parameters,'T2',.false.)
c      call set_rule('O2',ttype_op,CLONE_OP,'O2',1,1,
c     &              parameters,1,tgt_info)
      
*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*
      ! e.g.
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RESP_LAG'
      labels(2) = 'LRESP'
      labels(3) = 'H_0'
      labels(4) = 'T2'
      call add_target('RESP_LAG',ttype_frm,.false.,tgt_info)
      call set_dependency('RESP_LAG','LRESP',tgt_info)
      call set_dependency('RESP_LAG','H_0',tgt_info)
      call set_dependency('RESP_LAG','T2',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'response lagrange functional',0,'---')
      call set_rule('RESP_LAG',ttype_frm,DEF_EXP_FORMULA,
     &              labels,4,1,
     &              parameters,2,tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'MP2_FORM'
c      labels(2) = 'MP2_FORM'
c      labels(3) = 'FOCK'
c      labels(4) = 'H'
c      labels(5) = 'PHI'
c      labels(6) = 'H' 
c      call set_dependency('MP2_FORM','H',tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'MP2 lagrange functional',2,'---')
c      call set_rule('MP2_FORM',ttype_frm,REPLACE,
c     &     labels,6,1,
c     &     parameters,2,tgt_info)

c      labels(1:20)(1:len_target_name)= ' '
c      labels(1) = 'MP2_RES'
c      labels(2) = 'MP2_FORM'
c      labels(3) = 'O2'
c      labels(4) = 'T2^+'
c      labels(5) = ' '
c      call add_target('MP2_RES',ttype_frm,.false.,tgt_info)
c      call set_dependency('MP2_RES','MP2_FORM',tgt_info)
c      call set_dependency('MP2_RES','O2',tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'residual',1,'---')
c      call set_rule('MP2_RES',ttype_frm,DERIVATIVE,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)

c      labels(1:20)(1:len_target_name)= ' '
c      labels(1) = 'MP2_EN'
c      labels(2) = 'MP2_FORM'
c      labels(3) = 'MP2_RES'
c      call add_target('MP2_EN',ttype_frm,.false.,tgt_info)
c      call set_dependency('MP2_EN','MP2_FORM',tgt_info)
c      call set_dependency('MP2_EN','MP2_RES',tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'energy',1,'---')
c      call set_rule('MP2_EN',ttype_frm,FACTOR_OUT,
c     &              labels,3,1,
c     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! e.g.:
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'MP2_OPT'
c      labels(2) = 'MP2_RES'
c      labels(3) = 'MP2_EN'
c      ncat = 2  ! 2 formulae pasted into final formula
c      nint = 0  ! no intermediate to factor out so far ...
c      call add_target('MP2_OPT',ttype_frm,.false.,tgt_info)
c      call set_dependency('MP2_OPT','MP2_RES',tgt_info)
c      call set_dependency('MP2_OPT','MP2_EN',tgt_info)
c      call set_dependency('MP2_OPT','DEF_ME_LMP2',tgt_info)
c      call set_dependency('MP2_OPT','DEF_ME_O2',tgt_info)
c      call set_dependency('MP2_OPT','DEF_ME_T2',tgt_info)
c      call set_dependency('MP2_OPT','H0',tgt_info)      
c      call opt_parameters(-1,parameters,ncat,nint)
c      call set_rule('MP2_OPT',ttype_frm,OPTIMIZE,
c     &              labels,ncat+nint+1,1,
c     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! e.g.: 
c      call add_target('DEF_ME_LMP2',ttype_opme,.false.,tgt_info)
c      call set_dependency('DEF_ME_LMP2','LMP2',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'ME_LMP2'
c      labels(2) = 'LMP2'
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule('DEF_ME_LMP2',ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      
c      call add_target('DEF_ME_T2',ttype_opme,.false.,tgt_info)
c      call set_dependency('DEF_ME_T2','T2',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'ME_T2'
c      labels(2) = 'T2'
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule('DEF_ME_T2',ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)

c      call add_target('DEF_ME_O2',ttype_opme,.false.,tgt_info)
c      call set_dependency('DEF_ME_O2','O2',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'ME_O2'
c      labels(2) = 'O2'
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule('DEF_ME_O2',ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)

c      call add_target('DIAG',ttype_opme,.false.,tgt_info)
c      call set_dependency('DIAG','H0',tgt_info)
c      call set_dependency('DIAG','T2',tgt_info)
c      call cloneop_parameters(-1,parameters,
c     &     'T2',.false.)
c      call set_rule('DIAG',ttype_op,CLONE_OP,
c     &     'D2',1,1,
c     &     parameters,1,tgt_info)
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = 'DIAG'
c      labels(2) = 'D2'
c      call me_list_parameters(-1,parameters,
c     &     0,0,1,0,0,.false.)
c      call set_rule('DIAG',ttype_opme,DEF_ME_LIST,
c     &     labels,2,1,
c     &     parameters,1,tgt_info)
c      labels(1) = 'DIAG'
c      labels(2) = 'H0'
c      call set_rule('DIAG',ttype_opme,PRECONDITIONER,
c     &              labels,2,1,
c     &              parameters,0,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*
   
      call add_target('MY_TARGET',ttype_gen,.true.,tgt_info)
      call set_dependency('MY_TARGET','RESP_LAG',tgt_info)
c      call set_dependency('SOLVE_MP2','DIAG',tgt_info)
c      call solve_parameters(-1,parameters,2, 1,1,'DIA')
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'ME_T2'
c      labels(2) = 'ME_O2'
c      labels(3) = 'DIAG'
c      labels(4) = 'ME_LMP2'
c      labels(5) = 'MP2_OPT'
c      call set_rule('SOLVE_MP2',ttype_opme,SOLVENLEQ,
c     &     labels,5,2,
c     &     parameters,2,tgt_info)

      return
      end
