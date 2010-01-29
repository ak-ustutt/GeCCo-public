*----------------------------------------------------------------------*
      subroutine set_r12_test_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
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
      character(*), intent(in) ::
     &     env_type

      integer ::
     &     ndef,
     &     ncat, nint, ncntrtest, nfgentest, icnt,
     &     isym, ms, msc, sym_arr(8), 
     &     occ_def(ngastp,2,20)
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(12) ::
     &     approx, F_appr, K_appr
      integer, pointer ::
     &     testnr(:)


      integer, external ::
     &     idxlist

*----------------------------------------------------------------------*
      if (is_keyword_set('calculate.check').eq.0) return

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting test targets for R12 ...'

      msc = +1  ! assuming closed shell
*----------------------------------------------------------------------*
*     read input
*----------------------------------------------------------------------*
      ncntrtest = is_argument_set('calculate.check','contr_test')
      allocate(testnr(ncntrtest))
      do icnt = 1, ncntrtest
        call get_argument_value('calculate.check','contr_test',
     &       argcount=icnt,ival=testnr(icnt))
      end do
      nfgentest = is_argument_set('calculate.check','algebra_test')

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      call add_target('TEST_RES',ttype_op,.false.,tgt_info)
      call set_rule('TEST_RES',ttype_op,DEF_SCALAR,
     &              'TEST_RES',1,1,
     &              parameters,0,tgt_info)

      occ_def = 0
      ndef = 3
      occ_def(IPART,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IPART,1,2) = 2
      occ_def(IHOLE,2,2) = 2
      occ_def(IPART,1,3) = 1
      occ_def(IEXTR,1,3) = 1
      occ_def(IHOLE,2,3) = 2
      call add_target('TEST_R',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_R',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_R',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 2
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,2) = 2
      occ_def(IPART,1,3) = 2
      occ_def(IHOLE,2,4) = 2
      call add_target('TEST_R2',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/  0,   0,  0,   0/),ndef)
      call set_rule('TEST_R2',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_R2',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 1
      occ_def(IPART,1,1) = 2
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,1) = 3
      call add_target('TEST_RT',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_RT',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_RT',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 2
      occ_def(IEXTR,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 2
      occ_def(IEXTR,1,2) = 1
      occ_def(IPART,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      call add_target('TEST_R3',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_R3',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_R3',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 2
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,2) = 2
      occ_def(IPART,1,3) = 2
      occ_def(IHOLE,2,4) = 2
      call add_target('TEST_H2',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/  0,   0,  0,   0/),ndef)
      call set_rule('TEST_H2',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_H2',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 3
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IPART,2,1) = 2
      occ_def(IHOLE,1,2) = 1
      occ_def(IPART,1,2) = 1
      occ_def(IHOLE,2,2) = 1
      occ_def(IPART,2,2) = 1
      occ_def(IHOLE,1,3) = 1
      occ_def(IPART,1,3) = 1
      occ_def(IPART,2,3) = 2
      call add_target('TEST_H',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_H',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_H',1,1,
     &              parameters,2,tgt_info)

      ! prototype singles and doubles ex.
      occ_def = 0
      ndef = 2
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,1,2) = 2
      occ_def(IHOLE,2,2) = 2
      call add_target('TEST_A',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_A',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_A',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 2
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,1,2) = 2
      occ_def(IHOLE,2,2) = 2
      call add_target('TEST_B',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule('TEST_B',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_B',1,1,
     &              parameters,2,tgt_info)

      
*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*

      ! test formula I
      call add_target('FORM_TEST_I',ttype_frm,.false.,tgt_info)
      call set_dependency('FORM_TEST_I','TEST_RES',tgt_info)
      call set_dependency('FORM_TEST_I','TEST_A',tgt_info)
      call set_dependency('FORM_TEST_I','TEST_R',tgt_info)
      call set_dependency('FORM_TEST_I','TEST_H',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      labels(2) = 'TEST_RES'
      labels(3) = 'TEST_A^+'
      labels(4) = 'TEST_R^+'
      labels(5) = 'TEST_H'
      labels(6) = 'TEST_A'
      labels(7) = 'TEST_A'
      labels(8) = 'TEST_A'
      call expand_parameters(-1,
     &     parameters,3,
     &     'formula test 2',6,
     &     (/1,2,3,4,5,6/),
     &     (/1,1,1,1,1,1/),
     &     (/1,1,1,1,1,1/),
c     &     0,0,
     &     (/1,6, 3,4, 2,5, 3,5/),4,
     &     0,0,
     &     (/1,2,1,IPART, 2,3,1,IEXTR/),2)
      call set_rule('FORM_TEST_I',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,8,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_I',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! test formula II
      call add_target('FORM_TEST_II',ttype_frm,.false.,tgt_info)
      call set_dependency('FORM_TEST_II','TEST_RES',tgt_info)
      call set_dependency('FORM_TEST_II','TEST_A',tgt_info)
      call set_dependency('FORM_TEST_II','TEST_B',tgt_info)
      call set_dependency('FORM_TEST_II','TEST_R',tgt_info)
      call set_dependency('FORM_TEST_II','TEST_H',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_II'
      labels(2) = 'TEST_RES'
      labels(3) = 'TEST_A^+'
      labels(4) = 'TEST_R^+'
      labels(5) = 'TEST_H'
      labels(6) = 'TEST_A'
      labels(7) = 'TEST_B'
      labels(8) = 'TEST_B'
      call expand_parameters(-1,
     &     parameters,3,
     &     'formula test 2',6,
     &     (/1,2,3,4,5,6/),
     &     (/1,1,1,1,1,1/),
     &     (/1,1,1,1,1,1/),
     &     (/1,2, 2,3, 1,4, 3,4, 2,5, 3,5/),6,
     &     0,0,
     &     0,0)
      call set_rule('FORM_TEST_II',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,8,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_II',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! test formula III
      call add_target('FORM_TEST_III',ttype_frm,.false.,tgt_info)
      call set_dependency('FORM_TEST_III','TEST_RES',tgt_info)
      call set_dependency('FORM_TEST_III','TEST_A',tgt_info)
      call set_dependency('FORM_TEST_III','TEST_B',tgt_info)
      call set_dependency('FORM_TEST_III','TEST_R',tgt_info)
      call set_dependency('FORM_TEST_III','TEST_H',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_III'
      labels(2) = 'TEST_RES'
      labels(3) = 'TEST_A^+'
      labels(4) = 'TEST_H'
      labels(5) = 'TEST_B'
      call expand_parameters(-1,
     &     parameters,3,
     &     'formula test 3',3,
     &     (/1,2,3/),
     &     (/2,2,2/),
     &     (/2,2,2/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('FORM_TEST_III',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,5,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_III',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! formula generation test: RT = R.T
      call add_target('FORM_TEST_RT',ttype_frm,.false.,tgt_info)
      call set_dependency('FORM_TEST_RT','TEST_RT',tgt_info)
      call set_dependency('FORM_TEST_RT','TEST_A',tgt_info)
      call set_dependency('FORM_TEST_RT','TEST_R',tgt_info)
      call set_dependency('FORM_TEST_RT','TEST_H',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_RT'
      labels(2) = 'TEST_RT'
      labels(3) = 'TEST_RT'
      labels(4) = 'TEST_R'
      labels(5) = 'TEST_A'
      labels(6) = 'TEST_RT'
      call expand_parameters(-1,
     &     parameters,3,
     &     'formula test RT',4,
     &     (/1,2,3,1/),
     &     (/1,1,2,1/),
     &     (/1,1,2,1/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('FORM_TEST_RT',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,6,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_RT',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! formula generation test: R+.H.RT + expand RT
      needed = nfgentest.gt.0
      call add_target('FORM_TEST_RHR',ttype_frm,needed,tgt_info)
      call set_dependency('FORM_TEST_RHR','TEST_RES',tgt_info)
      call set_dependency('FORM_TEST_RHR','TEST_RT',tgt_info)
      call set_dependency('FORM_TEST_RHR','TEST_R',tgt_info)
      call set_dependency('FORM_TEST_RHR','TEST_H',tgt_info)
      call set_dependency('FORM_TEST_RHR','FORM_TEST_RT',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_RHR'
      labels(2) = 'TEST_RES'
      labels(3) = 'TEST_R^+'
      labels(4) = 'TEST_H'
      labels(5) = 'TEST_RT'
      call expand_parameters(-1,
     &     parameters,3,
     &     'formula test RHR',3,
     &     (/1,2,3/),
     &     (/3,3,1/),
     &     (/3,3,1/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('FORM_TEST_RHR',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,5,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_RHR',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_RHRX'
      labels(2) = 'FORM_TEST_RHR'
      labels(3) = 'FORM_TEST_RT'
      call form_parameters(-1,
     &     parameters,2,'Expanded',1,'---')
      call set_rule('FORM_TEST_RHR',ttype_frm,EXPAND,
     &              labels,3,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('FORM_TEST_RHR',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)
      
*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! test I
      ! (1)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.1',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.1','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,2,(/1,0/),2) !just for safety
      call set_rule('FOPT_TEST_I.1',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.1'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.1','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.1','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.1','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.1','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.1',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (2)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.2',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.2','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,3,5,2,2,1/),7)
      call set_rule('FOPT_TEST_I.2',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.2'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.2','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.2','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.2','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.2','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.2',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (3)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.3',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.3','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,3,6,5,2,1/),7)
      call set_rule('FOPT_TEST_I.3',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.3'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.3','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.3','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.3','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.3','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.3',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (4)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.4',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.4','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,1,1,1,2,1/),7)
      call set_rule('FOPT_TEST_I.4',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.4'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.4','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.4','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.4','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.4','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.4',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (5)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.5',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.5','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,7,5,1,1,1/),7)
      call set_rule('FOPT_TEST_I.5',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.5'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.5','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.5','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.5','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.5','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.5',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (6)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.6',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.6','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,7,6,3,2,1/),7)
      call set_rule('FOPT_TEST_I.6',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.6'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.6','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.6','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.6','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.6','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.6',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (7)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_I'
      call add_target('FOPT_TEST_I.7',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I.7','FORM_TEST_I',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,5,6,4,3,1/),7)
      call set_rule('FOPT_TEST_I.7',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_I.7'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_I.7','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I.7','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I.7','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I.7','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I.7',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! test II
      ! (1)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_II'
      call add_target('FOPT_TEST_II.1',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_II.1','FORM_TEST_II',tgt_info)
      call modify_parameters(-1,parameters,2,(/1,0/),2) !just for safety
      call set_rule('FOPT_TEST_II.1',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_II.1'
      labels(2) = 'FORM_TEST_II'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_II.1','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_II.1','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_II.1','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_II.1','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_II.1','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_II.1',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (2)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_II'
      call add_target('FOPT_TEST_II.2',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_II.2','FORM_TEST_II',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,3,4,6,2,1/),7)
      call set_rule('FOPT_TEST_II.2',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_II.2'
      labels(2) = 'FORM_TEST_II'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_II.2','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_II.2','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_II.2','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_II.2','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_II.2','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_II.2',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (3)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_II'
      call add_target('FOPT_TEST_II.3',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_II.3','FORM_TEST_II',tgt_info)
      call modify_parameters(-1,parameters,7,(/1,5,8,5,4,1,1/),7)
      call set_rule('FOPT_TEST_II.3',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_II.3'
      labels(2) = 'FORM_TEST_II'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_II.3','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_II.3','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_II.3','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_II.3','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_II.3','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_II.3',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! test III
      ! (1)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_III'
      call add_target('FOPT_TEST_III.1',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_III.1','FORM_TEST_III',tgt_info)
      call modify_parameters(-1,parameters,4,(/1,2,1,1/),4)
      call set_rule('FOPT_TEST_III.1',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_III.1'
      labels(2) = 'FORM_TEST_III'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_III.1','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_III.1','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_III.1','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_III.1','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_III.1','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_III.1',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (2)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_III'
      call add_target('FOPT_TEST_III.2',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_III.2','FORM_TEST_III',tgt_info)
      call modify_parameters(-1,parameters,4,(/1,2,2,2/),4)
      call set_rule('FOPT_TEST_III.2',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_III.2'
      labels(2) = 'FORM_TEST_III'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_III.2','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_III.2','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_III.2','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_III.2','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_III.2','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_III.2',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
      ! (3)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_TEST_III'
      call add_target('FOPT_TEST_III.3',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_III.3','FORM_TEST_III',tgt_info)
      call modify_parameters(-1,parameters,4,(/1,2,3,1/),4)
      call set_rule('FOPT_TEST_III.3',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      labels(1) = 'FOPT_TEST_III.3'
      labels(2) = 'FORM_TEST_III'
      ncat = 1
      nint = 0
      call set_dependency('FOPT_TEST_III.3','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_III.3','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_III.3','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_III.3','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_III.3','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_III.3',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      call add_target('DEF-ME_TEST_RES',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-ME_TEST_RES','TEST_RES',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_RES'
      labels(2) = 'TEST_RES'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-ME_TEST_RES',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target('ME_TEST_R',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_R','TEST_R',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_R'
      labels(2) = 'TEST_R'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_R',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_R'
      call import_parameters(-1,parameters,'F12_INT',env_type)
      call set_rule('ME_TEST_R',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('ME_TEST_R2',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_R2','TEST_R2',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_R2'
      labels(2) = 'TEST_R2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_R2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_R2'
      call import_parameters(-1,parameters,'F12_INT',env_type)
      call set_rule('ME_TEST_R2',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      ! 
      call add_target('ME_TEST_H',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_H','TEST_H',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_H'
      labels(2) = 'TEST_H'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_H',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_H'
      call import_parameters(-1,parameters,'G_INT',env_type)
      call set_rule('ME_TEST_H',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('ME_TEST_H2',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_H2','TEST_H2',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_H2'
      labels(2) = 'TEST_H2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_H2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      ! (b) import
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_H2'
      call import_parameters(-1,parameters,'G_INT',env_type)
      call set_rule('ME_TEST_H2',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('ME_TEST_A',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_A','TEST_A',tgt_info)
      call set_dependency('ME_TEST_A','ME_TEST_R2',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_A'
      labels(2) = 'TEST_A'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_A',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_A'
      labels(2) = 'TEST_A'
      labels(3) = 'TEST_A'
      labels(4) = 'TEST_R2'
      labels(5) = 'TEST_R2'
      labels(6) = 'TEST_A'
      call expand_parameters(-1,
     &     parameters,3,
     &     'A from R',4,
     &     (/1,2,2,1/),
     &     (/-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('ME_TEST_A',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,6,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('ME_TEST_A',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FOPT_A'
      labels(2) = 'FORM_A'
      call opt_parameters(-1,parameters,1,0)
      call set_rule('ME_TEST_A',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_rule('ME_TEST_A',ttype_opme,EVAL,
     &     'FOPT_A',1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST OPERATOR A: norm =',0,'NORM E20.10')
      call set_rule('ME_TEST_A',ttype_opme,PRINT_MEL,
     &     'ME_TEST_A',1,0,
     &     parameters,2,tgt_info)

      call add_target('ME_TEST_B',ttype_opme,.false.,tgt_info)
      call set_dependency('ME_TEST_B','TEST_B',tgt_info)
      call set_dependency('ME_TEST_B','ME_TEST_H2',tgt_info)
      ! (a) define
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_TEST_B'
      labels(2) = 'TEST_B'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('ME_TEST_B',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FORM_B'
      labels(2) = 'TEST_B'
      labels(3) = 'TEST_B'
      labels(4) = 'TEST_H2'
      labels(5) = 'TEST_H2'
      labels(6) = 'TEST_B'
      call expand_parameters(-1,
     &     parameters,3,
     &     'B from H',4,
     &     (/1,2,2,1/),
     &     (/-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1/),
     &     0,0,
     &     0,0,
     &     0,0)
      call set_rule('ME_TEST_B',ttype_frm,EXPAND_OP_PRODUCT,
     &     labels,6,1,
     &     parameters,3,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('ME_TEST_B',ttype_frm,PRINT_FORMULA,
     &              labels,2,1,
     &              parameters,2,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FOPT_B'
      labels(2) = 'FORM_B'
      call opt_parameters(-1,parameters,1,0)
      call set_rule('ME_TEST_B',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_rule('ME_TEST_B',ttype_opme,EVAL,
     &     'FOPT_B',1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST OPERATOR B: norm =',0,'NORM E20.10')
      call set_rule('ME_TEST_B',ttype_opme,PRINT_MEL,
     &     'ME_TEST_B',1,0,
     &     parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      needed = idxlist(1,testnr,ncntrtest,1).gt.0
      
      call add_target('CONTR_TEST_I',ttype_gen,needed,tgt_info)
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.1',tgt_info)
      labels(1) = 'FOPT_TEST_I.1'
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.1, opt. factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.2',tgt_info)
      labels(1) = 'FOPT_TEST_I.2'
c      call modify_parameters(-1,parameters,7,(/1,5,3,5,2,2,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.2, 35221 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.3',tgt_info)
      labels(1) = 'FOPT_TEST_I.3'
c      call modify_parameters(-1,parameters,7,(/1,5,3,6,5,2,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.3, 36521 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.4',tgt_info)
      labels(1) = 'FOPT_TEST_I.4'
c      call modify_parameters(-1,parameters,7,(/1,5,1,1,1,2,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.4, 11121 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.5',tgt_info)
      labels(1) = 'FOPT_TEST_I.5'
c      call modify_parameters(-1,parameters,7,(/1,5,7,5,1,1,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.5, 75111 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.6',tgt_info)
      labels(1) = 'FOPT_TEST_I.6'
c      call modify_parameters(-1,parameters,7,(/1,5,7,6,3,2,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.6, 76321 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      call set_dependency('CONTR_TEST_I','FOPT_TEST_I.7',tgt_info)
      labels(1) = 'FOPT_TEST_I.7'
c      call modify_parameters(-1,parameters,7,(/1,5,5,6,4,3,1/),7)
c      call set_rule('CONTR_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST I.7, 56431 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_I',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! test 2
      needed = idxlist(2,testnr,ncntrtest,1).gt.0
      
      call add_target('CONTR_TEST_II',ttype_gen,needed,tgt_info)
      call set_dependency('CONTR_TEST_II','FOPT_TEST_II.1',tgt_info)
      labels(1) = 'FOPT_TEST_II.1'
      call set_rule('CONTR_TEST_II',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST II.1, opt. factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_II',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_II',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call set_dependency('CONTR_TEST_II','FOPT_TEST_II.2',tgt_info)
      labels(1) = 'FOPT_TEST_II.2'
c      call modify_parameters(-1,parameters,7,(/1,5,3,4,6,2,1/),7)
c      call set_rule('CONTR_TEST_II',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_II',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST II.2, 34621 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_II',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_II',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call set_dependency('CONTR_TEST_II','FOPT_TEST_II.3',tgt_info)
      labels(1) = 'FOPT_TEST_II.3'
c      call modify_parameters(-1,parameters,7,(/1,5,8,5,4,1,1/),7)
c      call set_rule('CONTR_TEST_II',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_II',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST II.3, 85411 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_II',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_II',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! test 3
      needed = idxlist(3,testnr,ncntrtest,1).gt.0
      
      call add_target('CONTR_TEST_III',ttype_gen,needed,tgt_info)
      call set_dependency('CONTR_TEST_III','FOPT_TEST_III.1',tgt_info)
      labels(1) = 'FOPT_TEST_III.1'
c      call modify_parameters(-1,parameters,4,(/1,2,1,1/),4)
c      call set_rule('CONTR_TEST_III',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_III',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST III.1, 12 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_III',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_III',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call set_dependency('CONTR_TEST_III','FOPT_TEST_III.2',tgt_info)
      labels(1) = 'FOPT_TEST_III.2'
c      call modify_parameters(-1,parameters,4,(/1,2,2,2/),4)
c      call set_rule('CONTR_TEST_III',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_III',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST III.2, 23 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_III',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_III',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      call set_dependency('CONTR_TEST_III','FOPT_TEST_III.3',tgt_info)
      labels(1) = 'FOPT_TEST_III.3'
c      call modify_parameters(-1,parameters,4,(/1,2,3,1/),4)
c      call set_rule('CONTR_TEST_III',ttype_frm,MODIFY_FACTORIZATION,
c     &     labels,1,0,
c     &     parameters,1,tgt_info)
      call set_rule('CONTR_TEST_III',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     'TEST III.3, 31 factorization:',0,'SCAL E20.10')
      call set_rule('CONTR_TEST_III',ttype_opme,PRINT_MEL,
     &     'ME_TEST_RES',1,0,
     &     parameters,2,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('CONTR_TEST_III',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)


      deallocate(testnr)
      
      return

      end
