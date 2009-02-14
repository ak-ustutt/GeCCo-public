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

      integer ::
     &     hpvx_constr(2,ngastp,2),
     &     gas_constr(2,orb_info%ngas,2,2)

      integer ::
     &     min_rank, max_rank, ansatz, n_pp, ndef,
     &     min_rank_tp, min_rank_tpp,
     &     isim, ncat, nint, icnt, nlab, irank, idef,
     &     isym, ms, msc, sym_arr(8), extend, r12op,
     &     occ_def(ngastp,2,20),
     &     ntp_min, ntp_max, ntpp_min, ntpp_max, testnr
      logical ::
     &     needed, r12fix, set_tp, set_tpp
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(12) ::
     &     approx, F_appr, K_appr

      character(*), intent(in) ::
     &     env_type

*----------------------------------------------------------------------*
      if (is_keyword_set('calculate.check').eq.0) return

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting test targets for R12 ...'

      msc = +1  ! assuming closed shell
*----------------------------------------------------------------------*
*     read input
*----------------------------------------------------------------------*
      call get_argument_value('calculate.check','sign_test',ival=testnr)

      ! set approx string
      approx(1:12) = ' '
      F_appr(1:12) = ' '
      K_appr(1:12) = ' '
      call get_argument_value('method.R12','ansatz',ival=ansatz)
      call get_argument_value('method.R12','approx',str=approx)
      call get_argument_value('method.R12','F_appr',str=F_appr)
      call get_argument_value('method.R12','K_appr',str=K_appr)
      call get_argument_value('method.R12','min_tp',ival=min_rank_tp)
      call get_argument_value('method.R12','min_tpp',ival=min_rank_tpp)
      call get_argument_value('method.R12','minexc',ival=min_rank)
      call get_argument_value('method.R12','maxexc',ival=max_rank)
      call get_argument_value('method.R12','fixed',lval=r12fix)
      call get_argument_value('method.R12','extend',ival=extend)
      call get_argument_value('method.R12','r12op',ival=r12op)

      n_pp = 0  ! number of particle-particle interaction in R12
      set_tp = .false.
      set_tpp = .false.
      select case(extend)
      case(1) 
        ! T1'
        set_tp = .true.
        ntp_min=1
        ntp_max=1
        n_pp=1
      case(2)
        ! T0'
        set_tp = .true.
        ntp_min=0
        ntp_max=0
        n_pp=0
      case(3)
        ! T0' + T1'
        set_tp = .true.
        ntp_min=0
        ntp_max=1
        n_pp=1
      case(4)
        ! T1' + T2' (for CC)
        set_tp = .true.
        ntp_min=1
        ntp_max=2
        n_pp=1
      case(5)
        ! T1' + T2'
        set_tp = .true.
        ntp_min=1
        ntp_max=2
        n_pp=2
      case(6)
        ! T1' 
        set_tp = .true.
        ntp_min=1
        ntp_max=1
        n_pp=2
      case default 
        set_tp = .false.
        ntp_min=0
        ntp_max=0
        n_pp=0
      end select

      if (r12op.gt.0.and.extend.gt.0)
     &     call quit(1,'set_r12f_general_targets',
     &     'use either r12op or extend')
      if (r12op.gt.0) then
        set_tp  = .false.
        set_tpp = .false.
        ntp_min=0
        ntp_max=0
        ntpp_min=0
        ntpp_max=0
      end if
      select case(r12op)
      case(1)
        ! T' operators (singly p-connected to R12)
        set_tp = .true.
        ntp_min=min_rank_tp
        ntp_max=max_rank-1
        n_pp=1
      case(2)
        ! T'' operators (doubly p-connected to R12)
        set_tpp = .true.
        ntpp_min=min_rank_tpp
        ntpp_max=max_rank
        n_pp=2
      case(3,4)
        ! T' + T'' operators
        set_tp = .true.
        ntp_min=min_rank_tp
        ntp_max=max_rank-1
        n_pp=1
        set_tpp = .true.
        ntpp_min=min_rank_tpp
        ntpp_max=max_rank
        n_pp=2
      end select

      ! assemble approx string
      select case(trim(F_appr))
      case('none')
        write(luout,*) 'no approximations wrt. Fock made'
      case('no_Z')
        write(luout,*) 'Z matrix omitted'
        approx(4:6) = 'noZ'
      case('GBC','EBC')
        write(luout,*)
     &  'GBC/EBC are currently only possible be supplying the'
        write(luout,*)
     &  'suitable integrals. Make that sure and restart w/o'
        write(luout,*)
     &  'GBC/EBC flag'
        call quit(0,'set_r12_general_targets','GBC/EBC?')
      case default
        call quit(0,'set_r12_general_targets',
     &       'F_appr unknown: "'//trim(F_appr)//'"')
      end select

      select case(trim(K_appr))
      case('none')
        write(luout,*) 'no approximations wrt. Xchange made'
      case('HY1')
        write(luout,*) 'Y contribution omitted'
        approx(8:10) = 'HY1'
      case('HY2')
        write(luout,*) 'Y contribution approx with 1 CABS index'
        approx(8:10) = 'HY2'
      case default
        call quit(0,'set_r12_general_targets',
     &       'K_appr unknown: "'//trim(K_appr)//'"')
      end select

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! adapt: better define test-operators here
      call add_target('TEST_RES',ttype_op,.false.,tgt_info)
      call set_rule('TEST_RES',ttype_op,DEF_SCALAR,
     &              'TEST_RES',1,1,
     &              parameters,0,tgt_info)

      occ_def = 0
      ndef = 2
      occ_def(IPART,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,2,1) = 1
      occ_def(IPART,1,2) = 2
      occ_def(IHOLE,2,2) = 2
      call add_target('TEST_R',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/.true.,.true./),ndef)
      call set_rule('TEST_R',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_R',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,2) = 2
      call add_target('TEST_R2',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/.true.,.true./),ndef)
      call set_rule('TEST_R2',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_R2',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,2) = 2
      call add_target('TEST_H2',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,2,(/.true.,.true./),ndef)
      call set_rule('TEST_H2',ttype_op,DEF_OP_FROM_OCC,
     &              'TEST_H2',1,1,
     &              parameters,2,tgt_info)

      occ_def = 0
      ndef = 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IEXTR,1,1) = 1
      occ_def(IPART,2,1) = 2
      call add_target('TEST_H',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/.true.,.true./),ndef)
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
     &     occ_def,ndef,1,(/.true.,.true./),ndef)
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
     &     occ_def,ndef,1,(/.true.,.true./),ndef)
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

*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*
      ! test I
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FOPT_TEST_I'
      labels(2) = 'FORM_TEST_I'
      ncat = 1
      nint = 0
      call add_target('FOPT_TEST_I',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_I','FORM_TEST_I',tgt_info)
      call set_dependency('FOPT_TEST_I','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_I','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_I','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_I','ME_TEST_A',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_I',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

      ! test I
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'FOPT_TEST_II'
      labels(2) = 'FORM_TEST_II'
      ncat = 1
      nint = 0
      call add_target('FOPT_TEST_II',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_TEST_II','FORM_TEST_II',tgt_info)
      call set_dependency('FOPT_TEST_II','DEF-ME_TEST_RES',tgt_info)
      call set_dependency('FOPT_TEST_II','ME_TEST_H',tgt_info)
      call set_dependency('FOPT_TEST_II','ME_TEST_R',tgt_info)
      call set_dependency('FOPT_TEST_II','ME_TEST_A',tgt_info)
      call set_dependency('FOPT_TEST_II','ME_TEST_B',tgt_info)
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('FOPT_TEST_II',ttype_frm,OPTIMIZE,
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
     &     (/1,1,1,1/),
     &     (/1,1,1,1/),
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
     &     (/1,1,1,1/),
     &     (/1,1,1,1/),
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

*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*
      needed = testnr.eq.1
      
      call add_target('SIGN_TEST_I',ttype_gen,needed,tgt_info)
      call set_dependency('SIGN_TEST_I','FOPT_TEST_I',tgt_info)
      labels(1) = 'FOPT_TEST_I'
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,3,5,2,2,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,3,6,5,2,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,1,1,1,2,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,7,5,1,1,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,7,6,3,2,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_I'
      call modify_parameters(-1,parameters,7,(/1,5,5,6,4,3,1/),7)
      call set_rule('SIGN_TEST_I',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_I',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_I',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      ! test 2
      needed = testnr.eq.2
      
      call add_target('SIGN_TEST_II',ttype_gen,needed,tgt_info)
      call set_dependency('SIGN_TEST_II','FOPT_TEST_II',tgt_info)
      labels(1) = 'FOPT_TEST_II'
      call set_rule('SIGN_TEST_II',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_II',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      
      labels(1) = 'FOPT_TEST_II'
      call modify_parameters(-1,parameters,7,(/1,5,3,4,6,2,1/),7)
      call set_rule('SIGN_TEST_II',ttype_frm,MODIFY_FACTORIZATION,
     &     labels,1,0,
     &     parameters,1,tgt_info)
      call set_rule('SIGN_TEST_II',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)
      labels(1) = 'ME_TEST_RES'
      call set_rule('SIGN_TEST_II',ttype_opme,RES_ME_LIST,
     &     labels,1,0,
     &     parameters,0,tgt_info)

      return

      end
