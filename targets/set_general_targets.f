*----------------------------------------------------------------------*
      subroutine set_general_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set targets needed in more or less all kinds of CC calculations
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
      character(*), intent(in) ::
     &     env_type

      integer ::
     &     min_rank, max_rank,
     &     isim, ncat, nint, icnt, iformal, extern,
     &     isym, ms, msc, sym_arr(8),trunc_type,t1ext_mode,
     &     maxr12exc, icase, icaseF, norb
      logical ::
     &     needed, explicit, truncate
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(10)
      character(len_command_par) ::
     &     parameters, z2str
      integer, allocatable ::
     &     iRsys(:)

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting general targets ...'

      call get_argument_value('method.R12','trunc',ival=trunc_type)
      truncate = trunc_type.ge.0
      if (is_keyword_set('method.truncate').gt.0) then
        truncate = is_keyword_set('method.truncate').gt.0
        call get_argument_value('method.truncate','trunc_type',
     &       ival=trunc_type)
      end if
      call get_argument_value('method.R12','maxexc',ival=maxr12exc)
      call get_argument_value('method.R12','Z2_appr',str=z2str)

      call get_argument_value('method.CCPT','extern',ival=extern)
      call get_argument_value('method.CC','T1ext',ival=t1ext_mode)
      if (t1ext_mode.eq.0)
     &     call get_argument_value('method.R12','T1ext',ival=t1ext_mode)
      if (t1ext_mode.eq.0)
     &     call get_argument_value('method.ECC','T1ext',ival=t1ext_mode)

      msc = +1  ! assuming closed shell
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! Unit operator (purely formal)
      ! we use DEF_HAMILTONIAN with rank 0 to 0 and iformal=0
      call add_target2(op_unity,.false.,tgt_info)
      call set_rule2(op_unity,DEF_HAMILTONIAN,tgt_info)
      call set_arg(op_unity,DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/op_unity/))
      call set_arg(op_unity,DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/0/))
      call set_arg(op_unity,DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/0/))

      ! Hamiltonian
      iformal = 1
      explicit = is_keyword_set('method.R12').gt.0
      if (explicit.and.orb_info%caborb.gt.0.and.(.not.truncate
     &     .or.(truncate.and.trunc_type.eq.1)).or.extern.ge.2)
     &     iformal = 4
      if (explicit.and.orb_info%caborb.gt.0.and.truncate
     &     .and.trunc_type.ne.1.or.extern.eq.1)
     &     iformal = 3
      call add_target2(op_ham,.false.,tgt_info)
c patch for CCPT-R12 tests:
c      if (explicit.and.
c     &     (is_keyword_set('method.CCPT').gt.0.or.maxr12exc.ge.3))
      if (explicit.and.
     &     trim(z2str).eq.'direct'.or.trim(z2str).eq.'DIRECT')
     &     iformal = 4
c patch end
      if (t1ext_mode.gt.0) iformal = min(5,max(t1ext_mode+2,iformal))
c another patch
      if (is_keyword_set('method.ECC').gt.0)
     &     iformal = min(6,max(t1ext_mode+1,4))
c patch end
      call set_rule2(op_ham,DEF_HAMILTONIAN,tgt_info)
      call set_arg(op_ham,DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/op_ham/))
      call set_arg(op_ham,DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/0/))
      call set_arg(op_ham,DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/))
      call set_arg(op_ham,DEF_HAMILTONIAN,'FORMAL',1,tgt_info,
     &     val_int=(/iformal/))
      call set_arg(op_ham,DEF_HAMILTONIAN,'SET_X',1,tgt_info,
     &     val_log=(/explicit.or.extern.gt.0.or.t1ext_mode.gt.0/))
      call set_rule2(op_ham,SET_HERMIT,tgt_info)
      call set_arg(op_ham,SET_HERMIT,'LABEL',1,tgt_info,
     &     val_label=(/op_ham/))
      call set_arg(op_ham,SET_HERMIT,'CA_SYMMETRY',1,tgt_info,
     &     val_int=(/+1/))

      ! Fock-Operator
      call add_target2(op_fock,.false.,tgt_info)
      call set_rule2(op_fock,DEF_HAMILTONIAN,tgt_info)
      call set_arg(op_fock,DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/op_fock/))
      call set_arg(op_fock,DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_fock,DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(op_fock,DEF_HAMILTONIAN,'FORMAL',1,tgt_info,
     &     val_int=(/iformal/))
      call set_arg(op_fock,DEF_HAMILTONIAN,'SET_X',1,tgt_info,
     &     val_log=(/explicit.or.extern.gt.0.or.t1ext_mode.gt.0/))
      call set_rule2(op_fock,SET_HERMIT,tgt_info)
      call set_arg(op_fock,SET_HERMIT,'LABEL',1,tgt_info,
     &     val_label=(/op_fock/))
      call set_arg(op_fock,SET_HERMIT,'CA_SYMMETRY',1,tgt_info,
     &     val_int=(/+1/))

*----------------------------------------------------------------------*
*     Formulae
*----------------------------------------------------------------------*
      if (is_keyword_set('calculate.check.formulae').gt.0) then
        call add_target2(form_test,.true.,tgt_info)
        call set_rule2(form_test,CHECK_FORMGEN,tgt_info)
      end if
      
*----------------------------------------------------------------------*
*     Opt. Formulae
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! Hamilton list:
      call add_target2(mel_ham,.false.,tgt_info)
      call set_dependency(mel_ham,op_ham,tgt_info)
      ! (a) define
      call set_rule2(mel_ham,DEF_ME_LIST,tgt_info)
      call set_arg(mel_ham,DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/mel_ham/))
      call set_arg(mel_ham,DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/op_ham/))
      call set_arg(mel_ham,DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &     val_int=(/msc/))
      call set_arg(mel_ham,DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &     val_int=(/0/))
      call set_arg(mel_ham,DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg(mel_ham,DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      ! (b) import
      call set_rule2(mel_ham,IMPORT,tgt_info)
      call set_arg(mel_ham,IMPORT,'LIST',1,tgt_info,
     &     val_label=(/mel_ham/))
      call set_arg(mel_ham,IMPORT,'TYPE',1,tgt_info,
     &     val_str='H_INT')
      ! (c) modify, if requested
      if (is_keyword_set('orb_space.GEtest').gt.0) then
        call get_argument_dimension(norb,'orb_space.GEtest','Rsys')
        allocate(iRsys(max(1,norb)))
        call get_argument_value('orb_space.GEtest','Rsys',iarr=iRsys)
        call get_argument_value('orb_space.GEtest','case',ival=icase)
        call get_argument_value('orb_space.GEtest','caseF',ival=icaseF)
        call set_rule2(mel_ham,GETEST,tgt_info)
        call set_arg(mel_ham,GETEST,'LIST',1,tgt_info,
     &       val_label=(/mel_ham/))
        call set_arg(mel_ham,GETEST,'R-SYS',norb,tgt_info,
     &       val_int=iRsys) 
        call set_arg(mel_ham,GETEST,'CASE',1,tgt_info,
     &       val_int=(/icase/))
        call set_arg(mel_ham,GETEST,'SPLIT-FOCK',1,tgt_info,
     &       val_int=(/icaseF/)) 
        deallocate(iRsys)
      end if
*----------------------------------------------------------------------*
*     "phony" targets
*----------------------------------------------------------------------*

      return
      end
