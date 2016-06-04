*----------------------------------------------------------------------*
      subroutine set_unc_mrci_targets(tgt_info,orb_info,calc,
     &                                name_infile,name_orbinfo)
*----------------------------------------------------------------------*
*     set targets for CASSCF or uncontracted CI coefficient optimization
*
*     matthias, fall 2009
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'ifc_targets.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'
      include 'routes.h'

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      logical, intent(in) ::
     &     calc
      character(*), intent(in) ::
     &     name_infile, name_orbinfo

      integer ::
     &     ndef, occ_def(ngastp,2,124),!60),
     &     isym, msc, ims, ip, ih, cminexc, 
     &     cminh, cmaxh, cminp, cmaxp, cmaxexc, ciroot, maxroot, cmaxv,
     &     guess, refproj, spinexpec, n_states, i_state, optref,
     &     ncnt, ncnt2, icnt, arglen, i, nnn
      logical ::
     &     oldref, l_exist, writeF, multistate
      character(len_target_name) ::
     &     dia_label, dia_label2, labels(20), c_st
      character(100) ::
     &     label_title
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     opstr*4, mestr*7, filestr*15
      character(len_target_name), external ::
     &     state_label
      real(8), allocatable ::
     &     xscr(:)
      character(len=256) ::
     &     gecco_path

      if (iprlvl.gt.0)
     &     write(lulog,*) 'setting targets for multiref. wave function'

      ims = orb_info%ims
      if (ims.eq.0.and.mod(orb_info%imult-1,4).eq.0) then
        msc = 1
      else if (ims.eq.0.and.mod(orb_info%imult+1,4).eq.0) then
        msc = -1
      else
        msc = 0
      end if

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('method.MR','cminh',
     &     ival=cminh)
      call get_argument_value('method.MR','cmaxh',
     &     ival=cmaxh)
      call get_argument_value('method.MR','cminp',
     &     ival=cminp)
      call get_argument_value('method.MR','cmaxp',
     &     ival=cmaxp)
      call get_argument_value('method.MR','cmaxexc',
     &     ival=cmaxexc)
      call get_argument_value('method.MR','cminexc',
     &     ival=cminexc)
      call get_argument_value('method.MR','ciroot',
     &     ival=ciroot)
      call get_argument_value('method.MR','maxroot',
     &     ival=maxroot)
      if(maxroot.le.0) maxroot=ciroot
      call get_argument_value('method.MR','multistate',
     &     lval=multistate)
      call get_argument_value('method.MR','oldref',
     &     lval=oldref)
      call get_argument_value('method.MR','writeFock',
     &     lval=writeF)
      call get_argument_value('method.MR','guess',
     &     ival=guess)
      call get_argument_value('method.MR','refproj',
     &     ival=refproj)
      call get_argument_value('method.MR','spinexpec',
     &     ival=spinexpec)
      if (cmaxh.lt.0) cmaxh = cmaxexc
      if (cmaxp.lt.0) cmaxp = cmaxexc
      cmaxv = orb_info%norb_hpv(IVALE,1)*2
      if(multistate)then
       n_states = ciroot
      else
       n_states = 1
      end if

      call get_environment_variable( "GECCO_DIR", value=gecco_path)

      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)

      if (ntest.ge.100) then
        write(lulog,*) 'cminh      = ',cminh
        write(lulog,*) 'cmaxh      = ',cmaxh
        write(lulog,*) 'cminp      = ',cminp
        write(lulog,*) 'cmaxp      = ',cmaxp
        write(lulog,*) 'cmaxv      = ',cmaxv
        if (cminexc.gt.0) write(lulog,*) 'cminexc    = ',cminexc
        write(lulog,*) 'cmaxexc    = ',cmaxexc
        write(lulog,*) 'nactel     = ',orb_info%nactel
        write(lulog,*) 'ciroot     = ',ciroot
        write(lulog,*) 'maxroot    = ',maxroot
        write(lulog,*) 'multistate = ',multistate
        write(lulog,*) 'oldref     = ',oldref
        write(lulog,*) 'spinadapt  = ',spinadapt
        if (guess.gt.0) write(lulog,*) 'guess    =',guess
        if (refproj.gt.0) write(lulog,*) 'refproj  =',refproj
      end if

      if (guess.gt.0) then
        inquire(file='ME_C0start_list.da',exist=l_exist)
        if (.not.l_exist) call quit(1,'set_unc_mrci_targets',
     &           'Initial guess: Did not find file ME_C0start_list.da!')
        if (maxroot.gt.1) call quit(1,'set_unc_mrci_targets',
     &           'option guess only available for ciroot=1/maxroot=1')
      end if
      if (refproj.gt.0) then
!       if (refproj.gt.9.or.spinadapt.gt.0)
!    &     call quit(1,'set_unc_mrci_targets',
!    &          'refproj>9 or with spinadapt>0 not available yet')
        filestr='ME_C0_x_list.da'
        do ip = 1, refproj
          write(filestr(7:7),'(i1)') ip
          inquire(file=filestr,exist=l_exist)
          if (.not.l_exist) call quit(1,'set_unc_mrci_targets',
     &         'Reference projection: Did not find file '//filestr//'!')
        end do
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
          if (orb_info%nactel+ih-ip.lt.0.or.
     &        orb_info%nactel+ih-ip.gt.cmaxv) cycle
          if (max(ih,ip).lt.cminexc) cycle
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
      if (multistate)
     &     call set_multistate_operator(tgt_info,n_states,'C0',1)

      ! clone of C0 for spin projections
      call add_target2('C0_sp',.false.,tgt_info)
      call set_dependency('C0_sp','C0',tgt_info)
      call set_rule2('C0_sp',CLONE_OP,tgt_info)
      call set_arg('C0_sp',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'C0_sp'/))
      call set_arg('C0_sp',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'C0'/))

      ! clone: keep copy of unmodified coefficients C00
      call add_target3((/
     &     'target C00(                             ',
     &     '  CLONE_OPERATOR(label=C00,template=C0))'/),tgt_info)
 
      ! define product Jacobian times C0
      call add_target('A_C0',ttype_op,.false.,tgt_info)
      call set_dependency('A_C0','C0',tgt_info)
      call cloneop_parameters(-1,parameters,'C0',.false.)
      call set_rule('A_C0',ttype_op,CLONE_OP,'A_C0',1,1,
     &              parameters,1,tgt_info)
      if (multistate)
     &     call set_multistate_operator(tgt_info,n_states,'A_C0',2)

      ! Diagonal Preconditioner for reference
      call add_target(trim(op_dia)//'_C0',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(trim(op_dia)//'_C0','C0',tgt_info)
      call cloneop_parameters(-1,parameters,'C0',.false.)
      call set_rule(trim(op_dia)//'_C0',ttype_op,CLONE_OP,
     &              trim(op_dia)//'_C0',1,1,
     &              parameters,1,tgt_info)

      ! define Fock operator wrt reference function
      ! for all current methods, we only need the purely inactive part
      ! (valence part might actually mess up something if prc_type=3)
c      call add_target('FREF',ttype_op,.false.,tgt_info)
c      call hop_parameters(-1,parameters,0,1,1,.false.)
c      call set_rule('FREF',ttype_op,DEF_HAMILTONIAN,'FREF',
c     &              1,1,parameters,1,tgt_info)
      call add_target2('FREF',.false.,tgt_info)
      call set_rule2('FREF',DEF_HAMILTONIAN,tgt_info)
      call set_arg('FREF',DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/'FREF'/))
      call set_arg('FREF',DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/1/))
      if (writeF) then
        if (.not.calc) call quit(1,'set_unc_mrci_targets',
     &    'Please use writeFock=T only in pure CASSCF calculations')
        call set_arg('FREF',DEF_HAMILTONIAN,'X_SPCS',1,tgt_info,
     &       val_int=(/IEXTR/))
      else
        call set_arg('FREF',DEF_HAMILTONIAN,'X_SPCS',2,tgt_info,
     &       val_int=(/IVALE,IEXTR/))
      end if

c add special FREF0
      call add_target2('FREF0',.false.,tgt_info)
      call set_rule2('FREF0',DEF_HAMILTONIAN,tgt_info)
      call set_arg('FREF0',DEF_HAMILTONIAN,'LABEL',1,tgt_info,
     &     val_label=(/'FREF0'/))
      call set_arg('FREF0',DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('FREF0',DEF_HAMILTONIAN,'X_SPCS',2,tgt_info,
     &       val_int=(/IVALE,IEXTR/))


      ! Spin operators
      ! S+
      call add_target('S+',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 3
      occ_def(IHOLE,1,1) = 1
      occ_def(IHOLE,2,1) = 1
      occ_def(IPART,1,2) = 1
      occ_def(IPART,2,2) = 1
      occ_def(IVALE,1,3) = 1
      occ_def(IVALE,2,3) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('S+',ttype_op,DEF_OP_FROM_OCC,
     &              'S+',1,1,
     &              parameters,2,tgt_info)
      ! S-
      call add_target('S-',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('S-','S+',tgt_info)
      call cloneop_parameters(-1,parameters,'S+',.false.)
      call set_rule('S-',ttype_op,CLONE_OP,'S-',1,1,
     &              parameters,1,tgt_info)
      ! Sz
      call add_target('Sz',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('Sz','S+',tgt_info)
      call cloneop_parameters(-1,parameters,'S+',.false.)
      call set_rule('Sz',ttype_op,CLONE_OP,'Sz',1,1,
     &              parameters,1,tgt_info)

      ! define scalar spin expectation value
      call add_target('S(S+1)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('S(S+1)',ttype_op,DEF_HAMILTONIAN,'S(S+1)',
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
c dbg
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_A_C0',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)
c      call set_rule2('F_A_C0',ABORT,tgt_info)
c dbg end

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
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_FREF',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! The same but for C00
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_FREF0'
      labels(2) = 'FREF0'
      labels(3) = 'FREF0'
      labels(4) = 'C00^+'
      labels(5) = op_ham
      labels(6) = 'C00'
      labels(7) = 'FREF0'
      call add_target('F_FREF0',ttype_frm,.false.,tgt_info)
      call set_dependency('F_FREF0','FREF0',tgt_info)
      call set_dependency('F_FREF0',op_ham,tgt_info)
      call set_dependency('F_FREF0','C00',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'reference energy expression',5,
     &     (/1,2,3,4,1/),
     &     (/-1,-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1,-1/),
     &     0,0,
     &     (/1,4,2,5/),2,
     &     0,0)
      call set_rule('F_FREF0',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,7,1,
     &              parameters,3,tgt_info)

      ! expectation value of S^2
      call add_target2('F_REF_S(S+1)',.false.,tgt_info)
      call set_dependency('F_REF_S(S+1)','S(S+1)',tgt_info)
      call set_dependency('F_REF_S(S+1)','C0',tgt_info)
      call set_dependency('F_REF_S(S+1)','S+',tgt_info)
      call set_dependency('F_REF_S(S+1)','S-',tgt_info)
      ! (a) 1/2*(S+S- + S-S+) = {S+S-} + 1/2*(S+S-)_c + 1/2*(S-S+)_c
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','S+  ','S-  ','C0  '/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','S+  ','S-  ','C0  '/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'CONNECT',2,
     &     tgt_info,val_int=(/2,3/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','S-  ','S+  ','C0  '/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'CONNECT',2,
     &     tgt_info,val_int=(/2,3/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      if (ims.ne.0) then
        ! (b) + Sz^2
        call set_dependency('F_REF_S(S+1)','Sz',tgt_info)
        call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,val_label=(/'F_REF_S(S+1)'/))
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,val_label=(/'S(S+1)'/))
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &       tgt_info,
     &       val_label=(/'C0^+  ','Sz    ','Sz    ','C0    '/))
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,
     &       tgt_info,val_int=(/2,3,4,5/))
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'FIX_VTX',1,
     &       tgt_info,val_log=(/.true./)) !prevent "BCH factor"
      end if
c dbg
c      call set_rule2('F_REF_S(S+1)',PRINT_FORMULA,tgt_info)
c      call set_arg('F_REF_S(S+1)',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_REF_S(S+1)'/))
c dbgend

      ! expectation value of S^2 in operator form
      call add_target2('F_REPL_H_S2',.false.,tgt_info)
      call set_dependency('F_REPL_H_S2','S+',tgt_info)
      call set_dependency('F_REPL_H_S2','S-',tgt_info)
      call set_dependency('F_REPL_H_S2','Sz',tgt_info)
      call set_dependency('F_REPL_H_S2','H',tgt_info)
      ! (a) 1/2*(S+S- + S-S+)
      call set_rule2('F_REPL_H_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'H ','S+','S-','H '/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_rule2('F_REPL_H_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'H ','S-','S+','H '/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! (b) + Sz^2
      call set_rule2('F_REPL_H_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'H     ','Sz    ','Sz    ','H     '/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'FIX_VTX',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_REPL_H_S2',PRINT_FORMULA,tgt_info)
c      call set_arg('F_REPL_H_S2',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_REPL_H_S2'/))
c dbgend

      ! formula for spin projection: C0_sp = S^2 C0
      call add_target2('F_C0_sp',.false.,tgt_info)
      call set_dependency('F_C0_sp','C0_sp',tgt_info)
      call set_dependency('F_C0_sp','C0',tgt_info)
      call set_dependency('F_C0_sp','S+',tgt_info)
      call set_dependency('F_C0_sp','S-',tgt_info)
      ! (a) 1/2*(S+S- + S-S+) = {S+S-} + 1/2*(S+S-)_c + 1/2*(S-S+)_c
      call set_rule2('F_C0_sp',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'C0_sp','S+   ','S-   ','C0   ','C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,4,1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/2,3/))
      call set_rule2('F_C0_sp',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'C0_sp','S+   ','S-   ','C0   ','C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,4,1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'CONNECT',2,
     &     tgt_info,val_int=(/2,3/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_C0_sp',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'C0_sp','S-   ','S+   ','C0   ','C0_sp'/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/1,2,3,4,1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'N_CONNECT',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'CONNECT',2,
     &     tgt_info,val_int=(/2,3/))
      call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      if (ims.ne.0) then
        ! (b) + Sz^2
        call set_dependency('F_C0_sp','Sz',tgt_info)
        call set_rule2('F_C0_sp',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_C0_sp'/))
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'C0_sp'/))
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'C0_sp ','Sz    ','Sz    ','C0    ','C0_sp '/))
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,4,1/))
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        call set_arg('F_C0_sp',EXPAND_OP_PRODUCT,'FIX_VTX',1,tgt_info,
     &       val_log=(/.true./))
      end if
c dbg
c      call set_rule2('F_C0_sp',PRINT_FORMULA,tgt_info)
c      call set_arg('F_C0_sp',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_C0_sp'/))
c dbgend

      ! formula for reference projection: -|C1><C0|C1>-|C2><C0|C2>...
      call add_target2('F_C0_prj',.false.,tgt_info)
      call set_dependency('F_C0_prj','C0',tgt_info)
      call set_dependency('F_C0_prj','DEF_ME_C0_x',tgt_info)
      opstr='C0_x'
      do ip = 1, refproj
        write(opstr(4:4),'(i1)') ip
        call set_rule2('F_C0_prj',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_C0_prj'/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'C0'/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'C0  ',opstr,'C0^+',opstr,'C0  '/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,4,1/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &       tgt_info,val_int=(/2/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'AVOID',4,tgt_info,
     &       val_int=(/1,4,3,5/))
        call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &       val_rl8=(/-1d0/))
        if (ip.gt.1)
     &     call set_arg('F_C0_prj',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &          val_log=(/.false./))
      end do
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

      ! Fock operator wrt zeroth-order reference function
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_FREF0'
      labels(2) = 'F_FREF0'
      call add_target('FOPT_FREF0',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_FREF0','F_FREF0',tgt_info)
      call set_dependency('FOPT_FREF0',mel_ham,tgt_info)
      call set_dependency('FOPT_FREF0','DEF_ME_C00',tgt_info)
      call set_dependency('FOPT_FREF0','DEF_ME_FREF0',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_FREF0',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! spin expectation value expression
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_REF_S(S+1)'
      labels(2) = 'F_REF_S(S+1)'
      call add_target('FOPT_REF_S(S+1)',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_REF_S(S+1)','F_REF_S(S+1)',tgt_info)
      call set_dependency('FOPT_REF_S(S+1)','DEF_ME_S(S+1)',tgt_info)
      call set_dependency('FOPT_REF_S(S+1)','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_REF_S(S+1)','DEF_ME_S+',tgt_info)
      call set_dependency('FOPT_REF_S(S+1)','DEF_ME_S-',tgt_info)
      if (ims.ne.0) 
     &   call set_dependency('FOPT_REF_S(S+1)','DEF_ME_Sz',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_REF_S(S+1)',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! formula for spin projection
      call add_target2('FOPT_C0_sp',.false.,tgt_info)
      call set_dependency('FOPT_C0_sp','F_C0_sp',tgt_info)
      call set_dependency('FOPT_C0_sp','DEF_ME_C0_sp',tgt_info)
      call set_dependency('FOPT_C0_sp','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_C0_sp','DEF_ME_S+',tgt_info)
      call set_dependency('FOPT_C0_sp','DEF_ME_S-',tgt_info)
      if (ims.ne.0)
     &   call set_dependency('FOPT_C0_sp','DEF_ME_Sz',tgt_info)
      call set_rule2('FOPT_C0_sp',OPTIMIZE,tgt_info)
      call set_arg('FOPT_C0_sp',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_C0_sp'/))
      call set_arg('FOPT_C0_sp',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_C0_sp'/))

      ! formula for reference projection
      call add_target2('FOPT_C0_prj',.false.,tgt_info)
      call set_dependency('FOPT_C0_prj','F_C0_prj',tgt_info)
      call set_dependency('FOPT_C0_prj','DEF_ME_C0',
     &     tgt_info)
      call set_dependency('FOPT_C0_prj','DEF_ME_C0_x',tgt_info)
      call set_rule2('FOPT_C0_prj',OPTIMIZE,tgt_info)
      call set_arg('FOPT_C0_prj',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_C0_prj'/))
      call set_arg('FOPT_C0_prj',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_C0_prj'/))
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME_C0
      call add_target2('DEF_ME_C0',.false.,tgt_info)
      call set_dependency('DEF_ME_C0','C0',tgt_info)
      call set_rule2('DEF_ME_C0',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C0'/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C0'/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/maxroot/))
      if (multistate) then
       call set_arg('DEF_ME_C0',DEF_ME_LIST,'REC',1,tgt_info,
     &      val_int=(/1/))
      else
       call set_arg('DEF_ME_C0',DEF_ME_LIST,'REC',1,tgt_info,
     &      val_int=(/ciroot/))
      end if
!     ME_C0_<i_state>
      if (multistate) then
       do i_state=1,n_states,1
        c_st = state_label(i_state,.true.)
        call set_rule2('DEF_ME_C0',DEF_ME_LIST,tgt_info)
        call set_arg('DEF_ME_C0',DEF_ME_LIST,'LIST',1,tgt_info,
     &       val_label=(/'ME_C0'//trim(c_st)/))
        call set_arg('DEF_ME_C0',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &       val_label=(/'C0'//trim(c_st)/))
        call set_arg('DEF_ME_C0',DEF_ME_LIST,'2MS',1,tgt_info,
     &       val_int=(/ims/))
        call set_arg('DEF_ME_C0',DEF_ME_LIST,'IRREP',1,tgt_info,
     &       val_int=(/orb_info%lsym/))
        call set_arg('DEF_ME_C0',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &       val_int=(/msc/))
       end do
      end if

      ! ME_C0_sp
      call add_target2('DEF_ME_C0_sp',.false.,tgt_info)
      call set_dependency('DEF_ME_C0_sp','C0_sp',tgt_info)
      call set_rule2('DEF_ME_C0_sp',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_C0_sp',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C0_sp'/))
      call set_arg('DEF_ME_C0_sp',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C0_sp'/))
      call set_arg('DEF_ME_C0_sp',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('DEF_ME_C0_sp',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('DEF_ME_C0_sp',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME_C00
      call add_target2('DEF_ME_C00',.false.,tgt_info)
      call set_dependency('DEF_ME_C00','C00',tgt_info)
      call set_dependency('DEF_ME_C00','SOLVE_REF',tgt_info)
      call set_rule2('DEF_ME_C00',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_C00',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C00'/))
      call set_arg('DEF_ME_C00',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C00'/))
      call set_arg('DEF_ME_C00',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('DEF_ME_C00',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('DEF_ME_C00',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      if(multistate) then
       call set_arg('DEF_ME_C00',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &              val_int=(/1/))
       call set_arg('DEF_ME_C00',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &              val_int=(/n_states/))
       call set_arg('DEF_ME_C00',DEF_ME_LIST,'REC',1,tgt_info,
     &              val_int=(/1/))
      end if
      nnn = n_states
      if (oldref) nnn = 0  ! no copying for oldref==T
      do i_state = 1, nnn
        c_st = state_label(i_state,.false.)
        call set_rule2('DEF_ME_C00',SCALE_COPY,tgt_info)
        call set_arg('DEF_ME_C00',SCALE_COPY,'LIST_RES',1,tgt_info,
     &       val_label=(/'ME_C00'/))
        call set_arg('DEF_ME_C00',SCALE_COPY,'LIST_INP',1,tgt_info,
     &       val_label=(/'ME_C0'/))
        call set_arg('DEF_ME_C00',SCALE_COPY,'FAC',1,tgt_info,
     &       val_rl8=(/1d0/))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Saved C0, C00: ',0,'LIST')
c       call set_rule('DEF_ME_C00',ttype_opme,PRINT_MEL,
c     &     'ME_C00',1,0,
c     &     parameters,2,tgt_info)
c dbg end
      if(multistate)then
       call set_rule2('DEF_ME_C00',ADV_STATE,tgt_info)
       call set_arg('DEF_ME_C00',ADV_STATE,'LISTS',2,tgt_info,
     &      val_label=['ME_C0 ',
     &                 'ME_C00'])
       call set_arg('DEF_ME_C00',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
      end if
      end do

      ! ME_C0_x
      call add_target2('DEF_ME_C0_x',.false.,tgt_info)
      call set_dependency('DEF_ME_C0_x','C0',tgt_info)
      opstr='C0_x'
      mestr='ME_C0_x'
      do ip = 1, refproj
        write(opstr(4:4),'(i1)') ip
        write(mestr(7:7),'(i1)') ip
        call set_rule2('DEF_ME_C0_x',CLONE_OP,tgt_info)
        call set_arg('DEF_ME_C0_x',CLONE_OP,'LABEL',1,tgt_info,
     &       val_label=(/opstr/))
        call set_arg('DEF_ME_C0_x',CLONE_OP,'TEMPLATE',1,tgt_info,
     &       val_label=(/'C0'/))
        call set_rule2('DEF_ME_C0_x',DEF_ME_LIST,tgt_info)
        call set_arg('DEF_ME_C0_x',DEF_ME_LIST,'LIST',1,tgt_info,
     &               val_label=(/mestr/))
        call set_arg('DEF_ME_C0_x',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &               val_label=(/opstr/))
        call set_arg('DEF_ME_C0_x',DEF_ME_LIST,'2MS',1,tgt_info,
     &               val_int=(/ims/))
        call set_arg('DEF_ME_C0_x',DEF_ME_LIST,'IRREP',1,tgt_info,
     &               val_int=(/orb_info%lsym/))
        call set_arg('DEF_ME_C0_x',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &               val_int=(/msc/))
      end do

      ! ME_A_C0
      call add_target('DEF_ME_A_C0',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_A_C0','A_C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_A_C0'
      labels(2) = 'A_C0'
      if(optref.eq.-1.or.optref.eq.-2)then
       do i_state=2,n_states
        c_st = state_label(i_state,.false.)
        labels(1) = 'ME_A_C0'//trim(c_st)
        call me_list_parameters(-1,parameters,
     &       msc,0,orb_info%lsym,
     &       0,ims,.false.)
        call set_rule('DEF_ME_A_C0',ttype_opme,DEF_ME_LIST,
     &                labels,2,1,
     &                 parameters,1,tgt_info)
       end do
      end if
      labels(1) = 'ME_A_C0'
      call me_list_parameters(-1,parameters,
     &     msc,0,orb_info%lsym,
     &     0,ims,.false.)
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
     &                    trim(op_dia)//'_'//'C0',tgt_info)
      if(optref.eq.-1.or.optref.eq.-2)then
      do i_state = 2,n_states
      c_st = state_label(i_state,.false.)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C0'//trim(c_st)
      labels(2) = trim(op_dia)//'_'//'C0'
      call me_list_parameters(-1,parameters,
     &     0,0,orb_info%lsym,
     &     0,ims,.false.)
      call set_rule(trim(dia_label)//'C0',ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      labels(1) = trim(dia_label)//'C0'//trim(c_st)
      labels(2) = mel_ham
      call set_rule(trim(dia_label)//'C0',ttype_opme,
     &              PRECONDITIONER,
     &              labels,2,1,
     &              parameters,2,tgt_info)
      end do
      end if
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C0'
      labels(2) = trim(op_dia)//'_'//'C0'
      call me_list_parameters(-1,parameters,
     &     0,0,orb_info%lsym,
     &     0,ims,.false.)
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
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner :',0,'LIST')
c      call set_rule(trim(dia_label)//'C0',ttype_opme,PRINT_MEL,
c     &     trim(dia_label)//'C0',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME_FREF
      call add_target2('DEF_ME_FREF',.false.,tgt_info)
      call set_dependency('DEF_ME_FREF','FREF',tgt_info)
      call set_rule2('DEF_ME_FREF',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_FREF'/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'FREF'/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      if (ims.eq.0)
     &   call set_arg('DEF_ME_FREF',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &        val_int=(/1/))
      if (spinadapt.ge.2)
     &   call set_arg('DEF_ME_FREF',DEF_ME_LIST,'S2',1,tgt_info,
     &        val_int=(/0/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_FREF',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))

      ! ME_FREF0
      call add_target2('DEF_ME_FREF0',.false.,tgt_info)
      call set_dependency('DEF_ME_FREF0','FREF0',tgt_info)
      call set_rule2('DEF_ME_FREF0',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_FREF0'/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'FREF0'/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      if (ims.eq.0)
     &   call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &        val_int=(/1/))
      if (spinadapt.ge.2)
     &   call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'S2',1,tgt_info,
     &        val_int=(/0/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_FREF0',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))


      ! ME_S+
      call add_target2('DEF_ME_S+',.false.,tgt_info)
      call set_dependency('DEF_ME_S+','S+',tgt_info)
      call set_rule2('DEF_ME_S+',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_S+',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_S+'/))
      call set_arg('DEF_ME_S+',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'S+'/))
      call set_arg('DEF_ME_S+',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_S+',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/2/))
      call set_rule2('DEF_ME_S+',UNITY,tgt_info)
      call set_arg('DEF_ME_S+',UNITY,'LIST',1,tgt_info,
     &             val_label=(/'ME_S+'/))
      call set_arg('DEF_ME_S+',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('DEF_ME_S+',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.true./))
      ! ME_S-
      call add_target2('DEF_ME_S-',.false.,tgt_info)
      call set_dependency('DEF_ME_S-','S-',tgt_info)
      call set_rule2('DEF_ME_S-',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_S-',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_S-'/))
      call set_arg('DEF_ME_S-',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'S-'/))
      call set_arg('DEF_ME_S-',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_S-',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/-2/))
      call set_rule2('DEF_ME_S-',UNITY,tgt_info)
      call set_arg('DEF_ME_S-',UNITY,'LIST',1,tgt_info,
     &             val_label=(/'ME_S-'/))
      call set_arg('DEF_ME_S-',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('DEF_ME_S-',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.true./))
      ! ME_Sz
      call add_target2('DEF_ME_Sz',.false.,tgt_info)
      call set_dependency('DEF_ME_Sz','Sz',tgt_info)
      call set_rule2('DEF_ME_Sz',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Sz',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_Sz'/))
      call set_arg('DEF_ME_Sz',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'Sz'/))
      call set_arg('DEF_ME_Sz',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_Sz',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_rule2('DEF_ME_Sz',UNITY,tgt_info)
      call set_arg('DEF_ME_Sz',UNITY,'LIST',1,tgt_info,
     &             val_label=(/'ME_Sz'/))
      call set_arg('DEF_ME_Sz',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Sz',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.true./))

      ! ME_S(S+1)
      call add_target2('DEF_ME_S(S+1)',.false.,tgt_info)
      call set_dependency('DEF_ME_S(S+1)','S(S+1)',tgt_info)
      call set_rule2('DEF_ME_S(S+1)',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_S(S+1)'/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'S(S+1)'/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'S2',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &             val_log=(/.false./))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_S(S+1)',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! SOLVE reference eigenvalue equation
      call add_target2('SOLVE_REF',.false.,tgt_info)
      call set_dependency('SOLVE_REF','FOPT_A_C0',tgt_info)
      call me_list_label(dia_label,mel_dia,orb_info%lsym,
     &     0,0,0,.false.)
      dia_label2 = trim(dia_label)//'C0'
      call set_dependency('SOLVE_REF',trim(dia_label2),tgt_info)
      if (spinadapt.gt.0) then
        call set_dependency('SOLVE_REF','DEF_ME_C0_sp',tgt_info)
        call set_dependency('SOLVE_REF','FOPT_C0_sp',tgt_info)
      end if
      if (refproj.gt.0)
     &   call set_dependency('SOLVE_REF','FOPT_C0_prj',tgt_info)
      if (guess.gt.0)
     &   call set_dependency('SOLVE_REF','C0guess',tgt_info)
      if (multistate) then
       do i_state = 1,n_states
        c_st = state_label(i_state,.true.)
        labels(i_state) = "ME_C0"//trim(c_st)
       end do
      else
       labels(1) = "ME_C0"
      endif
!     use reference given in input
      ncnt = is_keyword_set('calculate.reference')
      if (ncnt.GT.0) then
       if (oldref) call quit(1,'set_unc_mrci_targets',
     &      'We cannot use oldref and read the reference '//
     &      'from input at the same time.')
       if (ncnt.NE.n_states) call quit(1,'set_unc_mrci_targets',
     &      'Give one set of reference coefficients '//
     &      'for each state.')

       call set_rule2('SOLVE_REF',PRINT_,tgt_info)
       call set_arg('SOLVE_REF',PRINT_,'STRING',1,tgt_info,
     &      val_str='Reference coefficients from input.')

       do icnt=1,ncnt
        if (is_argument_set
     &         ('calculate.reference','def',keycount=icnt).NE.1)
     &         call quit(0,'set_unc_mrci_targets',
     &         'There must be exactly one def for each reference!')
        call set_rule2('SOLVE_REF',SET_STATE,tgt_info)
        call set_arg('SOLVE_REF',SET_STATE,'LISTS',1,tgt_info,
     &       val_label=['ME_C0'])
        call set_arg('SOLVE_REF',SET_STATE,'ISTATE',1,tgt_info,
     &       val_int=[icnt])
        call get_argument_dimension(arglen,'calculate.reference','def',
     &       keycount=icnt)
        if (allocated(xscr)) deallocate(xscr)
        allocate(xscr(arglen))
        call get_argument_value('calculate.reference','def',
     &       keycount=icnt,xarr=xscr)
        call set_rule2('SOLVE_REF',SET_MEL,tgt_info)
        call set_arg('SOLVE_REF',SET_MEL,'LIST',1,tgt_info,
     &       val_label=['ME_C0'])
        call set_arg('SOLVE_REF',SET_MEL,'IDX_LIST',arglen,tgt_info,
     &       val_int=[(i,i=1,arglen)])
        call set_arg('SOLVE_REF',SET_MEL,'VAL_LIST',arglen,tgt_info,
     &       val_rl8=xscr(1:arglen))
       end do

       call set_rule2('SOLVE_REF',SET_STATE,tgt_info)
       call set_arg('SOLVE_REF',SET_STATE,'LISTS',1,tgt_info,
     &      val_label=['ME_C0'])
       call set_arg('SOLVE_REF',SET_STATE,'ISTATE',1,tgt_info,
     &      val_int=[1])

      else if (.not.oldref) then
        call set_rule2('SOLVE_REF',PRINT_,tgt_info)
        call set_arg('SOLVE_REF',PRINT_,'STRING',1,tgt_info,
     &      val_str='Now solving equations for reference state')
        call set_rule2('SOLVE_REF',SOLVEEVP,tgt_info)
        call set_arg('SOLVE_REF',SOLVEEVP,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_C0'/))
        call set_arg('SOLVE_REF',SOLVEEVP,'N_ROOTS',1,tgt_info,
     &       val_int=(/maxroot/))
        call set_arg('SOLVE_REF',SOLVEEVP,'TARG_ROOT',1,tgt_info,
     &       val_int=(/ciroot/))
        call set_arg('SOLVE_REF',SOLVEEVP,'OP_MVP',1,tgt_info,
     &       val_label=(/'A_C0'/))
        call set_arg('SOLVE_REF',SOLVEEVP,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label2)/))
        call set_arg('SOLVE_REF',SOLVEEVP,'OP_SVP',1,tgt_info,
     &       val_label=(/'C0'/))
        call set_arg('SOLVE_REF',SOLVEEVP,'FORM',1,tgt_info,
     &       val_label=(/'FOPT_A_C0'/))
        if (spinadapt.eq.0.and.refproj.eq.0) then
          call set_arg('SOLVE_REF',SOLVEEVP,'MODE',1,tgt_info,
     &         val_str='DIA')
        else if (spinadapt.ne.0.and.refproj.eq.0) then
          call set_arg('SOLVE_REF',SOLVEEVP,'MODE',1,tgt_info,
     &         val_str='SPP')
          call set_arg('SOLVE_REF',SOLVEEVP,'LIST_SPC',1,tgt_info,
     &         val_label=(/'ME_C0_sp'/))
          call set_arg('SOLVE_REF',SOLVEEVP,'FORM_SPC',1,tgt_info,
     &         val_label=(/'FOPT_C0_sp'/))
        else if (spinadapt.ne.0.and.refproj.ne.0) then
          call set_arg('SOLVE_REF',SOLVEEVP,'MODE',1,tgt_info,
     &         val_str='SRP')
          call set_arg('SOLVE_REF',SOLVEEVP,'LIST_SPC',1,tgt_info,
     &         val_label=(/'ME_C0_sp'/))
          call set_arg('SOLVE_REF',SOLVEEVP,'FORM_SPC',2,tgt_info,
     &         val_label=(/'FOPT_C0_prj','FOPT_C0_sp '/))
        else ! refproj.ne.0.and.spinadapt.eq.0
          call set_arg('SOLVE_REF',SOLVEEVP,'MODE',1,tgt_info,
     &         val_str='PRJ')
          call set_arg('SOLVE_REF',SOLVEEVP,'FORM_SPC',1,tgt_info,
     &         val_label=(/'FOPT_C0_prj'/))
        end if
      else
       call set_rule2('SOLVE_REF',PRINT_,tgt_info)
       call set_arg('SOLVE_REF',PRINT_,'STRING',1,tgt_info,
     &      val_str='Reference coefficients from old file.')
       inquire(file='ME_C0_list.da',exist=l_exist)
       if (.not.l_exist) then
        if (.not.multistate)
     &       call quit(1,'set_unc_mrci_targets',
     &       'File for CASSCF coefficients not found!')
        call set_rule2('SOLVE_REF',SET_STATE,tgt_info)
        call set_arg('SOLVE_REF',SET_STATE,'LISTS',1,tgt_info,
     &       val_label=['ME_C0'])
        call set_arg('SOLVE_REF',SET_STATE,'ISTATE',1,tgt_info,
     &       val_int=[1])
        do i_state=1,n_states,1
         c_st = state_label(i_state,.true.)
         inquire(file='ME_C0'//trim(c_st)//'_list.da',exist=l_exist)
         if (.not.l_exist) call quit(1,'set_unc_mrci_targets',
     &        'File for CASSCF coefficients not found! State:'
     &        //trim(c_st))
         call set_rule2('SOLVE_REF',SCALE_COPY,tgt_info)
         call set_arg('SOLVE_REF',SCALE_COPY,'LIST_RES',1,tgt_info,
     &        val_label=['ME_C0'])
         call set_arg('SOLVE_REF',SCALE_COPY,'LIST_INP',1,tgt_info,
     &        val_label=['ME_C0'//trim(c_st)])
         call set_arg('SOLVE_REF',SCALE_COPY,'FAC',1,tgt_info,
     &        val_rl8=[1.0d0])
         call set_rule2('SOLVE_REF',SPREAD_MEL,tgt_info)
         call set_arg('SOLVE_REF',SPREAD_MEL,'LIST_IN',1,tgt_info,
     &        val_label=(/'ME_C0'/))
         if (n_states.GT.20) call quit(1,'set_unc_mrci_targets',
     &       'Static vector labels does not suport more than 20 states')
         call set_arg('SOLVE_REF',SPREAD_MEL,'LIST_OUT',i_state,
     &        tgt_info,val_label=labels)
         call set_rule2('SOLVE_REF',ADV_STATE,tgt_info)
         call set_arg('SOLVE_REF',ADV_STATE,'LISTS',1,tgt_info,
     &        val_label=['ME_C0'])
         call set_arg('SOLVE_REF',ADV_STATE,'N_ROOTS',1,tgt_info,
     &        val_int=[n_states])
        end do
       end if
c dbg
c        call set_rule2('SOLVE_REF',SET_MEL,tgt_info)
c        call set_arg('SOLVE_REF',SET_MEL,'LIST',1,tgt_info,
c     &       val_label=(/'ME_C0'/))
c        call set_arg('SOLVE_REF',SET_MEL,'IDX_LIST',2,tgt_info,
c     &       val_int=(/1,2/))
c        call set_arg('SOLVE_REF',SET_MEL,'VAL_LIST',2,tgt_info,
c     &       val_rl8=(/1d0,1d0/))
c dbgend
      end if
!     Spread ME_C0 for the states
      if (multistate) then
       call set_rule2('SOLVE_REF',SPREAD_MEL,tgt_info)
       call set_arg('SOLVE_REF',SPREAD_MEL,'LIST_IN',1,tgt_info,
     &      val_label=(/'ME_C0'/))
       if (n_states.GT.20) call quit(1,'set_unc_mrci_targets',
     &      'Static vector labels does not suport more than 20 states')
       call set_arg('SOLVE_REF',SPREAD_MEL,'LIST_OUT',n_states,
     &      tgt_info,val_label=labels)
      endif
      if (.false..and.cmaxexc.eq.0) then ! false: no long output
       do i_state=1,n_states,1
        if (.NOT.multistate) then
         call form_parameters(-1,parameters,2,
     &        'CI coefficients:',0,'LIST')
        else
         write(label_title,
     &        '("CI coefficients for state ",i0," (",A,"):")')
     &        i_state,trim(labels(i_state))
         call form_parameters(-1,parameters,2,
     &        label_title,0,'LIST')
        endif
        call set_rule('SOLVE_REF',ttype_opme,PRINT_MEL,
     &       labels(i_state),1,0,
     &       parameters,2,tgt_info)
       enddo
      else
!        call set_rule2('SOLVE_REF',ANALYZE_MEL,tgt_info)
!        call set_arg('SOLVE_REF',ANALYZE_MEL,'LISTS',1,tgt_info,
!     &               val_label=(/'ME_C0'/))
!        call set_arg('SOLVE_REF',ANALYZE_MEL,'LISTS_CV',1,tgt_info,
!     &               val_label=(/'ME_C0'/))
      end if

      ! Evaluate reference energy
      call add_target('EVAL_E_REF',ttype_gen,calc,tgt_info)
      if (spinexpec.gt.0.or.cmaxexc.eq.0) then
        call set_dependency('EVAL_E_REF','EVAL_REF_S(S+1)',tgt_info)
      else
        call set_dependency('EVAL_E_REF','SOLVE_REF',tgt_info)
      end if
      call set_dependency('EVAL_E_REF','FOPT_REF',tgt_info)
      call set_rule('EVAL_E_REF',ttype_opme,EVAL,
     &     'FOPT_REF',1,0,
     &     parameters,0,tgt_info)
      call set_rule2('EVAL_E_REF',PRINT_MEL,tgt_info)
      call set_arg('EVAL_E_REF',PRINT_MEL,'LIST',1,tgt_info,
     &     val_label=(/'ME_E_REF'/))
      call set_arg('EVAL_E_REF',PRINT_MEL,'COMMENT',1,tgt_info,
     &     val_str='Reference energy :')
      call set_arg('EVAL_E_REF',PRINT_MEL,'FORMAT',1,tgt_info,
     &     val_str='SCAL F20.12')
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Hamiltonian :',0,'LIST')
c      call set_rule('EVAL_E_REF',ttype_opme,PRINT_MEL,
c     &     mel_ham,1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Evaluate Fock operator wrt reference function
      call add_target('EVAL_FREF',ttype_gen,writeF,tgt_info)
      call set_dependency('EVAL_FREF','FOPT_FREF',tgt_info)
      if (spinexpec.gt.0.or.cmaxexc.eq.0) then
        call set_dependency('EVAL_FREF','EVAL_REF_S(S+1)',tgt_info)
      else
        call set_dependency('EVAL_FREF','SOLVE_REF',tgt_info)
      end if
      do i_state = 1,n_states
      c_st = state_label(i_state,.false.)
      call set_rule('EVAL_FREF',ttype_opme,EVAL,
     &     'FOPT_FREF',1,0,
     &     parameters,0,tgt_info)
      if (writeF) then ! write to file 'FEFF'
        call form_parameters(-1,parameters,2,
     &       'effective Fock operator:',0,'FEFF')
        call set_rule('EVAL_FREF',ttype_opme,PRINT_MEL,
     &       'ME_FREF',1,0,
     &       parameters,2,tgt_info)
      end if
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'effective Fock operator:',0,'LIST')
c      call set_rule('EVAL_FREF',ttype_opme,PRINT_MEL,
c     &     'ME_FREF',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      if(multistate)then
       call set_rule2('EVAL_FREF',ADV_STATE,tgt_info)
       call set_arg('EVAL_FREF',ADV_STATE,'LISTS',2,tgt_info,
     &      val_label=['ME_C0  ',
     &                 'ME_FREF'])
       call set_arg('EVAL_FREF',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
      end if
      end do

      ! Evaluate Fock operator wrt zeroth order reference function
      call add_target('EVAL_FREF0',ttype_gen,writeF,tgt_info)
      call set_dependency('EVAL_FREF0','FOPT_FREF0',tgt_info)
C      if (spinexpec.gt.0.or.cmaxexc.eq.0) then
C        call set_dependency('EVAL_FREF0','EVAL_REF_S(S+1)',tgt_info)
C      else
C        call set_dependency('EVAL_FREF0','SOLVE_REF',tgt_info)
C      end if
      do i_state = 1,n_states
      c_st = state_label(i_state,.false.)
      call set_rule('EVAL_FREF0',ttype_opme,EVAL,
     &     'FOPT_FREF0',1,0,
     &     parameters,0,tgt_info)
      if (writeF) then ! write to file 'FEFF'
        call form_parameters(-1,parameters,2,
     &       'effective Fock operator:',0,'FEFF0')
        call set_rule('EVAL_FREF0',ttype_opme,PRINT_MEL,
     &       'ME_FREF0',1,0,
     &       parameters,2,tgt_info)
      end if
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'effective Fock operator:',0,'LIST')
c      call set_rule('EVAL_FREF0',ttype_opme,PRINT_MEL,
c     &     'ME_FREF0',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      if(multistate)then
       call set_rule2('EVAL_FREF0',ADV_STATE,tgt_info)
       call set_arg('EVAL_FREF0',ADV_STATE,'LISTS',2,tgt_info,
     &      val_label=['ME_C0   ',
     &                 'ME_FREF0'])
       call set_arg('EVAL_FREF0',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
      end if
      end do

      ! Evaluate spin expectation value
      call add_target('EVAL_REF_S(S+1)',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_REF_S(S+1)','SOLVE_REF',tgt_info)
      call set_dependency('EVAL_REF_S(S+1)','FOPT_REF_S(S+1)',tgt_info)
      do i_state = 1,n_states
      c_st = state_label(i_state,.false.)
      call set_rule('EVAL_REF_S(S+1)',ttype_opme,EVAL,
     &     'FOPT_REF_S(S+1)',1,0,
     &     parameters,0,tgt_info)
      call set_rule2('EVAL_REF_S(S+1)',PRINT_MEL,tgt_info)
      call set_arg('EVAL_REF_S(S+1)',PRINT_MEL,'LIST',1,tgt_info,
     &     val_label=(/'ME_S(S+1)'/))
      call set_arg('EVAL_REF_S(S+1)',PRINT_MEL,'COMMENT',1,tgt_info,
     &     val_str='Spin expectation value <C0| S^2 |C0> :')
      call set_arg('EVAL_REF_S(S+1)',PRINT_MEL,'FORMAT',1,tgt_info,
     &     val_str='SCAL F20.12')
      call set_arg('EVAL_REF_S(S+1)',PRINT_MEL,'CHECK_THRESH',1,
     &     tgt_info,val_rl8=(/1d-2/))
      call set_arg('EVAL_REF_S(S+1)',PRINT_MEL,'EXPECTED',1,tgt_info,
     &     val_rl8=(/(dble(orb_info%imult**2)-1d0)/4d0/))
      if(multistate)then
       call set_rule2('EVAL_REF_S(S+1)',ADV_STATE,tgt_info)
       call set_arg('EVAL_REF_S(S+1)',ADV_STATE,'LISTS',2,tgt_info,
     &      val_label=['ME_C0    ',
     &                 'ME_S(S+1)'])
       call set_arg('EVAL_REF_S(S+1)',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
      end if
      end do

      ! prepare initial guess from CAS-CI
      call add_target2('C0guess',.false.,tgt_info)
      call set_dependency('C0guess','DEF_ME_C0',tgt_info)
      ! (a) define operator and ME list for old CAS-CI coefficients
      occ_def = 0
      occ_def(IVALE,1,1) = orb_info%nactel
      call set_rule2('C0guess',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('C0guess',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'C0start'/))
      call set_arg('C0guess',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('C0guess',DEF_OP_FROM_OCC,'OCC',1,tgt_info,
     &     val_occ=occ_def(1:ngastp,1:2,1))
      call set_rule2('C0guess',DEF_ME_LIST,tgt_info)
      call set_arg('C0guess',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C0start'/))
      call set_arg('C0guess',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C0start'/))
      call set_arg('C0guess',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('C0guess',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('C0guess',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_arg('C0guess',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('C0guess',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/guess/))
      call set_arg('C0guess',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/guess/))
      ! (b) define, optimize and evaluate formula C0 = C0start
      call set_rule2('C0guess',DEF_FORMULA,tgt_info)
      call set_arg('C0guess',DEF_FORMULA,'LABEL',1,tgt_info,
     &             val_label=(/'F_C0guess'/))
      call set_arg('C0guess',DEF_FORMULA,'FORMULA',1,tgt_info,
     &             val_str='C0=C0start')
      call set_rule2('C0guess',OPTIMIZE,tgt_info)
      call set_arg('C0guess',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_C0guess'/))
      call set_arg('C0guess',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_C0guess'/))
      call set_rule2('C0guess',EVAL,tgt_info)
      call set_arg('C0guess',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_C0guess'/))

      ! another target for restart purposes (of non-linear equations)
      call add_target2('C0rst',.false.,tgt_info)
      call set_dependency('C0rst','DEF_ME_C0',tgt_info)
      call set_dependency('C0rst','DEF_ME_C00',tgt_info) ! save C0 from prev. run
      ! (a) define operator and ME list for old CAS-CI coefficients
      occ_def = 0
      occ_def(IVALE,1,1) = orb_info%nactel
      call set_rule2('C0rst',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('C0rst',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'C0rst'/))
      call set_arg('C0rst',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('C0rst',DEF_OP_FROM_OCC,'OCC',1,tgt_info,
     &     val_occ=occ_def(1:ngastp,1:2,1))
      call set_rule2('C0rst',DEF_ME_LIST,tgt_info)
      call set_arg('C0rst',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C0rst'/))
      call set_arg('C0rst',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C0rst'/))
      call set_arg('C0rst',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('C0rst',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('C0rst',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_arg('C0rst',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('C0rst',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('C0rst',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))
      ! (b) define, optimize and evaluate formula C0 = C0rst
      call set_rule2('C0rst',DEF_FORMULA,tgt_info)
      call set_arg('C0rst',DEF_FORMULA,'LABEL',1,tgt_info,
     &             val_label=(/'F_C0rst'/))
      call set_arg('C0rst',DEF_FORMULA,'FORMULA',1,tgt_info,
     &             val_str='C0=C0rst')
      call set_rule2('C0rst',OPTIMIZE,tgt_info)
      call set_arg('C0rst',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_C0rst'/))
      call set_arg('C0rst',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_C0rst'/))
      call set_rule2('C0rst',EVAL,tgt_info)
      call set_arg('C0rst',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_C0rst'/))

      return
      end
