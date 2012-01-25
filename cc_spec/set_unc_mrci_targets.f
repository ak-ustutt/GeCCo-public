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
      include 'ifc_targets.h'
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
     &     ndef, occ_def(ngastp,2,124),!60),
     &     isym, msc, ims, ip, ih, 
     &     cminh, cmaxh, cminp, cmaxp, cmaxexc, ciroot, cmaxv
      logical ::
     &     oldref, l_exist
      character(len_target_name) ::
     &     dia_label, labels(20)
      character(len_command_par) ::
     &     parameters(3)

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting targets for multiref. wave function'

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
      call get_argument_value('method.MR','ciroot',
     &     ival=ciroot)
      call get_argument_value('method.MR','oldref',
     &     lval=oldref)
      if (cmaxh.lt.0) cmaxh = cmaxexc
      if (cmaxp.lt.0) cmaxp = cmaxexc
      cmaxv = orb_info%norb_hpv(IVALE,1)*2

      if (ntest.ge.100) then
        write(luout,*) 'cminh   = ',cminh
        write(luout,*) 'cmaxh   = ',cmaxh
        write(luout,*) 'cminp   = ',cminp
        write(luout,*) 'cmaxp   = ',cmaxp
        write(luout,*) 'cmaxv   = ',cmaxv
        write(luout,*) 'cmaxexc = ',cmaxexc
        write(luout,*) 'nactel  = ',orb_info%nactel
        write(luout,*) 'ciroot  = ',ciroot
        write(luout,*) 'oldref  = ',oldref
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
      call set_arg('FREF',DEF_HAMILTONIAN,'X_SPCS',2,tgt_info,
     &     val_int=(/IVALE,IEXTR/))

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
      ! Sz_dum
      call add_target('Sz_dum',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('Sz_dum','S+',tgt_info)
      call cloneop_parameters(-1,parameters,'S+',.false.)
      call set_rule('Sz_dum',ttype_op,CLONE_OP,'Sz_dum',1,1,
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
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_FREF',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! expectation value of S^2
      call add_target2('F_REF_S(S+1)',.false.,tgt_info)
      call set_dependency('F_REF_S(S+1)','S(S+1)',tgt_info)
      call set_dependency('F_REF_S(S+1)','C0',tgt_info)
      call set_dependency('F_REF_S(S+1)','S+',tgt_info)
      call set_dependency('F_REF_S(S+1)','S-',tgt_info)
      call set_dependency('F_REF_S(S+1)','Sz',tgt_info)
      call set_dependency('F_REF_S(S+1)','Sz_dum',tgt_info)
      ! (a) 1/2*(S+S- + S-S+)
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','S+','S-','C0'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','S-','S+','C0'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! (b) + Sz^2 (Sz_dum is used to circumvent automatic "BCH" factor)
      call set_rule2('F_REF_S(S+1)',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','Sz','Sz_dum','C0'/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
      call set_arg('F_REF_S(S+1)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_REF_S(S+1)',REPLACE,tgt_info)
      call set_arg('F_REF_S(S+1)',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_REF_S(S+1)'/))
      call set_arg('F_REF_S(S+1)',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'Sz_dum','Sz'/))
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
      call set_dependency('F_REPL_H_S2','Sz_dum',tgt_info)
      call set_dependency('F_REPL_H_S2','H',tgt_info)
      ! (a) 1/2*(S+S- + S-S+)
      call set_rule2('F_REPL_H_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'H','S+','S-','H'/))
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
     &     val_label=(/'H','S-','S+','H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! (b) + Sz^2 (Sz_dum is used to circumvent automatic "BCH" factor)
      call set_rule2('F_REPL_H_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'H','Sz','Sz_dum','H'/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('F_REPL_H_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_REPL_H_S2',REPLACE,tgt_info)
      call set_arg('F_REPL_H_S2',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_REPL_H_S2',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'Sz_dum','Sz'/))
c dbg
c      call set_rule2('F_REPL_H_S2',PRINT_FORMULA,tgt_info)
c      call set_arg('F_REPL_H_S2',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_REPL_H_S2'/))
c dbgend
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
      call set_dependency('FOPT_REF_S(S+1)','DEF_ME_Sz',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_REF_S(S+1)',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

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
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/ims/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/orb_info%lsym/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/ciroot/))
      call set_arg('DEF_ME_C0',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/ciroot/))

      ! ME_A_C0
      call add_target('DEF_ME_A_C0',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_A_C0','A_C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_A_C0'
      labels(2) = 'A_C0'
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
     &                    op_dia//'_'//'C0',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C0'
      labels(2) = op_dia//'_'//'C0'
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
      call set_arg('DEF_ME_S+',DEF_ME_LIST,'MS',1,tgt_info,
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
      call set_arg('DEF_ME_S-',DEF_ME_LIST,'MS',1,tgt_info,
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
      call set_arg('DEF_ME_Sz',DEF_ME_LIST,'MS',1,tgt_info,
     &     val_int=(/0/))
      call set_rule2('DEF_ME_Sz',UNITY,tgt_info)
      call set_arg('DEF_ME_Sz',UNITY,'LIST',1,tgt_info,
     &             val_label=(/'ME_Sz'/))
      call set_arg('DEF_ME_Sz',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Sz',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.true./))

      ! ME_S(S+1)
      call add_target('DEF_ME_S(S+1)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_S(S+1)','S(S+1)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_S(S+1)'
      labels(2) = 'S(S+1)'
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DEF_ME_S(S+1)',ttype_opme,DEF_ME_LIST,
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
      call solve_parameters(-1,parameters,2,1,ciroot,'DIA')
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_C0'
      labels(2) = trim(dia_label)//'C0'
      labels(3) = 'A_C0'
      labels(4) = 'C0'
      labels(5) = 'FOPT_A_C0'
      if (.not.oldref) then
        call set_rule('SOLVE_REF',ttype_opme,SOLVEEVP,
     &     labels,5,1,
     &     parameters,2,tgt_info)
      else
        inquire(file='ME_C0_list.da',exist=l_exist)
        if (.not.l_exist) call quit(1,'set_unc_mrci_targets',
     &           'File for CASSCF coefficients not found!')
      end if
      if (cmaxexc.eq.0) then
        call form_parameters(-1,parameters,2,
     &       'CI coefficients :',0,'LIST')
        call set_rule('SOLVE_REF',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
      end if

      ! Evaluate reference energy
      call add_target('EVAL_E_REF',ttype_gen,calc,tgt_info)
      call set_dependency('EVAL_E_REF','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('EVAL_E_REF','FOPT_REF',tgt_info)
      call set_rule('EVAL_E_REF',ttype_opme,EVAL,
     &     'FOPT_REF',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Hamiltonian :',0,'LIST')
c      call set_rule('EVAL_E_REF',ttype_opme,PRINT_MEL,
c     &     mel_ham,1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Evaluate Fock operator wrt reference function
      call add_target('EVAL_FREF',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_FREF','FOPT_FREF',tgt_info)
      call set_dependency('EVAL_FREF','EVAL_REF_S(S+1)',tgt_info)
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

      ! Evaluate spin expectation value
      call add_target('EVAL_REF_S(S+1)',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_REF_S(S+1)','SOLVE_REF',tgt_info)
      call set_dependency('EVAL_REF_S(S+1)','FOPT_REF_S(S+1)',tgt_info)
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

      return
      end
