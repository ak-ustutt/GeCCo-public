*----------------------------------------------------------------------*
      subroutine set_experimental_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets for experiments with GeCCo
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'opdim.h'

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
     &     isym, ms, msc, sym_arr(8),
     &     occ_def(ngastp,2,20), ndef, n_pp
      logical ::
     &     needed
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(12) ::
     &     approx

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting experimental targets ...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*
      ! e.g. 
      call add_target('MY_OP',ttype_op,.false.,tgt_info)
      call set_rule('MY_OP',ttype_op,DEF_SCALAR,
     &              'MY_OP',1,1,
     &              parameters,0,tgt_info)
      ! cf. set_xxxx_targets.f routines in e.g. cc_special for
      ! further examples
      call add_target('FOCK',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   1,1,1,.false.)
      call set_rule('FOCK',ttype_op,DEF_HAMILTONIAN,
     &              'FOCK',1,1,
     &              parameters,1,tgt_info)
      call add_target('PHI',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   2,2,1,.false.)
      call set_rule('PHI',ttype_op,DEF_HAMILTONIAN,
     &              'PHI',1,1,
     &              parameters,1,tgt_info)
      call add_target('V(1)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   1,1,1,.false.)
      call set_rule('V(1)',ttype_op,DEF_HAMILTONIAN,
     &              'V(1)',1,1,
     &              parameters,1,tgt_info)
      call add_target('VXX(1)',ttype_op,.false.,tgt_info)
      if (orb_info%ngas.eq.4) then
        call genop_parameters(-1,parameters,3,
     &                  1,1,0,
c     &                 -1,1,3,(/.true.,.false./),
     &                 -1,1,3,(/.false.,.false./),
     &                (/0,2, 0,2, 0,0, 0,2/),
     &                (/0,1, 0,1, 0,0, 0,1,
     &                  0,1, 0,1, 0,0, 0,1/),4,
     &                (/0,1, 0,1, 0,1, 0,1,
     &                  0,1, 0,1, 0,1, 0,1,
     &                  0,0, 0,0, 0,0, 0,0,
     &                  0,0, 0,0, 0,0, 0,0/),4
     &              )
      else if (orb_info%ngas.eq.3) then
        call genop_parameters(-1,parameters,3,
     &                  1,1,0,
     &                 -1,1,3,(/.false.,.false./),
     &                (/0,2, 0,2, 0,0, 0,2/),
     &                (/0,1, 0,1, 0,0, 0,1,
     &                  0,1, 0,1, 0,0, 0,1/),4,
     &                (/0,1, 0,1, 0,1,
     &                  0,1, 0,1, 0,1,
     &                  0,0, 0,0, 0,0,
     &                  0,0, 0,0, 0,0/),3
     &              )
      else
        call quit(1,'experiment: ','ngas.ne.3/4: test does not work')
      end if
      call set_rule('VXX(1)',ttype_op,DEF_GENERAL_OPERATOR,
     &              'VXX(1)',1,1,
     &              parameters,3,tgt_info)

      call add_target('VX(1)',ttype_op,.false.,tgt_info)
      if (orb_info%ngas.eq.4) then
        call genop_parameters(-1,parameters,3,
     &                  1,1,0,
     &                 -1,1,3,(/.true.,.false./),
c     &                 -1,1,3,(/.false.,.false./),
     &                (/0,2, 0,2, 0,0, 0,2/),
     &                (/0,1, 0,1, 0,0, 0,1,
     &                  0,1, 0,1, 0,0, 0,1/),4,
     &                (/0,0, 0,1, 0,1, 0,1,
     &                  0,1, 0,1, 0,1, 0,1,
     &                  0,0, 0,0, 0,0, 0,0,
     &                  0,0, 0,0, 0,0, 0,0/),4
     &              )
      else if (orb_info%ngas.eq.3) then
        call genop_parameters(-1,parameters,3,
     &                  1,1,0,
     &                 -1,1,3,(/.false.,.false./),
     &                (/0,2, 0,2, 0,0, 0,2/),
     &                (/0,1, 0,1, 0,0, 0,1,
     &                  0,1, 0,1, 0,0, 0,1/),4,
     &                (/0,1, 0,1, 0,1,
     &                  0,1, 0,1, 0,1,
     &                  0,0, 0,0, 0,0,
     &                  0,0, 0,0, 0,0/),3
     &              )
      else
        call quit(1,'experiment: ','ngas.ne.3/4: test does not work')
      end if
      call set_rule('VX(1)',ttype_op,DEF_GENERAL_OPERATOR,
     &              'VX(1)',1,1,
     &              parameters,3,tgt_info)

      call add_target('V(A)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   1,1,1,.false.)
      call set_rule('V(A)',ttype_op,DEF_HAMILTONIAN,
     &              'V(A)',1,1,
     &              parameters,1,tgt_info)
      call add_target('V(B)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,
     &                   1,1,1,.false.)
      call set_rule('V(B)',ttype_op,DEF_HAMILTONIAN,
     &              'V(B)',1,1,
     &              parameters,1,tgt_info)
      call add_target('H(0)',ttype_op,.false.,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     'H',.false.)
      call set_rule('H(0)',ttype_op,CLONE_OP,
     &              'H(0)',1,1,
     &              parameters,1,tgt_info)

      call add_target('T(A)',ttype_op,.false.,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     'T',.false.)
      call set_rule('T(A)',ttype_op,CLONE_OP,
     &              'T(A)',1,1,
     &              parameters,1,tgt_info)
      call add_target('T(B)',ttype_op,.false.,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     'T',.false.)
      call set_rule('T(B)',ttype_op,CLONE_OP,
     &              'T(B)',1,1,
     &              parameters,1,tgt_info)

      call add_target('LMP2',ttype_op,.false.,tgt_info)
      call set_rule('LMP2',ttype_op,DEF_SCALAR,
     &              'LMP2',1,1,
     &              parameters,0,tgt_info)
      call add_target('EMP2',ttype_op,.false.,tgt_info)
      call set_rule('EMP2',ttype_op,DEF_SCALAR,
     &              'EMP2',1,1,
     &              parameters,0,tgt_info)


      call add_target('T2',ttype_op,.false.,tgt_info)
      call xop_parameters(-1,parameters,
     &                   .false.,2,2,0,1)
      call set_rule('T2',ttype_op,DEF_EXCITATION,
     &              'T2',1,1,
     &              parameters,1,tgt_info)
      call add_target('O2',ttype_op,.false.,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     'T2',.false.)
      call set_rule('O2',ttype_op,CLONE_OP,
     &              'O2',1,1,
     &              parameters,1,tgt_info)


      call add_target('BV',ttype_op,.false.,tgt_info)
      call set_dependency('BV',op_b_inter,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_b_inter,.false.)
      call set_rule('BV',ttype_op,CLONE_OP,
     &              'BV',1,1,
     &              parameters,1,tgt_info)
      
      call add_target('RDGB(V)',ttype_op,.false.,tgt_info)
      call set_dependency('RDGB(V)',op_rint,tgt_info)
      call cloneop_parameters(-1,parameters,
     &     op_rint,.false.)
      call set_rule('RDGB(V)',ttype_op,CLONE_OP,
     &              'RDGB(V)',1,1,
     &              parameters,1,tgt_info)
      

      n_pp = 0
      occ_def = 0
      ndef = 2
      ! 1
      occ_def(IHOLE,1,1) = 1
      occ_def(IPART,1,1) = 1
      occ_def(IHOLE,2,1) = 2
      ! 2
      occ_def(IHOLE,1,2) = 1
      occ_def(IEXTR,1,2) = 1
      occ_def(IHOLE,2,2) = 2
      if (n_pp.eq.1) then
        ndef = 4
        ! 3
        occ_def(IHOLE,1,3) = 1
        occ_def(IPART,1,3) = 1
        occ_def(IHOLE,2,3) = 1
        occ_def(IPART,2,3) = 1
        ! 4
        occ_def(IHOLE,1,4) = 1
        occ_def(IEXTR,1,4) = 1
        occ_def(IHOLE,2,4) = 1
        occ_def(IPART,2,4) = 1
      end if
      call add_target('RBRV(V)',ttype_op,.false.,tgt_info)
      call op_from_occ_parameters(-1,parameters,2,
     &     occ_def,ndef,1,(/.false.,.true./),ndef)
      call set_rule('RBRV(V)',ttype_op,DEF_OP_FROM_OCC,
     &              'RBRV(V)',1,1,
     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*
      ! e.g.
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_FORM'
      labels(2) = 'LMP2'
      labels(3) = 'FOCK'
      labels(4) = 'PHI'
      labels(5) = 'T2'
      call add_target('MP2_FORM',ttype_frm,.false.,tgt_info)
      call set_dependency('MP2_FORM','LMP2',tgt_info)
      call set_dependency('MP2_FORM','FOCK',tgt_info)
      call set_dependency('MP2_FORM','PHI',tgt_info)
      call set_dependency('MP2_FORM','T2',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'title of my new formula',0,'mode string')
      call set_rule('MP2_FORM',ttype_frm,DEF_EXP_FORMULA,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_FORM'
      labels(2) = 'MP2_FORM'
      labels(3) = 'FOCK'
      labels(4) = 'H'
      labels(5) = 'PHI'
      labels(6) = 'H'
      call set_dependency('MP2_FORM','H',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'same title?',2,'---')
      call set_rule('MP2_FORM',ttype_frm,REPLACE,
     &     labels,6,1,
     &     parameters,2,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_RES'
      labels(2) = 'MP2_FORM'
      labels(3) = 'O2'
      labels(4) = 'T2^+'
      labels(5) = ' '
      call add_target('MP2_RES',ttype_frm,.false.,tgt_info)
      call set_dependency('MP2_RES','MP2_FORM',tgt_info)
      call set_dependency('MP2_RES','O2',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'residual',1,'---')
      call set_rule('MP2_RES',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call set_rule('MP2_RES',ttype_frm,TEX_FORMULA,
c     &              labels,1,0,
c     &              'mp2res.tex',1,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_EN'
      labels(2) = 'MP2_FORM'
      labels(3) = 'MP2_RES'
      call add_target('MP2_EN',ttype_frm,.false.,tgt_info)
      call set_dependency('MP2_EN','MP2_FORM',tgt_info)
      call set_dependency('MP2_EN','MP2_RES',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'energy',1,'---')
      call set_rule('MP2_EN',ttype_frm,FACTOR_OUT,
     &              labels,3,1,
     &              parameters,2,tgt_info)
     

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'H_FORM'
      call add_target('H_FORM',ttype_frm,.false.,tgt_info)
      call set_dependency('H_FORM','H',tgt_info)
      call set_dependency('H_FORM','H(0)',tgt_info)
      call set_dependency('H_FORM','V(1)',tgt_info)
      call def_form_parameters(-1,
     &     parameters,2,'H=H(0)+V(1)','pert exp of H')
      call set_rule('H_FORM',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'VV_FORM'
      call add_target('VV_FORM',ttype_frm,.false.,tgt_info)
      call set_dependency('VV_FORM','V(1)',tgt_info)
      call set_dependency('VV_FORM','V(A)',tgt_info)
      call set_dependency('VV_FORM','V(B)',tgt_info)
      call def_form_parameters(-1,
     &     parameters,2,'V(1)=V(A)+V(B)','pert exp of H')
      call set_rule('VV_FORM',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'TT_FORM'
      call add_target('TT_FORM',ttype_frm,.false.,tgt_info)
      call set_dependency('TT_FORM','T',tgt_info)
      call set_dependency('TT_FORM','T(A)',tgt_info)
      call set_dependency('TT_FORM','T(B)',tgt_info)
      call def_form_parameters(-1,
     &     parameters,2,'T=T(A)+T(B)','pert exp of T')
      call set_rule('TT_FORM',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_FORM2'
      labels(2) = 'MP2_FORM'
      labels(3) = 'H_FORM'
      call add_target('MP2_FORM2',ttype_frm,.false.,tgt_info)
      call set_dependency('MP2_FORM2','MP2_FORM',tgt_info)
      call set_dependency('MP2_FORM2','H_FORM',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'title of my new formula',1,'mode string')
      call set_rule('MP2_FORM2',ttype_frm,EXPAND,
     &              labels,3,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('MP2_FORM2',ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)

c      call add_target('CC_FORM1',ttype_frm,.false.,tgt_info)
c      call set_dependency('CC_FORM1',form_cclg0,tgt_info)
c      call set_dependency('CC_FORM1','H_FORM',tgt_info)
c      call set_dependency('CC_FORM1','VV_FORM',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'CC_FORM1'
c      labels(2) = form_cclg0
c      labels(3) = 'H_FORM'
c      call form_parameters(-1,
c     &     parameters,2,'title of my new formula',1,'mode string')
c      call set_rule('CC_FORM1',ttype_frm,EXPAND,
c     &              labels,3,1,
c     &              parameters,2,tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'CC_FORM1'
c      labels(2) = 'CC_FORM1'
c      labels(3) = 'VV_FORM'
c      call form_parameters(-1,
c     &     parameters,2,'title of my new formula',1,'mode string')
c      call set_rule('CC_FORM1',ttype_frm,EXPAND,
c     &              labels,3,1,
c     &              parameters,2,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule('CC_FORM1',ttype_frm,PRINT_FORMULA,
c     &              labels,1,0,
c     &              parameters,2,tgt_info)

c      call add_target('CC_FORM2',ttype_frm,.false.,tgt_info)
c      call set_dependency('CC_FORM2',form_cclg0,tgt_info)
c      call set_dependency('CC_FORM2','H_FORM',tgt_info)
c      call set_dependency('CC_FORM2','TT_FORM',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'CC_FORM2'
c      labels(2) = form_cctbar_a
c      labels(3) = 'TT_FORM'
c      call form_parameters(-1,
c     &     parameters,2,'title of my new formula',1,'mode string')
c      call set_rule('CC_FORM2',ttype_frm,EXPAND,
c     &              labels,3,1,
c     &              parameters,2,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule('CC_FORM2',ttype_frm,PRINT_FORMULA,
c     &              labels,1,0,
c     &              parameters,2,tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'CC_FORM2'
c      labels(2) = 'CC_FORM2'
c      labels(3) = 'VV_FORM'
c      call form_parameters(-1,
c     &     parameters,2,'title of my new formula',1,'mode string')
c      call set_rule('CC_FORM2',ttype_frm,EXPAND,
c     &              labels,3,1,
c     &              parameters,2,tgt_info)
c      call form_parameters(-1,
c     &     parameters,2,'stdout',0,'---')
c      call set_rule('CC_FORM2',ttype_frm,PRINT_FORMULA,
c     &              labels,1,0,
c     &              parameters,2,tgt_info)


      call add_target('BV_formal',ttype_frm,.false.,tgt_info)
      call set_dependency('BV_formal','VX(1)',tgt_info)
      call set_dependency('BV_formal','BV',tgt_info)
      call set_dependency('BV_formal',op_r12,tgt_info)
      labels(1) = 'BV_formal'
      labels(2) = 'BV'
      labels(3) = op_r12
      labels(4) = 'VX(1)'
      call form_parameters(-1,
     &     parameters,2,'BV formal',0,'Bp')
      call set_rule('BV_formal',ttype_frm,DEF_R12INTM_FORMAL,
     &              labels,4,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('BV_formal',ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)

      call add_target('RDGB(V)_FRM',ttype_frm,.false.,tgt_info)
      call set_dependency('RDGB(V)_FRM','VXX(1)',tgt_info)
      call set_dependency('RDGB(V)_FRM',op_rint,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RDGB(V)_FRM'
      labels(2) = 'RDGB(V)'
      labels(3) = op_rint
      labels(4) = '-'
      labels(5) = 'VXX(1)'
      call form_parameters(-1,
     &     parameters,2,'RBAR+(V)',3,'RD')
      call set_rule('RDGB(V)_FRM',ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('RDGB(V)_FRM',ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)

      call add_target('RBRV(V)_FRM',ttype_frm,.false.,tgt_info)
      call set_dependency('RBRV(V)_FRM','VXX(1)',tgt_info)
      call set_dependency('RBRV(V)_FRM',op_rint,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RBRV(V)_FRM'
      labels(2) = 'RBRV(V)'
      labels(3) = op_rint
      labels(4) = '-'
      labels(5) = 'VXX(1)'
      call form_parameters(-1,
     &     parameters,2,'RBRV(V)',3,'RV')
      call set_rule('RBRV(V)_FRM',ttype_frm,DEF_R12INTM_CABS,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('RBRV(V)_FRM',ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)

      call add_target('BV_CABS',ttype_frm,.false.,tgt_info)
      call set_dependency('BV_CABS','VX(1)',tgt_info)
      call set_dependency('BV_CABS','BV',tgt_info)
      call set_dependency('BV_CABS',op_rint,tgt_info)
      call set_dependency('BV_CABS','RDGB(V)',tgt_info)
      call set_dependency('BV_CABS','RBRV(V)',tgt_info)
      call set_dependency('BV_CABS','FF-X',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'BV_CABS'
      labels(2) = 'BV'
      labels(3) = op_rint
      labels(4) = 'RDGB(V)'
      labels(5) = 'FF-X'
      labels(6) = 'VX(1)'
      labels(7) = 'RBRV(V)'
      approx(1:12) = ' '
      approx(1:1) = 'C'
      approx(12:12) = 'S'
      call form_parameters(-1,
     &     parameters,2,'BV CABS',3,'BV'//approx)
      call set_rule('BV_CABS',ttype_frm,DEF_R12INTM_CABS,
     &              labels,7,1,
     &              parameters,2,tgt_info)
      call form_parameters(-1,
     &     parameters,2,'stdout',0,'---')
      call set_rule('BV_CABS',ttype_frm,PRINT_FORMULA,
     &              labels,1,0,
     &              parameters,2,tgt_info)

c*----------------------------------------------------------------------*
c*     Opt. Formulae 
c*----------------------------------------------------------------------*
c
c      ! e.g.:
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'MP2_OPT'
      labels(2) = 'MP2_RES'
      labels(3) = 'MP2_EN'
      ncat = 2  ! 2 formulae pasted into final formula
      nint = 0  ! no intermediate to factor out so far ...
      call add_target('MP2_OPT',ttype_frm,.false.,tgt_info)
      call set_dependency('MP2_OPT','MP2_EN',tgt_info)
      call set_dependency('MP2_OPT','MP2_RES',tgt_info)
      call set_dependency('MP2_OPT','DEF_ME_LMP2',tgt_info)
      call set_dependency('MP2_OPT','DEF_ME_O2',tgt_info)
      call set_dependency('MP2_OPT','DEF_ME_T2',tgt_info)
      call set_dependency('MP2_OPT','H0',tgt_info)
c      call set_dependency('MY_OPT','MY_MELR',tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('MP2_OPT',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)
c

      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'B(VZ)_OPT'
      labels(2) = 'RDGB(V)_FRM'
      labels(3) = 'RBRV(V)_FRM'
      labels(4) = 'BV_CABS'
      ncat = 3  ! 3 formulae pasted into final formula
      nint = 0  ! no intermediate to factor out so far ...
      call add_target('B(VZ)_OPT',ttype_frm,.false.,tgt_info)
      call set_dependency('B(VZ)_OPT','BV_CABS',tgt_info)
      call set_dependency('B(VZ)_OPT','RDGB(V)_FRM',tgt_info)
      call set_dependency('B(VZ)_OPT','RBRV(V)_FRM',tgt_info)
      call set_dependency('B(VZ)_OPT','DEF-B(VZ)',tgt_info)
      call set_dependency('B(VZ)_OPT','DEF-RDGB(VZ)',tgt_info)
      call set_dependency('B(VZ)_OPT','DEF-RBRV(VZ)',tgt_info)
      call set_dependency('B(VZ)_OPT',mel_rint,tgt_info)
      call set_dependency('B(VZ)_OPT','VX(1)LIST',tgt_info)
      call set_dependency('B(VZ)_OPT','VXX(1)LIST',tgt_info)
      call set_dependency('B(VZ)_OPT','FF-X-INT',tgt_info)
c      call set_dependency('MY_OPT','MY_MELR',tgt_info)      
      call opt_parameters(-1,parameters,ncat,nint)
      call set_rule('B(VZ)_OPT',ttype_frm,OPTIMIZE,
     &              labels,ncat+nint+1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      call add_target('DEF_ME_LMP2',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_LMP2','LMP2',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_LMP2'
      labels(2) = 'LMP2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF_ME_LMP2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target('DEF_ME_T2',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_T2','T2',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_T2'
      labels(2) = 'T2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF_ME_T2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call add_target('DEF_ME_O2',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_O2','O2',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_O2'
      labels(2) = 'O2'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF_ME_O2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      call add_target('DIAG',ttype_opme,.false.,tgt_info)
      call set_dependency('DIAG','H0',tgt_info)
      call set_dependency('DIAG','T2',tgt_info)
      call cloneop_parameters(-1,parameters,
     &     'T2',.false.)
      call set_rule('DIAG',ttype_op,CLONE_OP,
     &              'D2',1,1,
     &              parameters,1,tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'DIAG'
      labels(2) = 'D2'
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DIAG',ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      labels(1) = 'DIAG'
      labels(2) = 'H0'
      call set_rule('DIAG',ttype_opme,PRECONDITIONER,
     &              labels,2,1,
     &              parameters,0,tgt_info)

      call add_target('V(1)LIST',ttype_opme,.false.,tgt_info)
      call set_dependency('V(1)LIST','V(1)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'V(1)LIST'
      labels(2) = 'V(1)'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('V(1)LIST',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'V(1)LIST'
      call import_parameters(-1,parameters,'ZDIPLEN','DALTON')
      call set_rule('V(1)LIST',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('VXX(1)LIST',ttype_opme,.false.,tgt_info)
      call set_dependency('VXX(1)LIST','VXX(1)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'VXX(1)LIST'
      labels(2) = 'VXX(1)'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('VXX(1)LIST',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'VXX(1)LIST'
      call import_parameters(-1,parameters,'ZDIPLEN','DALTON')
c      call import_parameters(-1,parameters,'H_INT','DALTON_SPECIAL')
      call set_rule('VXX(1)LIST',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('VX(1)LIST',ttype_opme,.false.,tgt_info)
      call set_dependency('VX(1)LIST','VX(1)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'VX(1)LIST'
      labels(2) = 'VX(1)'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('VX(1)LIST',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'VX(1)LIST'
      call import_parameters(-1,parameters,'ZDIPLEN','DALTON')
c      call import_parameters(-1,parameters,'H_INT','DALTON_SPECIAL')
      call set_rule('VX(1)LIST',ttype_opme,IMPORT,
     &              labels,1,1,
     &              parameters,1,tgt_info)

      call add_target('DEF-RDGB(VZ)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-RDGB(VZ)','RDGB(V)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RDGB(VZ)'
      labels(2) = 'RDGB(V)'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-RDGB(VZ)',ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)

      call add_target('DEF-RBRV(VZ)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-RDGB(VZ)','RBRV(V)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'RBRV(VZ)'
      labels(2) = 'RBRV(V)'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-RBRV(VZ)',ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)

      call add_target('DEF-B(VZ)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF-B(VZ)','BV',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'B(VZ)'
      labels(2) = 'BV'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,0,0,.false.)
      call set_rule('DEF-B(VZ)',ttype_opme,DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)

c      call add_target('V(1)LISTX',ttype_opme,.false.,tgt_info)
c      call set_dependency('V(1)LISTX','V(1)',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'V(1)LISTX'
c      labels(2) = 'V(1)'
c      call me_list_parameters(-1,parameters,
c     &     msc,0,2,0,0,.false.)
c      call set_rule('V(1)LISTX',ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'V(1)LISTX'
c      call import_parameters(-1,parameters,'XDIPLEN','DALTON')
c      call set_rule('V(1)LISTX',ttype_opme,IMPORT,
c     &              labels,1,1,
c     &              parameters,1,tgt_info)
      


c
c      ! e.g.: 
c      call add_target('MY_MEL1',ttype_opme,.false.,tgt_info)
c      call set_dependency('MY_MEL1','MY_OP1',tgt_info)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'MY_MEL1'
c      labels(2) = 'MY_OP1'
c      call me_list_parameters(-1,parameters,
c     &     msc,0,1,0,0,.false.)
c      call set_rule('MY_MEL1',ttype_opme,DEF_ME_LIST,
c     &              labels,2,1,
c     &              parameters,1,tgt_info)
c      
*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      call add_target('EVAL-B(VZ)',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL-B(VZ)','B(VZ)_OPT',tgt_info)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'B(VZ)_OPT'
      call set_rule('EVAL-B(VZ)',ttype_opme,EVAL,
     &     labels,1,0,
     &     parameters,0,tgt_info)


      call add_target('MY_TARGET',ttype_gen,.true.,tgt_info)
c      call set_dependency('MY_TARGET','FOCK',tgt_info)
c      call set_dependency('MY_TARGET','PHI',tgt_info)
c      call set_dependency('MY_TARGET','EMP2',tgt_info)
c      call set_dependency('MY_TARGET','LMP2',tgt_info)
c      call set_dependency('MY_TARGET','T2',tgt_info)
c      call set_dependency('MY_TARGET','O2',tgt_info)
c      call set_dependency('MY_TARGET','MP2_FORM2',tgt_info)
c      call set_dependency('MY_TARGET','V(1)LIST',tgt_info)
c      call set_dependency('MY_TARGET','VX(1)LIST',tgt_info)
      call set_dependency('MY_TARGET','BV_formal',tgt_info)
      call set_dependency('MY_TARGET','EVAL-B(VZ)',tgt_info)
c      call set_dependency('MY_TARGET','BV_CABS',tgt_info)
c      call set_dependency('MY_TARGET','RDGB(V)_FRM',tgt_info)
c      call set_dependency('MY_TARGET','RBRV(V)_FRM',tgt_info)
c      call set_dependency('MY_TARGET','CC_FORM1',tgt_info)
c      call set_dependency('MY_TARGET','CC_FORM2',tgt_info)
c      call set_dependency('MY_TARGET','MP2_RES',tgt_info)
c      call set_dependency('MY_TARGET','MP2_EN',tgt_info)
c      call set_dependency('MY_TARGET','MP2_OPT',tgt_info)
c      call set_dependency('MY_TARGET','DEF_ME_LMP2',tgt_info)
c      call set_dependency('MY_TARGET','DEF_ME_T2',tgt_info)
c      call set_dependency('MY_TARGET','DIAG',tgt_info)

c      call solve_parameters(-1,parameters,2, 1,1,'DIA')
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'ME_T2'
c      labels(2) = 'ME_O2'
c      labels(3) = 'DIAG'
c      labels(4) = 'ME_LMP2'
c      labels(5) = 'MP2_OPT'
c      call set_rule('MY_TARGET',ttype_opme,SOLVENLEQ,
c     &     labels,5,2,
c     &     parameters,2,tgt_info)

      return
      end
