*----------------------------------------------------------------------*
      subroutine set_ic_mrci_targets(tgt_info,orb_info,
     &                               excrestr,maxh,maxp,use_met)
*----------------------------------------------------------------------*
*     contains targets for internally contracted MRCI
*
*     matthias, 2009/2010
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
      integer, intent(in) ::
     &     maxh, maxp, excrestr(0:maxh,0:maxp,1:2)
      logical, intent(in) ::
     &     use_met

      integer ::
     &     ndef, occ_def(ngastp,2,124),!60),
     &     ncat, msc, ip, ih,
     &     nlabels, nroots, gno, idef, iexc, jexc, prc_type
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      logical ::
     &     skip

      if (iprlvl.gt.0) write(luout,*) 'setting icMRCI targets'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1
      if (orb_info%ims.ne.0) msc = 0

      call get_argument_value('method.MR','GNO',
     &     ival=gno)
      call get_argument_value('method.MR','prc_type',
     &     ival=prc_type)

      call get_argument_value('method.MRCI','nroots',
     &     ival=nroots)
      skip = (is_keyword_set('calculate.skip_E').gt.0)

      if (ntest.ge.100) then
        print *,'nroots  = ',nroots
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving excitation operator C (for icMRCI)
      call add_target('C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('C',ttype_op,DEF_OP_FROM_OCC,
     &              'C',1,1,
     &              parameters,2,tgt_info)

      ! "transformed" C needed for generation of formula for effective H
      call add_target('Ctr',ttype_op,.false.,tgt_info)
      call set_dependency('Ctr','C',tgt_info)
      call cloneop_parameters(-1,parameters,'C',.false.)
      call set_rule('Ctr',ttype_op,CLONE_OP,'Ctr',1,1,
     &              parameters,1,tgt_info)
      ! transposed "transformed" C needed for generation of effective H
      call add_target('Cdag',ttype_op,.false.,tgt_info)
      call set_dependency('Cdag','C',tgt_info)
      call cloneop_parameters(-1,parameters,'C',.true.)
      call set_rule('Cdag',ttype_op,CLONE_OP,'Cdag',1,1,
     &              parameters,1,tgt_info)

      if (gno.eq.1) then
        ! subset of C for non-redundant valence-only metric
        call add_target('c',ttype_op,.false.,tgt_info)
        occ_def = 0
        ndef = 0
        do ip = 0, maxp
          do ih = 0, maxh
            do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
              ! exclude blocks with one or more hole-part. excitations
              if (min(ih,ip).ge.1) cycle
              if (max(ih,ip).eq.0.and.iexc.eq.0) cycle
              ndef = ndef + 1
              occ_def(IHOLE,2,ndef) = ih
              occ_def(IPART,1,ndef) = ip
              occ_def(IVALE,1,ndef) = iexc - ip
              occ_def(IVALE,2,ndef) = iexc - ih
            end do
          end do
        end do
        call op_from_occ_parameters(-1,parameters,2,
     &                occ_def,ndef,1,(/0,0/),ndef)
        call set_rule('c',ttype_op,DEF_OP_FROM_OCC,
     &                'c',1,1,
     &                parameters,2,tgt_info)
      end if

      ! define Jacobian
      call add_target('A_C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = iexc - ip
            occ_def(IVALE,2,ndef*2-1) = iexc - ih
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('A_C',ttype_op,DEF_OP_FROM_OCC,
     &              'A_C',1,1,
     &              parameters,2,tgt_info)

      ! Metric times icCI vector product
      call add_target('SC',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('SC','A_C',tgt_info)
      call cloneop_parameters(-1,parameters,'A_C',.false.)
      call set_rule('SC',ttype_op,CLONE_OP,'SC',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner
      call add_target(trim(op_dia)//'_C',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(trim(op_dia)//'_C','C',tgt_info)
      call cloneop_parameters(-1,parameters,'C',.false.)
      call set_rule(trim(op_dia)//'_C',ttype_op,CLONE_OP,
     &              trim(op_dia)//'_C',1,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! multireference energy expression
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'F_E(MR)'
c      labels(2) = 'E(MR)'
c      labels(3) = 'C0^+'
c      labels(4) = 'C^+'
c      labels(5) = op_ham
c      labels(6) = 'C'
c      labels(7) = 'C0'
c      call add_target('F_E(MR)',ttype_frm,.false.,tgt_info)
c      call set_dependency('F_E(MR)','E(MR)',tgt_info)
c      call set_dependency('F_E(MR)',op_ham,tgt_info)
c      call set_dependency('F_E(MR)','C0',tgt_info)
c      call set_dependency('F_E(MR)','C',tgt_info)
c      call expand_parameters(-1,
c     &     parameters,3,
c     &     'multireference energy expression',5,
c     &     (/2,3,4,5,6/),
c     &     (/-1,-1,-1,-1,-1/),
c     &     (/-1,-1,-1,-1,-1/),
c     &     0,0,
c     &     0,0,
c     &     0,0)
c      call set_rule('F_E(MR)',ttype_frm,EXPAND_OP_PRODUCT,
c     &              labels,7,1,
c     &              parameters,3,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_E(MR)',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)
      call add_target2('F_preE(MR)',.false.,tgt_info)
      call set_dependency('F_preE(MR)','E(MR)',tgt_info)
      call set_dependency('F_preE(MR)','C',tgt_info)
      call set_dependency('F_preE(MR)','Cdag',tgt_info)
      call set_dependency('F_preE(MR)',op_ham,tgt_info)
      if (gno.eq.0) then
        call set_dependency('F_preE(MR)','C0',tgt_info)
        call set_rule2('F_preE(MR)',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'C0^+','Cdag',op_ham(1:4),'C   ','C0  '/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/2,3,4,5,6/))
      else if (gno.eq.1) then
        call set_dependency('F_preE(MR)','DENS',tgt_info)
        call set_rule2('F_preE(MR)',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'Cdag',op_ham(1:4),'C   '/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/2,3,4/))
        call set_rule2('F_preE(MR)',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'DENS','Cdag',op_ham(1:4),'C   ','DENS'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/2,3,4,5,2/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &       tgt_info,val_int=(/1/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &       val_int=(/1,5/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        ! b) delete terms which are anyway forbidden by contraction rules
        call set_rule2('F_preE(MR)',SELECT_SPECIAL,tgt_info)
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'OPERATORS',3,tgt_info,
     &       val_label=(/op_ham(1:4),'C0  ','DENS'/)) ! C0 is dummy
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC')
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='no_con')
        ! c) expand reduced densities in terms of cumulants
        call set_dependency('F_preE(MR)','F_DENS',tgt_info)
        call set_rule2('F_preE(MR)',EXPAND,tgt_info)
        call set_arg('F_preE(MR)',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_DENS'/))
        ! d) select only terms allowed according to contraction rules
        call set_rule2('F_preE(MR)',SELECT_SPECIAL,tgt_info)
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'OPERATORS',3,tgt_info,
     &       val_label=(/op_ham(1:4),'C0  ','CUM '/)) ! C0 is dummy
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC')
        call set_arg('F_preE(MR)',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='no_con')
        ! e) factor out hole matrices by
        !    a) inserting valence unity 1v 
        !    b) replacing 1v with 1
        !    c) factoring out HOLE
        ! e1) Cdag - H contractions
        call set_dependency('F_preE(MR)','1v',tgt_info)
        call set_rule2('F_preE(MR)',INSERT,tgt_info)
        call set_arg('F_preE(MR)',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1v'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/'Cdag',op_ham(1:4)/))
        call set_dependency('F_preE(MR)','1',tgt_info)
        call set_rule2('F_preE(MR)',REPLACE,tgt_info)
        call set_arg('F_preE(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1v','1 '/))
        call set_dependency('F_preE(MR)','F_HOLE',tgt_info)
        call set_rule2('F_preE(MR)',FACTOR_OUT,tgt_info)
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_HOLE'/))
        ! e2) Cdag - C contractions
        call set_rule2('F_preE(MR)',INSERT,tgt_info)
        call set_arg('F_preE(MR)',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1v'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/'Cdag','C   '/))
        call set_rule2('F_preE(MR)',REPLACE,tgt_info)
        call set_arg('F_preE(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1v','1 '/))
        call set_rule2('F_preE(MR)',FACTOR_OUT,tgt_info)
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_HOLE'/))
        ! e3) H - C contractions
        call set_rule2('F_preE(MR)',INSERT,tgt_info)
        call set_arg('F_preE(MR)',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1v'/))
        call set_arg('F_preE(MR)',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/op_ham,'C'/))
        call set_rule2('F_preE(MR)',REPLACE,tgt_info)
        call set_arg('F_preE(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1v','1 '/))
        call set_rule2('F_preE(MR)',FACTOR_OUT,tgt_info)
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_HOLE'/))
      end if

      ! replace Cdag by C^+
      call add_target2('F_E(MR)',.false.,tgt_info)
      call set_dependency('F_E(MR)','F_preE(MR)',tgt_info)
      call set_dependency('F_E(MR)','F_E(MR)_diag',tgt_info)
      call set_rule2('F_E(MR)',REPLACE,tgt_info)
      call set_arg('F_E(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MR)'/))
      call set_arg('F_E(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MR)'/))
      call set_arg('F_E(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'Cdag','C^+ '/))
      call set_rule2('F_E(MR)',PRINT_FORMULA,tgt_info)
      call set_arg('F_E(MR)',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_E(MR)'/))

      ! transformed C needed for generation of formula for effective H
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_C'
      labels(2) = 'C'
      labels(3) = 'C'
      labels(4) = 'Dtr'
      labels(5) = 'Ctr'
      labels(6) = 'Dtr'
      labels(7) = 'C'
      call add_target('F_C',ttype_frm,.false.,tgt_info)
      call set_dependency('F_C','Ctr',tgt_info)
      call set_dependency('F_C','C',tgt_info)
      call set_dependency('F_C','Dtr',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'reference energy expression',5,
     &     (/1,2,3,2,1/),
     &     (/-1,-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1,-1/),
     &     0,0,
     &     (/3,5,2,4,2,5,1,4/),4,
     &     0,0)
      call set_rule('F_C',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,7,1,
     &              parameters,3,tgt_info)
      call set_rule2('F_C',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_C',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_C'/))
      call set_arg('F_C',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'C'/))
      call set_arg('F_C',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'C  ','Ctr','C  '/))
      call set_arg('F_C',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_C',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! delete terms with active open lines from the Ctr operator
      labels(2:20)(1:len_target_name) = ' '
      labels(2) = 'F_C'
      labels(3) = 'C'
      labels(4) = 'Ctr'
      call form_parameters(-1,
     &       parameters,2,'C formula',3,'no_ext')
      call set_rule('F_C',ttype_frm,SELECT_LINE,
     &              labels,4,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_C',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! transposed transformed C
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_Cdag'
      labels(2) = 'Cdag'
      labels(3) = 'Cdag'
      labels(4) = 'Dtr^+'
      labels(5) = 'Ctr^+'
      labels(6) = 'Dtr^+'
      labels(7) = 'Cdag'
      call add_target('F_Cdag',ttype_frm,.false.,tgt_info)
      call set_dependency('F_Cdag','Cdag',tgt_info)
      call set_dependency('F_Cdag','Ctr',tgt_info)
      call set_dependency('F_Cdag','Dtr',tgt_info)
      call expand_parameters(-1,
     &     parameters,3,
     &     'reference energy expression',5,
     &     (/1,2,3,2,1/),
     &     (/-1,-1,-1,-1,-1/),
     &     (/-1,-1,-1,-1,-1/),
     &     0,0,
     &     (/1,3,2,4,2,5,1,4/),4,
     &     0,0)
      call set_rule('F_Cdag',ttype_frm,EXPAND_OP_PRODUCT,
     &              labels,7,1,
     &              parameters,3,tgt_info)
      call set_rule2('F_Cdag',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_Cdag',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_Cdag'/))
      call set_arg('F_Cdag',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'Cdag'/))
      call set_arg('F_Cdag',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'Cdag ','Ctr^+','Cdag '/))
      call set_arg('F_Cdag',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_Cdag',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! delete terms with active open lines from the C operator
      labels(2:20)(1:len_target_name) = ' '
      labels(2) = 'F_Cdag'
      labels(3) = 'Cdag'
      labels(4) = 'Ctr^+'
      call form_parameters(-1,
     &       parameters,2,'Cdag formula',3,'no_ext')
      call set_rule('F_Cdag',ttype_frm,SELECT_LINE,
     &              labels,4,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_Cdag',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! only diagonal terms of multireference energy expression (for PREC)
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'F_E(MR)_diag'
c      labels(2) = 'E(MR)'
c      labels(3) = 'C0^+'
c      labels(4) = 'Cdag'
c      labels(5) = op_ham
c      labels(6) = 'C'
c      labels(7) = 'C0'
      call add_target2('F_E(MR)_diag',.false.,tgt_info)
      call set_dependency('F_E(MR)_diag','E(MR)',tgt_info)
c      call set_dependency('F_E(MR)_diag',op_ham,tgt_info)
c      call set_dependency('F_E(MR)_diag','C0',tgt_info)
      call set_dependency('F_E(MR)_diag','F_preE(MR)',tgt_info)
      call set_dependency('F_E(MR)_diag','F_C',tgt_info)
      call set_dependency('F_E(MR)_diag','F_Cdag',tgt_info)
c      call expand_parameters(-1,
c     &     parameters,3,
c     &     'multireference energy expression',5,
c     &     (/2,3,4,5,6/),
c     &     (/-1,-1,-1,-1,-1/),
c     &     (/-1,-1,-1,-1,-1/),
c     &     0,0,
c     &     0,0,
c     &     0,0)
c      call set_rule('F_E(MR)_diag',ttype_frm,EXPAND_OP_PRODUCT,
c     &              labels,7,1,
c     &              parameters,3,tgt_info)
      ! insert (particle/hole) unit operator to allow for differentiation
      if (prc_type.lt.3) then
        call set_dependency('F_E(MR)_diag','1ph',tgt_info)
        call set_rule2('F_E(MR)_diag',INSERT,tgt_info)
        call set_arg('F_E(MR)_diag',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)_diag'/))
        call set_arg('F_E(MR)_diag',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_E(MR)_diag',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_E(MR)_diag',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1ph'/))
        call set_arg('F_E(MR)_diag',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/'Cdag','C   '/))
        ! replace 1ph by 1
        call set_dependency('F_E(MR)_diag','1',tgt_info)
        call set_rule2('F_E(MR)_diag',REPLACE,tgt_info)
        call set_arg('F_E(MR)_diag',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)_diag'/))
        call set_arg('F_E(MR)_diag',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)_diag'/))
        call set_arg('F_E(MR)_diag',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1ph','1  '/))
      else
        call set_rule2('F_E(MR)_diag',INVARIANT,tgt_info)
        call set_arg('F_E(MR)_diag',INVARIANT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)_diag'/))
        call set_arg('F_E(MR)_diag',INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_E(MR)_diag',INVARIANT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_E(MR)_diag',INVARIANT,'OPERATORS',0,tgt_info,
     &       val_label=(/''/))
        call set_arg('F_E(MR)_diag',INVARIANT,'TITLE',1,tgt_info,
     &       val_str='---')
      end if
      ! expand C --> Dtr Ctr
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_E(MR)_diag'
      labels(2) = 'F_E(MR)_diag'
      labels(3) = 'F_C'
      labels(4) = 'F_Cdag'
      call form_parameters(-1,
     &     parameters,2,'full precursor for orthog. eff. H',2,'---')
      call set_rule('F_E(MR)_diag',ttype_frm,EXPAND,
     &              labels,4,1,
     &              parameters,2,tgt_info)
      labels(3:20)(1:len_target_name) = ' '
      labels(3) = 'E(MR)'
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_E(MR)_diag',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! effective hamiltonian times icCI coefficient vector
      call add_target2('F_A_C',.false.,tgt_info)
      call set_dependency('F_A_C','F_E(MR)',tgt_info)
      call set_dependency('F_A_C','A_C',tgt_info)
      call set_dependency('F_A_C','C',tgt_info)
      call set_rule2('F_A_C',DERIVATIVE,tgt_info)
      call set_arg('F_A_C',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_A_C'/))
      call set_arg('F_A_C',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MR)'/))
      call set_arg('F_A_C',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C'/))
      call set_arg('F_A_C',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'C^+'/))
c      call set_rule2('F_A_C',PRINT_FORMULA,tgt_info)
c      call set_arg('F_A_C',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_A_C'/))

      ! diagonal terms of effective hamiltonian times icCI coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'F_A_C_diag'
      labels(2) = 'F_E(MR)_diag'
      labels(3) = 'OMGtr'
      labels(4) = 'Ctr^+'
      labels(5) = ' '
      call add_target('F_A_C_diag',ttype_frm,.false.,tgt_info)
      call set_dependency('F_A_C_diag','F_E(MR)_diag',tgt_info)
      call set_dependency('F_A_C_diag','OMGtr',tgt_info)
      call set_dependency('F_A_C_diag','Ctr',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'diag of Jacobian',1,'---')
      call set_rule('F_A_C_diag',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_A_C_diag',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! diagonal terms of Jacobian
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'F_A_diag'
      labels(2) = 'F_A_C_diag'
      labels(3) = 'A'
      labels(4) = 'Ctr'
      labels(5) = ' '
      call add_target('F_A_diag',ttype_frm,.false.,tgt_info)
      call set_dependency('F_A_diag','F_A_C_diag',tgt_info)
      call set_dependency('F_A_diag','A',tgt_info)
      call set_dependency('F_A_diag','Ctr',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'diag of Hessian',1,'---')
      call set_rule('F_A_diag',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_A_diag',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! factor out reduced density matrices
c      if (use_met) then
c        call add_target('F_NORM_fact',ttype_frm,.false.,tgt_info)
c        call set_dependency('F_NORM_fact','F_NORM',tgt_info)
c        labels(1:10)(1:len_target_name) = ' '
c        labels(1) = 'F_NORM_fact'
c        labels(2) = 'F_NORM'
c        labels(3) = 'F_D'
c        call set_dependency('F_NORM_fact','F_D',tgt_info)
c        call form_parameters(-1,
c     &       parameters,2,'norm using density matrices',1,'---')
c        call set_rule('F_NORM_fact',ttype_frm,FACTOR_OUT,
c     &                labels,3,1,
c     &                parameters,2,tgt_info)
c        ! output needed for test suite (factoring out is checked)
c        call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c        call set_rule('F_NORM_fact',ttype_frm,PRINT_FORMULA,
c     &                  labels,2,1,parameters,2,tgt_info)
      ! a) set up norm expression
      call add_target2('F_NORM_fact',.false.,tgt_info)
      call set_dependency('F_NORM_fact','NORM',tgt_info)
      call set_dependency('F_NORM_fact','C',tgt_info)
      call set_dependency('F_NORM_fact','D',tgt_info)
      call set_rule2('F_NORM_fact',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'D  ','C^+','D  ','C  ','D  '/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/2,3,2,4,2/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/5/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'AVOID',10,tgt_info,
     &     val_int=(/1,3,1,4,1,5,2,5,3,5/))
      call set_rule2('F_NORM_fact',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &     tgt_info,
     &     val_label=(/'C^+','C  '/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &     val_int=(/2,3/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! no active lines between C^+ and C
      call set_rule2('F_NORM_fact',SELECT_LINE,tgt_info)
      call set_arg('F_NORM_fact',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'OP_INCL',2,tgt_info,
     &     val_label=(/'C^+','C  '/))
      call set_arg('F_NORM_fact',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_NORM_fact',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='delete')
      call set_rule2('F_NORM_fact',PRINT_FORMULA,tgt_info)
      call set_arg('F_NORM_fact',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
c      end if

      ! metric times icCI coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'F_SC'
      labels(2) = 'F_NORM'
      labels(3) = 'SC'
      labels(4) = 'C^+'
      labels(5) = ' '
      call add_target('F_SC',ttype_frm,.false.,tgt_info)
      if (use_met) then
        call set_dependency('F_SC','F_NORM_fact',tgt_info)
        labels(2) = 'F_NORM_fact'
      else
        call set_dependency('F_SC','F_NORM',tgt_info)
      end if
      call set_dependency('F_SC','SC',tgt_info)
      call set_dependency('F_SC','C',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'metric times right hand vector',1,'---')
      call set_rule('F_SC',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_SC',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! effective hamiltonian times icCI coefficient vector
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_A_C'
      labels(2) = 'F_A_C'
      labels(3) = 'F_SC'
      ncat = 2
      call add_target('FOPT_A_C',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_A_C','F_A_C',tgt_info)
      call set_dependency('FOPT_A_C',mel_ham,tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_C',tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_A_C',tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_A_C','F_SC',tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_SC',tgt_info)
      call set_dependency('FOPT_A_C','DEF_ME_1',tgt_info)
      if (use_met) then
        call set_dependency('FOPT_A_C','DEF_ME_D',tgt_info)
        labels(4) = 'F_C'
        ncat = 3
        call set_dependency('FOPT_A_C','F_C',tgt_info)
        call set_dependency('FOPT_A_C','DEF_ME_Dtr',tgt_info)
        call set_dependency('FOPT_A_C','DEF_ME_Ctr',tgt_info)
      end if
      call opt_parameters(-1,parameters,ncat,0)
      call set_rule('FOPT_A_C',ttype_frm,OPTIMIZE,
     &              labels,ncat+1,1,
     &              parameters,1,tgt_info)

      ! multireference energy expression
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_E(MR)'
      labels(2) = 'F_E(MR)'
      call add_target('FOPT_E(MR)',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_E(MR)','F_E(MR)',tgt_info)
      call set_dependency('FOPT_E(MR)',mel_ham,tgt_info)
      call set_dependency('FOPT_E(MR)','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_E(MR)','DEF_ME_C',tgt_info)
      call set_dependency('FOPT_E(MR)','DEF_ME_E(MR)',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_E(MR)',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! diagonal terms of Hessian
      labels(1:20)(1:len_target_name)= ' '
      call add_target('FOPT_A_diag',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_A_diag','F_A_diag',tgt_info)
      call set_dependency('FOPT_A_diag',mel_ham,tgt_info)
      call set_dependency('FOPT_A_diag','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_A_diag','DEF_ME_A',tgt_info)
      call set_dependency('FOPT_A_diag','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_A_diag','DEF_ME_Dtr',tgt_info)
      labels(1) = 'FOPT_A_diag'
      labels(2) = 'F_A_diag'
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_A_diag',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! norm
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_NORM'
      labels(2) = 'F_NORM_fact'
      call add_target('FOPT_NORM',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_NORM','F_NORM_fact',tgt_info)
c      call set_dependency('FOPT_NORM','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_NORM','DEF_ME_DENS',tgt_info)
      call set_dependency('FOPT_NORM','DEF_ME_C',tgt_info)
      call set_dependency('FOPT_NORM','DEF_ME_NORM',tgt_info)
      call set_dependency('FOPT_NORM','DEF_ME_1',tgt_info)
c      if (gno.eq.1) then
c        call set_dependency('FOPT_NORM','DEF_ME_HOLE',tgt_info)
c        call set_dependency('FOPT_NORM','DEF_ME_CUM',tgt_info)
c      end if
      if (use_met.or..true.)
     &        call set_dependency('FOPT_NORM','DEF_ME_D',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_NORM',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME_C
      call add_target2('DEF_ME_C',.false.,tgt_info)
      call set_dependency('DEF_ME_C','C',tgt_info)
      call set_rule2('DEF_ME_C',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_C',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_C'/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'C'/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_C',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/nroots/))

      ! ME_Ctr
      call add_target2('DEF_ME_Ctr',.false.,tgt_info)
      call set_dependency('DEF_ME_Ctr','Ctr',tgt_info)
      call set_rule2('DEF_ME_Ctr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ctr'/))
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Ctr'/))
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/nroots/))

      ! ME_A_C
      call add_target('DEF_ME_A_C',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_A_C','A_C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_A_C'
      labels(2) = 'A_C'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_A_C',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_SC
      call add_target('DEF_ME_SC',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_SC','SC',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_SC'
      labels(2) = 'SC'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_SC',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner
      call me_list_label(dia_label,mel_dia,1,
     &     0,0,0,.false.)
      call add_target(trim(dia_label)//'C',ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label)//'C','EVAL_FREF',tgt_info)
      call set_dependency(trim(dia_label)//'C',
     &                    trim(op_dia)//'_'//'C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C'
      labels(2) = trim(op_dia)//'_'//'C'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule(trim(dia_label)//'C',ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      ! use effective Fock op. (needed only for pure inactive exc.)
      labels(1) = trim(dia_label)//'C'
      labels(2) = 'ME_FREF'
      if (prc_type.ne.3) then
        call set_rule(trim(dia_label)//'C',ttype_opme,
     &                PRECONDITIONER,
     &                labels,2,1,
     &                parameters,3,tgt_info)
      else ! no scalar part (taken care of by active part)
        call set_rule(trim(dia_label)//'C',ttype_opme,
     &                PRECONDITIONER,
     &                labels,2,1,
     &                parameters,0,tgt_info)
      end if

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Evaluate norm
      call add_target('EVAL_NORM',ttype_gen,.not.skip,tgt_info)
      call set_dependency('EVAL_NORM','FOPT_NORM',tgt_info)
      call set_dependency('EVAL_NORM','SOLVE_ICCI',tgt_info)
      call set_rule('EVAL_NORM',ttype_opme,EVAL,
     &     'FOPT_NORM',1,0,
     &     parameters,0,tgt_info)

      ! Evaluate diagonal elements of Jacobian
      call add_target('EVAL_A_diag',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_A_diag','FOPT_A_diag',tgt_info)
      call set_dependency('EVAL_A_diag','EVAL_REF_S(S+1)',tgt_info)
      if (gno.eq.1)
     &   call set_dependency('EVAL_A_diag','H_GNO',tgt_info)
      call set_rule('EVAL_A_diag',ttype_opme,EVAL,
     &     'FOPT_A_diag',1,0,
     &     parameters,0,tgt_info)
      ! put diagonal elements to preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      call set_dependency('EVAL_A_diag',trim(dia_label)//'C',tgt_info)
      call set_rule2('EVAL_A_diag',EXTRACT_DIAG,tgt_info)
      call set_arg('EVAL_A_diag',EXTRACT_DIAG,'LIST_RES',1,tgt_info,
     &             val_label=(/trim(dia_label)//'C'/))
      call set_arg('EVAL_A_diag',EXTRACT_DIAG,'LIST_IN',1,tgt_info,
     &             val_label=(/'ME_A'/))
      if (prc_type.ge.3)
     &  call set_arg('EVAL_A_diag',EXTRACT_DIAG,'EXTEND',1,tgt_info,
     &               val_log=(/.true./))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (b) :',0,'LIST')
c      call set_rule('EVAL_A_diag',ttype_opme,PRINT_MEL,
c     &     trim(dia_label)//'C',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! SOLVE icCI eigenvalue equation
      call add_target('SOLVE_ICCI',ttype_gen,.false.,tgt_info)
      call set_dependency('SOLVE_ICCI','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('SOLVE_ICCI','FOPT_A_C',tgt_info)
      if (gno.eq.1)
     &   call set_dependency('SOLVE_ICCI','H_GNO',tgt_info)
      call me_list_label(dia_label,mel_dia,1,
     &     0,0,0,.false.)
      if (use_met)
     &       call set_dependency('SOLVE_ICCI','EVAL_A_diag',tgt_info)
      call set_dependency('SOLVE_ICCI',trim(dia_label)//'C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      if (use_met) then
        labels(1) = 'ME_D'
        labels(2) = 'D'
        call set_rule('SOLVE_ICCI',ttype_opme,ASSIGN_ME2OP,
     &       labels,2,1,
     &       parameters,0,tgt_info)
      end if
      call solve_parameters(-1,parameters,2,1,nroots,'TRF')
      labels(1) = 'ME_C'
      labels(2) = trim(dia_label)//'C'
      labels(3) = 'A_C'
      labels(4) = 'SC'
      labels(5) = 'FOPT_A_C'
      nlabels = 5
      if (use_met) then
        call set_dependency('SOLVE_ICCI','EVAL_D',tgt_info)
        call set_dependency('SOLVE_ICCI','DEF_ME_Dtrdag',tgt_info)
        labels(6) = 'ME_Ctr'
        labels(7) = 'ME_Dtr'
        labels(8) = 'ME_Dtrdag'
        nlabels = 8
      end if
      call set_rule('SOLVE_ICCI',ttype_opme,SOLVEEVP,
     &     labels,nlabels,1,
     &     parameters,2,tgt_info)

      ! Evaluate multireference energy
      call add_target('EVAL_E(MR)',ttype_gen,.not.skip,tgt_info)
      call set_dependency('EVAL_E(MR)','SOLVE_ICCI',tgt_info)
      call set_dependency('EVAL_E(MR)','FOPT_E(MR)',tgt_info)
      if (gno.eq.1)
     &   call set_dependency('EVAL_E(MR)','H_GNO',tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'icCI coefficients :',0,'LIST')
c      call set_rule('EVAL_E(MR)',ttype_opme,PRINT_MEL,
c     &     'ME_C',1,0,
c     &     parameters,2,tgt_info)
      call set_rule('EVAL_E(MR)',ttype_opme,EVAL,
     &     'FOPT_E(MR)',1,0,
     &     parameters,0,tgt_info)

      return
      end
