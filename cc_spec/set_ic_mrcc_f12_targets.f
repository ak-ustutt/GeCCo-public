*----------------------------------------------------------------------*
      subroutine set_ic_mrcc_f12_targets(tgt_info,orb_info,
     &                               excrestr,maxh,maxp)
*----------------------------------------------------------------------*
*     set targets for internally contracted MRCC_F12
*
*     matthias 2010, wenlan 2012
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

      integer ::
     &     ndef, occ_def(ngastp,2,124),!60),
     &     icnt, maxexc, igastp,
     &     msc, ip, ih, ivv, iv, ivv2, jvv,
     &     maxcom, maxcom_en, maxcom_h1bar, h1bar_maxp,
     &     n_t_cls, i_cls,
     &     n_tred_cls, len_form, optref, idef, ciroot,
     &     version(60), ivers, stndT(2,60), stndD(2,60), nsupT, nsupD,
     &     G_level, iexc, jexc, maxtt, iblk, jblk, kblk, prc_type,
     &     tred, nremblk, remblk(60), igasreo(3), ngas, lblk, ntrunc,
     &     tfix, maxit, t1ord, maxcum, cum_appr_mode, update_prc
      logical ::
     &     skip, preopt, project, first, Op_eqs,
     &     h1bar, htt, svdonly, fact_tt, ex_t3red, trunc, l_exist,
     &     oldref, solve, notrunc
      character(len_target_name) ::
     &     dia_label, dia_label2,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     op_ht*3, f_ht*5, op_ht0to*6, f_ht0to*8, form_str*50,
     &     def_ht*10
      real(8) ::
     &     x_ansatz, prc_shift

      if (iprlvl.gt.0) write(luout,*) 'setting icMRCC_F12 targets'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1
      if (orb_info%ims.ne.0) msc = 0

c      ! get some keywords
c      call get_argument_value('method.MR','maxexc',
c     &     ival=maxexc)
c      call get_argument_value('method.MR','ciroot',
c     &     ival=ciroot)
c      call get_argument_value('method.MR','prc_type',
c     &     ival=prc_type)
c      call get_argument_value('method.MR','prc_shift',
c     &     xval=prc_shift)
      call get_argument_value('method.MR','svdonly',
     &     lval=svdonly)
      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)
      call get_argument_value('calculate.solve.non_linear','update_prc',
     &     ival=update_prc)
c      call get_argument_value('calculate.solve.non_linear','preopt',
c     &     lval=preopt)
c      call get_argument_value('method.MR','project',
c     &     lval=project)
      call get_argument_value('method.MRCC','maxcom_res',
     &     ival=maxcom)
      call get_argument_value('method.MRCC','maxcom_en',
     &     ival=maxcom_en)
c      call get_argument_value('method.MRCC','maxcom_h1bar',
c     &     ival=maxcom_h1bar)
c      call get_argument_value('method.MRCC','h1bar_maxp',
c     &     ival=h1bar_maxp)
c      if (h1bar_maxp.lt.0.and.maxexc.le.2) h1bar_maxp = 2
c      if (h1bar_maxp.lt.0.and.maxexc.gt.2) h1bar_maxp = 3
      call get_argument_value('method.MRCC','H1bar',
     &     lval=h1bar)
c      call get_argument_value('method.MRCC','x_ansatz',
c     &     xval=x_ansatz)
c      call get_argument_value('method.MRCC','trunc_order',
c     &     ival=ntrunc)
c      call get_argument_value('method.MRCC','Tfix',
c     &     ival=tfix)
c      call get_argument_value('method.MRCC','T1ord',
c     &     ival=t1ord)
c      call get_argument_value('method.MR','oldref',
c     &     lval=oldref)
c      call get_argument_value('method.MR','maxcum',
c     &     ival=maxcum)
c      call get_argument_value('method.MR','cum_appr_mode',
c     &     ival=cum_appr_mode)
      call get_argument_value('method.R12','notrunc',lval=notrunc)
c      call get_argument_value('calculate.solve','maxiter',
c     &     ival=maxit)
c      if (is_argument_set('calculate.solve.non_linear','maxiter').gt.0)
c     &     call get_argument_value('calculate.solve.non_linear',
c     &     'maxiter',ival=maxit)
c      trunc = ntrunc.ge.0
c      solve = .not.svdonly.and.(tfix.eq.0.or.maxit.gt.1)
      skip = (is_keyword_set('calculate.skip_E').gt.0)
      solve = .not.svdonly.and..not.skip
      if (h1bar) call quit(1,'set_ic_mrcc_f12_targets',
     &                     'H1bar not available yet')
      if (optref.ne.0.and.optref.ne.-3)
     &   call quit(1,'set_ic_mrcc_f12_targets','use optref=0 or -3')
      
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving excitation operator S (for icMRCC)
      call add_target2('S',.false.,tgt_info)
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
      do ih = 0, maxh
        do ip = 1, 0, -1
          ndef = ndef + 1
          occ_def(IHOLE,2,ndef) = ih
          occ_def(IPART,1,ndef) = ip
          occ_def(IVALE,2,ndef) = 2 - ih
          occ_def(IEXTR,1,ndef) = 2 - ip
        end do
      end do
      call set_rule2('S',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('S',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'S'/))
      call set_arg('S',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &     val_int=(/ndef/))
      call set_arg('S',DEF_OP_FROM_OCC,'OCC',ndef,tgt_info,
     &     val_occ=occ_def(1:ngastp,1:2,1:ndef))

      ! Sbar operator
      call add_target2('Sbar',.false.,tgt_info)
      call set_dependency('Sbar','S',tgt_info)
      call set_rule2('Sbar',CLONE_OP,tgt_info)
      call set_arg('Sbar',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Sbar'/))
      call set_arg('Sbar',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'S'/))
      call set_arg('Sbar',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.true./))

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! multireference CC lagrangian
      ! a) set up
      call add_target2('F_MRCC_F12_LAG',.false.,tgt_info)
      call set_dependency('F_MRCC_F12_LAG','C0',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','S',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','Sbar',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','H',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','NORM',tgt_info)
      call set_rule2('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'LABEL',1,
     &     tgt_info,val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &     tgt_info,val_label=(/'NORM'/))
c      if (h1bar) then
c        call set_dependency('F_MRCC_LAG','H1bar',tgt_info)
c        call set_dependency('F_MRCC_LAG','T-T1',tgt_info)
c        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
c     &       tgt_info,val_label=(/'Sbar    ','H1bar','T-T1 ','C0   '/))
c      else
        call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &       tgt_info,val_label=(/'Sbar','H   ','S   ','C0  '/))
c      end if
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/maxcom/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MODE',1,
     &     tgt_info,val_str='---')
        call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &     tgt_info,val_str='ic-MRCC-F12 Lagrangian')
c      if (optref.eq.-1.or.optref.eq.-2) then
c        call set_dependency('F_MRCC_F12_LAG','E(MR)',tgt_info)
c        call set_rule2('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,tgt_info)
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'LABEL',1,
c     &       tgt_info,
c     &       val_label=(/'F_MRCC_F12_LAG'/))
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'OP_RES',1,
c     &       tgt_info,
c     &       val_label=(/'NORM'/))
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &       tgt_info,
c     &       val_label=(/'E(MR)','C0^+ ','C0   '/))
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'IDX_SV',3,
c     &       tgt_info,
c     &       val_int=(/2,3,4/))
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'FAC',1,
c     &       tgt_info,val_rl8=(/-1d0/))
c        call set_arg('F_MRCC_F12_LAG',EXPAND_OP_PRODUCT,'NEW',1,
c     &       tgt_info,val_log=(/.false./))
c      end if
c      if (h1bar) then
c        call set_rule2('F_MRCC_LAG',REPLACE,tgt_info)
c        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_MRCC_LAG'/))
c        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_MRCC_LAG'/))
c        call set_arg('F_MRCC_LAG',REPLACE,'OP_LIST',2,tgt_info,
c     &       val_label=(/'T-T1','T   '/))
c      end if
      call set_rule2('F_MRCC_F12_LAG',SELECT_SPECIAL,tgt_info)
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &     tgt_info,val_label=(/'F_MRCC_F12_LAG'/))
c      if (h1bar) then
c        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,
c     &     tgt_info,val_label=(/'H1bar','T    '/))
c      else
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'OPERATORS',2,
     &     tgt_info,val_label=(/'H','S'/))
c      end if
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='CHECK_FAC')
c      if (h1bar) then
c        call set_dependency('F_MRCC_LAG','F_H1bar',tgt_info)
c        if (h1bar_maxp.lt.4) then
c          ! expand formal part of H1bar (store H1bar blks. only up to P^2)
c          call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
c          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
c     &         val_label=(/'F_MRCC_LAG'/))
c          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
c     &         val_label=(/'F_MRCC_LAG'/))
c          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
c     &         val_label=(/'F_H1barformal'/))
c          call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
c     &         val_int=(/2/))
c        end if
c      end if

      call set_rule2('F_MRCC_F12_LAG',REPLACE,tgt_info)
      call set_dependency('F_MRCC_F12_LAG','T',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','L',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','R12',tgt_info)
      call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',REPLACE,'OP_LIST',8,tgt_info,
     &     val_label=(/'S    ','R12  ','Sbar ','R12^+',
     &                 'S    ','T    ','Sbar ','L    '/))

      call set_rule2('F_MRCC_F12_LAG',SELECT_SPECIAL,tgt_info)
      call set_dependency('F_MRCC_F12_LAG','Favg',tgt_info)
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &     tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'OPERATORS',5,
     &     tgt_info,
     &     val_label=(/'L   ','H   ','T   ',
     &                 'R12 ','Favg'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC_F12')
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='MRCC')

c dbg
c      call set_rule2('F_MRCC_F12_LAG',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_F12_LAG',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_F12_LAG'/))
c dbgend

      ! identify F12 intermediates
      call set_dependency('F_MRCC_F12_LAG','VINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','BINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','BhINT_R12',tgt_info)
      !call set_dependency('F_MRCC_F12_LAG','XhINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','XINT_R12',tgt_info)
      if(.not.notrunc) then
      call set_dependency('F_MRCC_F12_LAG','CINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','C1_formal',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','Vring_formal',tgt_info)
      end if

      call set_rule2('F_MRCC_F12_LAG',FACTOR_OUT,tgt_info)
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      if(notrunc) then
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'INTERM',5,tgt_info,
     &     val_label=(/'BINT_R12','BhINT_R12','XINT_R12',
     &                 'VINT_R12','VINT_R12^+'/))
      else
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'INTERM',10,tgt_info,
     &     val_label=(/'BINT_R12','BhINT_R12','XINT_R12',
     &                 'VINT_R12','VINT_R12^+',
     &                 'CINT_R12','CINT_R12^+',
     &                 'Vring_formal','Vring_formal^+',
     &                 'C1_formal'/))
      end if

      if(notrunc) then
        call set_rule2('F_MRCC_F12_LAG',REPLACE,tgt_info)
        call set_dependency('F_MRCC_F12_LAG','R12-INT',tgt_info)
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'OP_LIST',4,tgt_info,
     &       val_label=(/'R12      ','R12-INT  ',
     &                   'R12^+    ','R12-INT^+'/))
      else
        Call set_rule2('F_MRCC_F12_LAG',INVARIANT,tgt_info)
        call set_arg('F_MRCC_F12_LAG',INVARIANT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',INVARIANT,'OP_RES',1,tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_MRCC_F12_LAG',INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/'R12  ','R12^+'/))
        call set_arg('F_MRCC_F12_LAG',INVARIANT,'TITLE',1,tgt_info,
     &       val_str='MRCC-R12 Lagrangian for pert. eval.')
      end if

      ! precursor for combining selected PP contractions into one step
      call add_target2('F_prePPint_F12',.false.,tgt_info)
      call set_dependency('F_prePPint_F12','F_OMG_F12',tgt_info)
      call set_dependency('F_prePPint_F12','H_PP',tgt_info)
      call set_rule2('F_prePPint_F12',REPLACE,tgt_info)
      call set_arg('F_prePPint_F12',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_prePPint_F12',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG_F12'/))
c      if (h1bar) then
c        call set_arg('F_prePPint_F12',REPLACE,'OP_LIST',2,tgt_info,
c     &       val_label=(/'H1bar','H_PP '/))
c      else
        call set_arg('F_prePPint_F12',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_PP'/))
c      end if
      call set_rule2('F_prePPint_F12',INVARIANT,tgt_info)
      call set_arg('F_prePPint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_prePPint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_prePPint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
c      if (h1bar) then
c        call set_arg('F_prePPint_F12',INVARIANT,'OPERATORS',2,tgt_info,
c     &       val_label=(/'H1bar','H    '/))
c      else
        call set_arg('F_prePPint_F12',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
c      end if
      call set_arg('F_prePPint_F12',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Precursor for INT_PP')

c dbg
c      call set_rule2('F_prePPint_F12',PRINT_FORMULA,tgt_info)
c      call set_arg('F_prePPint_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_prePPint_F12'/))
c dbgend

      call add_target2('F_PPint_F12',.false.,tgt_info)
      call set_dependency('F_PPint_F12','F_prePPint_F12',tgt_info)
      call set_dependency('F_PPint_F12','INT_PP',tgt_info)
      call set_rule2('F_PPint_F12',DERIVATIVE,tgt_info)
      call set_arg('F_PPint_F12',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PPint_F12'/))
      call set_arg('F_PPint_F12',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_PPint_F12',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP'/))
      call set_arg('F_PPint_F12',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'H_PP'/))
      call set_rule2('F_PPint_F12',INVARIANT,tgt_info)
      call set_arg('F_PPint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PPint_F12'/))
      call set_arg('F_PPint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_PPint_F12'/))
      call set_arg('F_PPint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP'/))
      call set_arg('F_PPint_F12',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('F_PPint_F12',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_PPint_F12',PRINT_FORMULA,tgt_info)
c      call set_arg('F_PPint_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG_F12'/))
c dbgend

      ! precursor for combining selected HH contractions into one step
      call add_target2('F_preHHint_F12',.false.,tgt_info)
      call set_dependency('F_preHHint_F12','F_OMG_F12',tgt_info)
      call set_dependency('F_preHHint_F12','H_HH',tgt_info)
      call set_rule2('F_preHHint_F12',REPLACE,tgt_info)
      call set_arg('F_preHHint_F12',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_preHHint_F12',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG_F12'/))
c      if (h1bar) then
c        call set_arg('F_preHHint_F12',REPLACE,'OP_LIST',2,tgt_info,
c     &       val_label=(/'H1bar','H_HH '/))
c      else 
        call set_arg('F_preHHint_F12',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_HH'/))
c      end if
      call set_rule2('F_preHHint_F12',INVARIANT,tgt_info)
      call set_arg('F_preHHint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_preHHint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_preHHint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
c      if (h1bar) then
c        call set_arg('F_preHHint_F12',INVARIANT,'OPERATORS',2,tgt_info,
c     &       val_label=(/'H    ','H1bar'/))
c      else 
        call set_arg('F_preHHint_F12',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
c      end if
      call set_arg('F_preHHint_F12',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Precursor for INT_HH')
c dbg
c      call set_rule2('F_preHHint_F12',PRINT_FORMULA,tgt_info)
c      call set_arg('F_preHHint_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_preHHint_F12'/))
c dbgend   
     
      call add_target2('F_HHint_F12',.false.,tgt_info)
      call set_dependency('F_HHint_F12','F_preHHint_F12',tgt_info)
      call set_dependency('F_HHint_F12','INT_HH',tgt_info)
      call set_rule2('F_HHint_F12',DERIVATIVE,tgt_info)
      call set_arg('F_HHint_F12',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_HHint_F12'/))
      call set_arg('F_HHint_F12',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_HHint_F12',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_HH'/))
      call set_arg('F_HHint_F12',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'H_HH'/))
      call set_rule2('F_HHint_F12',INVARIANT,tgt_info)
      call set_arg('F_HHint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_HHint_F12'/))
      call set_arg('F_HHint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_HHint_F12'/))
      call set_arg('F_HHint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_HH'/))
      call set_arg('F_HHint_F12',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('F_HHint_F12',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_HHint_F12',PRINT_FORMULA,tgt_info)
c      call set_arg('F_HHint_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HHint_F12'/))
c dbgend

      ! Lagrangian without Lambda...
      call add_target2('F_E_F12_C0',.false.,tgt_info)
      call set_dependency('F_E_F12_C0','F_MRCC_F12_LAG',tgt_info)
      call set_dependency('F_E_F12_C0','L',tgt_info)
      call set_rule2('F_E_F12_C0',INVARIANT,tgt_info)
      call set_arg('F_E_F12_C0',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E_F12_C0'/))
      call set_arg('F_E_F12_C0',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_E_F12_C0',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_E_F12_C0',INVARIANT,'OPERATORS',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_E_F12_C0',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC-F12 energy expression for C0 equations')

      ! Residual for C0
      call add_target2('F_OMG_F12_C0',.false.,tgt_info)
      call set_dependency('F_OMG_F12_C0','F_E_F12_C0',tgt_info)
      call set_dependency('F_OMG_F12_C0','A_C0',tgt_info)
      call set_rule2('F_OMG_F12_C0',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_F12_C0',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_F12_C0'/))
      call set_arg('F_OMG_F12_C0',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E_F12_C0'/))
      call set_arg('F_OMG_F12_C0',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C0'/))
      call set_arg('F_OMG_F12_C0',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'C0^+'/))

      ! Residual for C0
      call add_target2('FOPT_OMG_F12_C0',.false.,tgt_info)
      call set_dependency('FOPT_OMG_F12_C0','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG_F12_C0','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_OMG_F12_C0','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG_F12_C0',mel_ham,tgt_info)
      call set_dependency('FOPT_OMG_F12_C0','F_OMG_F12_C0',tgt_info)
      call set_dependency('FOPT_OMG_F12_C0','DEF_ME_A_C0',tgt_info)
      call set_rule2('FOPT_OMG_F12_C0',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_OMG_F12_C0',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_A_C0'/))
      call set_arg('FOPT_OMG_F12_C0',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'A_C0'/))
      call set_rule2('FOPT_OMG_F12_C0',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG_F12_C0',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG_F12_C0'/))
      call set_arg('FOPT_OMG_F12_C0',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_OMG_F12_C0'/))

c      ! Residual part of Lagrangian
c      call add_target2('F_LAG_L',.false.,tgt_info)
c      call set_dependency('F_LAG_L','F_MRCC_LAG',tgt_info)
c      call set_dependency('F_LAG_L','NORM',tgt_info)
c      call set_dependency('F_LAG_L','Sbar',tgt_info)
c      call set_rule2('F_LAG_L',SELECT_TERMS,tgt_info)
c      call set_arg('F_LAG_L',SELECT_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_LAG_L'/))
c      call set_arg('F_LAG_L',SELECT_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_LAG_L',SELECT_TERMS,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_LAG_L',SELECT_TERMS,'OP_INCL',1,tgt_info,
c     &     val_label=(/'Sbar'/))
c      ! Cumulant approximation?
c      if (maxcum.gt.0) then
c        ! Factor out reduced density matrices
c        call set_rule2('F_LAG_L',FACTOR_OUT,tgt_info)
c        call set_dependency('F_LAG_L','F_DENS0',tgt_info)
c        call set_arg('F_LAG_L',FACTOR_OUT,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        call set_arg('F_LAG_L',FACTOR_OUT,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        call set_arg('F_LAG_L',FACTOR_OUT,'INTERM',1,tgt_info,
c     &       val_label=(/'F_DENS0'/))
c        ! Expand density matrices in terms of cumulants
c        call set_rule2('F_LAG_L',EXPAND,tgt_info)
c        call set_arg('F_LAG_L',EXPAND,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        call set_arg('F_LAG_L',EXPAND,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        if (maxcom.le.2.and.cum_appr_mode.ge.2.and.
c     &      maxcum.ge.2.and.maxcum.le.4) then
c          call set_dependency('F_LAG_L','F_DENS_appr',tgt_info)
c          call set_arg('F_LAG_L',EXPAND,'INTERM',1,tgt_info,
c     &         val_label=(/'F_DENS_appr'/))
c        else
c          call set_dependency('F_LAG_L','F_DENS',tgt_info)
c          call set_arg('F_LAG_L',EXPAND,'INTERM',1,tgt_info,
c     &         val_label=(/'F_DENS'/))
c        end if
c        call set_arg('F_LAG_L',EXPAND,'IMODE',1,tgt_info,
c     &       val_int=(/2/)) ! keep low-rank RDMs
cc        ! Delete disconnected terms (cum. only connected to L)
cc        call set_rule2('F_LAG_L',SELECT_SPECIAL,tgt_info)
cc        call set_arg('F_LAG_L',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
cc     &       val_label=(/'F_LAG_L'/))
cc        call set_arg('F_LAG_L',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
cc     &       val_label=(/'F_LAG_L'/))
cc        call set_arg('F_LAG_L',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
cc     &       val_label=(/'CUM','L'/))
cc        call set_arg('F_LAG_L',SELECT_SPECIAL,'TYPE',1,tgt_info,
cc     &       val_str='MRCC3')
c      end if
cc dbg
cc        call set_rule2('F_LAG_L',PRINT_FORMULA,tgt_info)
cc        call set_arg('F_LAG_L',PRINT_FORMULA,'LABEL',1,tgt_info,
cc     &       val_label=(/'F_LAG_L'/))
cc dbgend

      ! Residual
      call add_target2('F_OMG_F12',.false.,tgt_info)
      call set_dependency('F_OMG_F12','F_MRCC_F12_LAG',tgt_info)
      call set_dependency('F_OMG_F12','OMG',tgt_info)
      call set_rule2('F_OMG_F12',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_F12',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_F12'/))
c      if (maxcum.gt.0) then
c        call set_dependency('F_OMG','F_LAG_L',tgt_info)
c        call set_arg('F_OMG',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c      else
        call set_arg('F_OMG_F12',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
c      end if
      call set_arg('F_OMG_F12',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_OMG_F12',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))
c      if (tfix.eq.0) then ! tfix>0: contains both T and Tfix
c        call set_rule2('F_OMG',SELECT_SPECIAL,tgt_info)
c        call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_OMG'/))
c        call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_OMG'/))
c        call set_arg('F_OMG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
c     &       val_label=(/'H','T'/))
c        call set_arg('F_OMG',SELECT_SPECIAL,'TYPE',1,tgt_info,
c     &       val_str='MRCC2')
c        call set_arg('F_OMG',SELECT_SPECIAL,'MODE',1,tgt_info,
c     &       val_str='CHECK    X')
c      else if (maxit.eq.1) then ! T=0
c        call set_rule2('F_OMG',INVARIANT,tgt_info)
c        call set_arg('F_OMG',INVARIANT,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_OMG'/))
c        call set_arg('F_OMG',INVARIANT,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_OMG'/))
c        call set_arg('F_OMG',INVARIANT,'OP_RES',1,tgt_info,
c     &       val_label=(/'OMG'/))
c        call set_arg('F_OMG',INVARIANT,'OPERATORS',1,tgt_info,
c     &       val_label=(/'T'/))
c        call set_arg('F_OMG',INVARIANT,'TITLE',1,tgt_info,
c     &       val_str='Higher-order residual for first iteration')
c      end if
c      if (tfix.gt.0) then
c        call set_rule2('F_OMG_F12',PRINT_FORMULA,tgt_info)
c        call set_arg('F_OMG_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &       val_label=(/'F_OMG_F12'/))
c      end if

      ! Energy
      call add_target2('F_MRCC_F12_E',.false.,tgt_info)
      call set_dependency('F_MRCC_F12_E','F_MRCC_F12_LAG',tgt_info)
      call set_dependency('F_MRCC_F12_E','E(MR)',tgt_info)
      call set_dependency('F_MRCC_F12_E','L',tgt_info)
      call set_rule2('F_MRCC_F12_E',INVARIANT,tgt_info)
      call set_arg('F_MRCC_F12_E',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_E'/))
      call set_arg('F_MRCC_F12_E',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_E',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
c      if (tfix.eq.0) then
c        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',2,tgt_info,
c     &       val_label=(/'L    ','E(MR)'/))
        call set_arg('F_MRCC_F12_E',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'L'/))
c      else if (maxit.gt.1) then
c        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',3,tgt_info,
c     &       val_label=(/'L     ','Tfix^+','E(MR) '/))
c      else
c        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',4,tgt_info,
c     &       val_label=(/'L     ','Tfix^+','E(MR) ','T     '/))
c      end if
      call set_arg('F_MRCC_F12_E',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC-F12 energy expression')
c      if (tfix.ne.0) then
c        call set_rule2('F_MRCC_E',SELECT_SPECIAL,tgt_info)
c        call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_MRCC_E'/))
c        call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_MRCC_E'/))
c        call set_arg('F_MRCC_E',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
c     &       val_label=(/'H','T'/))
c        call set_arg('F_MRCC_E',SELECT_SPECIAL,'TYPE',1,tgt_info,
c     &       val_str='MRCC2')
c        call set_arg('F_MRCC_E',SELECT_SPECIAL,'MODE',1,tgt_info,
c     &       val_str='CHECK    X')
c      end if
c dbg
      call set_rule2('F_MRCC_F12_E',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_F12_E',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_E'/))
c group the energy equation for different density matrices
      call add_target2('F_1',.false.,tgt_info)
      call set_dependency('F_1','C0',tgt_info)
      call set_dependency('F_1','1scal',tgt_info)
      call set_rule2('F_1',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_1',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_1'/))
      call set_arg('F_1',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'1scal'/))
      call set_arg('F_1',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &     tgt_info,
     &     val_label=(/'C0^+','C0  '/))
      call set_arg('F_1',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &     val_int=(/1,2/))

      call add_target2('FORM_E_F12',.false.,tgt_info)
      call set_dependency('FORM_E_F12','F_MRCC_F12_E',tgt_info)
c      call set_dependency('FORM_E_F12','F_DENS0',tgt_info)
      call set_dependency('FORM_E_F12','F_1',tgt_info)
      call set_rule2('FORM_E_F12',INVARIANT,tgt_info)
      call set_arg('FORM_E_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'FORM_E_F12'/))
      call set_arg('FORM_E_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_E'/))
      call set_arg('FORM_E_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
        call set_arg('FORM_E_F12',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
      call set_arg('FORM_E_F12',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='F12 energy expression')
      call set_rule2('FORM_E_F12',FACTOR_OUT,tgt_info)
      call set_arg('FORM_E_F12',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'FORM_E_F12'/))
      call set_arg('FORM_E_F12',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'FORM_E_F12'/))
      call set_arg('FORM_E_F12',FACTOR_OUT,'INTERM',1,tgt_info,
     &     val_label=(/'F_1    '/))

      call set_rule2('FORM_E_F12',PRINT_FORMULA,tgt_info)
      call set_arg('FORM_E_F12',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'FORM_E_F12'/))
      call set_rule2('FORM_E_F12',TEX_FORMULA,tgt_info)
      call set_arg('FORM_E_F12',TEX_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'FORM_E_F12'/))
      call set_arg('FORM_E_F12',TEX_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='formula.tex')
c dbgend

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! Residual
      call add_target2('FOPT_OMG_F12',.false.,tgt_info)
      call set_dependency('FOPT_OMG_F12','EVAL-R12-INTER',tgt_info)
      call set_dependency('FOPT_OMG_F12','Vring-EVAL',tgt_info)
      if(.not.notrunc) 
     &   call set_dependency('FOPT_OMG_F12','C1-EVAL',tgt_info)
      call set_dependency('FOPT_OMG_F12','F_OMG_F12',tgt_info)
      call set_dependency('FOPT_OMG_F12','F_MRCC_F12_E',tgt_info)
      call set_dependency('FOPT_OMG_F12','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG_F12','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_OMG_F12','DEF_ME_OMG',tgt_info)
      call set_dependency('FOPT_OMG_F12','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG_F12',mel_ham,tgt_info)
c      if (optref.eq.-1.or.optref.eq.-2) then
c        call set_dependency('FOPT_OMG','F_OMG_C0',tgt_info)
c        call set_dependency('FOPT_OMG','DEF_ME_A_C0',tgt_info)
cc      call set_dependency('FOPT_OMG','DEF_ME_1v',tgt_info)
c        call set_rule2('FOPT_OMG',ASSIGN_ME2OP,tgt_info)
c        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &             val_label=(/'ME_A_C0'/))
c        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &             val_label=(/'A_C0'/))
c      end if
c      if (tfix.gt.0)
c     &    call set_dependency('FOPT_OMG','DEF_ME_Tfix',tgt_info)
      call set_rule2('FOPT_OMG_F12',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG_F12',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG_F12'/))
c      if ((maxp.ge.2.or.maxh.ge.2).and.tfix.eq.0) then
      if (maxp.ge.2.or.maxh.ge.2) then
        labels(1:20)(1:len_target_name) = ' '
        ndef = 0
c        if (maxp.ge.2.and.h1bar_maxp.lt.4) then
c          call set_dependency('FOPT_OMG_F12','F_PP0int_F12',tgt_info)
c          call set_dependency('FOPT_OMG_F12','DEF_ME_INT_PP0',tgt_info)
c          ndef = ndef + 1
c          labels(ndef) = 'F_PP0int'
c        end if
c        if (maxp.ge.2.and.(maxh.gt.0.or.h1bar_maxp.gt.2)) then
        if (maxp.ge.2.and.maxh.gt.0) then
          call set_dependency('FOPT_OMG_F12','F_PPint_F12',tgt_info)
          call set_dependency('FOPT_OMG_F12','DEF_ME_INT_PP',tgt_info)
          ndef = ndef + 1
          labels(ndef) = 'F_PPint_F12'
        end if
        if (maxh.ge.2) then
          call set_dependency('FOPT_OMG_F12','F_HHint_F12',tgt_info)
          call set_dependency('FOPT_OMG_F12','DEF_ME_INT_HH',tgt_info)
          ndef = ndef + 1
          labels(ndef) = 'F_HHint_F12'
        end if
c dbg
!        call set_dependency('FOPT_OMG','F_INT_HT2',tgt_info)
!        call set_dependency('FOPT_OMG','DEF_ME_INT_HT2',tgt_info)
!        call set_dependency('FOPT_OMG','F_INT_T2H',tgt_info)
!        call set_dependency('FOPT_OMG','DEF_ME_INT_T2H',tgt_info)
!        call set_dependency('FOPT_OMG','F_INT_D',tgt_info)
!        call set_dependency('FOPT_OMG','DEF_ME_INT_D',tgt_info)
c        labels(ndef+1) = 'F_INT_HT2'
c        labels(ndef+2) = 'F_INT_T2H' 
!        labels(ndef+1) = 'F_INT_D'
!        ndef = ndef + 1!3
c dbg
        call set_arg('FOPT_OMG_F12',OPTIMIZE,'INTERM',ndef,tgt_info,
     &               val_label=labels(1:ndef))
      end if
      labels(1:20)(1:len_target_name) = ' '
      ndef = 0
c      if (maxp.ge.2.and.tfix.eq.0) then
c      if (h1bar) then
c        call set_dependency('FOPT_OMG','F_H1bar',tgt_info)
c        call set_dependency('FOPT_OMG','DEF_ME_H1bar',tgt_info)
c        labels(ndef+1) = 'F_H1bar'
c        ndef = ndef + 1
c      end if
      labels(ndef+1) = 'F_MRCC_F12_E'
      ndef = ndef + 1
      labels(ndef+1) = 'F_OMG_F12'
      ndef = ndef + 1
c      if (optref.eq.-1.or.optref.eq.-2) then
c        labels(ndef+1) = 'F_OMG_C0'
c        ndef = ndef + 1
c      end if
      call set_arg('FOPT_OMG_F12',OPTIMIZE,'LABELS_IN',ndef,tgt_info,
     &             val_label=labels(1:ndef))

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Solve MR coupled cluster equations
      call add_target2('SOLVE_MRCC_F12',solve,tgt_info)
c      call add_target2('SOLVE_MRCC_F12',.false.,tgt_info)
      call set_dependency('SOLVE_MRCC_F12','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','FOPT_OMG_F12',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('SOLVE_MRCC_F12',trim(dia_label),tgt_info)
      call set_dependency('SOLVE_MRCC_F12','EVAL_D',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','DEF_ME_Dtrdag',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','FOPT_T',tgt_info)
c      select case(prc_type)
c      case(-1) !do nothing: use old preconditioner file!
c        call warn('set_ic_mrcc_targets','Using old preconditioner file')
c      case(0,3)
c dbg hybrid preconditioner
c        call warn('set_ic_mrcc_targets','Using hybrid preconditioner')
c dbgend
        call set_dependency('SOLVE_MRCC_F12','EVAL_Atr',tgt_info)
c      case(1)
c        call set_dependency('SOLVE_MRCC','EVAL_A_Ttr',tgt_info)
c      case(2)
c        call set_dependency('SOLVE_MRCC','PREC_diag',tgt_info)
c        ! we misused ME_Dtrdag. Now get correct one:
c        call set_rule2('SOLVE_MRCC',REORDER_MEL,tgt_info)
c        call set_arg('SOLVE_MRCC',REORDER_MEL,'LIST_RES',1,tgt_info,
c     &               val_label=(/'ME_Dtrdag'/))
c        call set_arg('SOLVE_MRCC',REORDER_MEL,'LIST_IN',1,tgt_info,
c     &               val_label=(/'ME_Dinv'/))
c        call set_arg('SOLVE_MRCC',REORDER_MEL,'FROMTO',1,tgt_info,
c     &               val_int=(/13/))
c        call set_arg('SOLVE_MRCC',REORDER_MEL,'ADJOINT',1,tgt_info,
c     &               val_log=(/.true./))
c      case default
c        call quit(1,'set_ic_mrcc_targets','unknown prc_type')
c      end select
      if (optref.ne.0) then
        call me_list_label(dia_label2,mel_dia,orb_info%lsym,
     &                     0,0,0,.false.)
        dia_label2 = trim(dia_label2)//'C0'
        call set_dependency('SOLVE_MRCC_F12',trim(dia_label2),tgt_info)
        if (optref.ne.-1.and.optref.ne.-2) 
     &     call set_dependency('SOLVE_MRCC_F12','FOPT_OMG_C0'
     &          ,tgt_info)
        call set_dependency('SOLVE_MRCC_F12','DEF_ME_Dproj',tgt_info)
      end if
c      do icnt = 1, max(1,optref)
      call set_rule2('SOLVE_MRCC_F12',SOLVENLEQ,tgt_info)
c      if (optref.lt.0) then
c        if (preopt) then
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
c     &         val_label=(/'ME_T'/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
c     &         val_str='TRF')
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
c     &         val_label=(/'ME_OMG'/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
c     &         val_label=(/trim(dia_label)/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
c     &       val_label=(/'ME_E(MR)'/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
c     &       val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag'/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',1,tgt_info,
c     &       val_label=(/'FOPT_T'/))
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
c     &         val_label=(/'FOPT_OMG'/))
c          call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
c        end if
c      end if
c      if (optref.eq.-1.or.optref.eq.-2) then
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',2,tgt_info,
c     &       val_label=(/'ME_T ','ME_C0'/))
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
c     &       val_str='TRF/NRM')
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',2,tgt_info,
c     &       val_label=(/'ME_OMG ','ME_A_C0'/))
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',2,tgt_info,
c     &       val_label=(/trim(dia_label),trim(dia_label2)/))
c      else
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_T'/))
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF')
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &       val_label=(/'ME_OMG'/))
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label)/))
c      end if
      call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_E',1,tgt_info,
     &     val_label=(/'ME_E(MR)'/))
      if (optref.ne.0.and.update_prc.gt.0) then
c        if (tred.eq.0) then
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag',
     &                 'ME_Dproj ',
     &                 'ME_D     ','ME_Dinv  ',
     &                 'ME_A     '/))
c        else if (ex_t3red) then
c        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
c        call set_dependency('SOLVE_MRCC','FOPT_T(3)red',tgt_info)
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',9,tgt_info,
c     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
c     &                 'ME_Dproj  ',
c     &                 'ME_D      ','ME_Dinv   ',
c     &                 'ME_A      ','ME_T(2)red','ME_T(3)red'/))
c        else
c        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
c        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',8,tgt_info,
c     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
c     &                 'ME_Dproj  ',
c     &                 'ME_D      ','ME_Dinv   ',
c     &                 'ME_A      ','ME_T(2)red'/))
c        end if
      else if (optref.ne.0.and.update_prc.le.0) then
cc      if (optref.ne.0) then
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',6,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   '/))
      else
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag '/))
      end if
      if (optref.ne.0) then
        if (update_prc.gt.0) then
c          if (tred.eq.0) then
          call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'FORM_SPC',3,tgt_info,
     &         val_label=(/'FOPT_T  ','FOPT_D  ','FOPT_Atr'/))
c          else if (ex_t3red) then
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',5,tgt_info,
c     &         val_label=(/'FOPT_T      ','FOPT_D      ',
c     &                     'FOPT_Atr    ',
c     &                     'FOPT_T(2)red','FOPT_T(3)red'/))
c          else
c          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',4,tgt_info,
c     &         val_label=(/'FOPT_T      ','FOPT_D      ',
c     &                     'FOPT_Atr    ',
c     &                     'FOPT_T(2)red'/))
c          end if
        else
          call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'FORM_SPC',2,tgt_info,
     &         val_label=(/'FOPT_T','FOPT_D'/))
        end if
      else
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'FORM_SPC',1,tgt_info,
     &       val_label=(/'FOPT_T'/))
      end if
      call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'FORM',1,tgt_info,
     &     val_label=(/'FOPT_OMG_F12'/))
c      if (optref.gt.0.and.icnt.ne.optref) then !not in last iteration
c        call set_rule2('SOLVE_MRCC',SOLVEEVP,tgt_info)
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_OPT',1,tgt_info,
c     &       val_label=(/'ME_C0'/))
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'MODE',1,tgt_info,
c     &       val_str='DIA')
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'N_ROOTS',1,tgt_info,
c     &       val_int=(/ciroot/))
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_MVP',1,tgt_info,
c     &       val_label=(/'A_C0'/))
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_PRC',1,tgt_info,
c     &       val_label=(/trim(dia_label2)/))
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_SVP',1,tgt_info,
c     &     val_label=(/'C0'/))
c        call set_arg('SOLVE_MRCC',SOLVEEVP,'FORM',1,tgt_info,
c     &       val_label=(/'FOPT_OMG_C0'/))
cc dbg
c        call form_parameters(-1,parameters,2,
c     &       'CI coefficients :',0,'LIST')
c        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &       'ME_C0',1,0,
c     &       parameters,2,tgt_info)
c        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &       'FOPT_REF_S(S+1)',1,0,
c     &       parameters,0,tgt_info)
c        call set_rule2('SOLVE_MRCC',PRINT_MEL,tgt_info)
c        call set_arg('SOLVE_MRCC',PRINT_MEL,'LIST',1,tgt_info,
c     &       val_label=(/'ME_S(S+1)'/))
c        call set_arg('SOLVE_MRCC',PRINT_MEL,'COMMENT',1,tgt_info,
c     &       val_str='Spin expectation value <C0| S^2 |C0> :')
c        call set_arg('SOLVE_MRCC',PRINT_MEL,'FORMAT',1,tgt_info,
c     &       val_str='SCAL F20.12')
c        call set_arg('SOLVE_MRCC',PRINT_MEL,'CHECK_THRESH',1,tgt_info,
c     &       val_rl8=(/1d-2/))
c        call set_arg('SOLVE_MRCC',PRINT_MEL,'EXPECTED',1,tgt_info,
c     &       val_rl8=(/(dble(orb_info%imult**2)-1d0)/4d0/))
cc        call form_parameters(-1,parameters,2,
cc     &       'Spin expectation value <C0| S^2 |C0> :',0,'SCAL F20.12')
cc        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
cc     &       'S(S+1)',1,0,
cc     &       parameters,2,tgt_info)
cc dbgend
c      end if
c      end do
cc dbg
      if (optref.ne.0) then
        call form_parameters(-1,parameters,2,
     &       'final CI coefficients :',0,'LIST')
        call set_rule('SOLVE_MRCC_F12',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call form_parameters(-1,parameters,2,
     &       'Spin expectation value <C0| S^2 |C0> :',0,'SCAL F20.12')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_S(S+1)',1,0,
     &       parameters,2,tgt_info)
      end if
c      call set_dependency('SOLVE_MRCC','FOPT_T_S2',tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
c     &     'ME_S(S+1)',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &     'FOPT_T_S2',1,0,
c     &     parameters,0,tgt_info)
c      call set_dependency('SOLVE_MRCC','FOPT_T_NORM',tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
c     &     'ME_NORM',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &     'FOPT_T_NORM',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule2('SOLVE_MRCC',SCALE_COPY,tgt_info)
c      call set_arg('SOLVE_MRCC',SCALE_COPY,'LIST_RES',1,tgt_info,
c     &             val_label=(/'ME_S(S+1)'/))
c      call set_arg('SOLVE_MRCC',SCALE_COPY,'LIST_INP',1,tgt_info,
c     &             val_label=(/'ME_NORM'/))
c      call set_arg('SOLVE_MRCC',SCALE_COPY,'FAC',1,tgt_info,
c     &             val_rl8=(/1d0/))
c      call set_arg('SOLVE_MRCC',SCALE_COPY,'MODE',1,tgt_info,
c     &             val_str='precond')
c      call form_parameters(-1,parameters,2,
c     &     'Spin expectation value <C0|T^+ S^2 T|C0>/<C0|T^+ T|C0> :',
c     &     0,'SCAL F20.12')
c      call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &     'ME_S(S+1)',1,0,
c     &     parameters,2,tgt_info)
c      if (.not.h1bar.and.tfix.eq.0.and.maxcum.le.0) then
c        call set_dependency('SOLVE_MRCC','FOPT_MRCC_S(S+1)',tgt_info)
c        call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
c     &       'ME_S(S+1)',1,0,
c     &       parameters,0,tgt_info)
c        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &       'FOPT_MRCC_S(S+1)',1,0,
c     &       parameters,0,tgt_info)
c        call form_parameters(-1,parameters,2,
c     &       'Spin expectation value <C0| T^+ e^-T S^2 e^T |C0> :',
c     &       0,'SCAL F20.12')
c        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &       'ME_S(S+1)',1,0,
c     &       parameters,2,tgt_info)
c      end if
cc dbgend
cc dbg
cc        call form_parameters(-1,parameters,2,
cc     &       'final T amplitudes :',0,'LIST')
cc        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
cc     &       'ME_T',1,0,
cc     &       parameters,2,tgt_info)
cc dbgend

      return
      end
