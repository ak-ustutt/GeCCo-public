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
     &     ntest =  00
      logical, parameter ::
     &     test_semi = .false.

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
     &     n_tred_cls, len_form, optref, isemi, ciroot,
     &     version(60), ivers, stndT(2,60), stndD(2,60), nsupT, nsupD,
     &     G_level, iexc, jexc, maxtt, iblk, jblk, kblk, prc_type,
     &     tred, nremblk, remblk(60), igasreo(3), ngas, lblk, ntrunc,
     &     tfix, maxit, t1ord, maxcum, cum_appr_mode, update_prc
      logical ::
     &     skip, preopt, project, first, Op_eqs,
     &     h1bar, htt, svdonly, fact_tt, ex_t3red, trunc, l_exist,
     &     oldref, solve, notrunc, eval_dens3, 
     &     semi_r12, semi_red, semi_intm,
     &     restart, prc_traf, use_u3
      character(len_target_name) ::
     &     dia_label, dia_label2,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     op_ht*3, f_ht*5, op_ht0to*6, f_ht0to*8, form_str*50,
     &     def_ht*10
      character(256) :: descr
      real(8) ::
     &     x_ansatz, prc_shift

      if (iprlvl.gt.0) write(lulog,*) 'setting icMRCC_F12 targets'

      msc = +1
      if (orb_info%ims.ne.0) msc = 0

      call get_argument_value('method.MR','svdonly',
     &     lval=svdonly)
      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)
      call get_argument_value('calculate.solve.non_linear','update_prc',
     &     ival=update_prc)
      call get_argument_value('calculate.solve.non_linear','restart',
     &     lval=restart)
      call get_argument_value('method.MRCC','maxcom_res',
     &     ival=maxcom)
      call get_argument_value('method.MRCC','maxcom_en',
     &     ival=maxcom_en)
      call get_argument_value('method.MRCC','H1bar',
     &     lval=h1bar)
      call get_argument_value('method.MRCC','eval_dens3',
     &     lval=eval_dens3)
      call get_argument_value('method.MR','prc_traf',
     &     lval=prc_traf)
      call get_argument_value('method.R12','notrunc',lval=notrunc)
      call get_argument_value('method.R12','use_U3',lval=use_u3)
      call get_argument_value('method.R12','semi_r12',lval=semi_r12)
      call get_argument_value('method.R12','semi_int',ival=isemi)
      if (semi_r12) call warn('set_ic_mrcc_f12_targets',
     &                'semi_r12 is obsolete!')
      if (isemi.eq.0) then
        ! semi_r12 = .false.  ! for the moment, we leave semi_r12
        semi_red = .false.
        semi_intm = .false. 
      else if (isemi.eq.1) then
        semi_r12 = .true.
        semi_red = .true.
        semi_intm = .true.
      else if (isemi.eq.2) then
        semi_r12 = .true.
        semi_red = .true.
        semi_intm = .false.
      else
        semi_r12 = .true.
        semi_red = .false.
        semi_intm = .false.
      end if
      skip = (is_keyword_set('calculate.skip_E').gt.0)
      solve = .not.svdonly.and..not.skip
      if (h1bar) call quit(1,'set_ic_mrcc_f12_targets',
     &                   'H1bar not yet available for F12 calculations')
      if (optref.ne.0.and.optref.ne.-3)
     &   call quit(1,'set_ic_mrcc_f12_targets','use optref=0 or -3')
      
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving excitation operator S (for icMRCC-F12)
      call add_target2('S',.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            if (iexc.eq.3.and.eval_dens3) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
          end do
        end do
      end do
      do ih = 0, min(maxh,2)
        do ip = 1, 0, -1
          ndef = ndef + 1
          occ_def(IHOLE,2,ndef) = ih
          occ_def(IPART,1,ndef) = ip
          occ_def(IVALE,2,ndef) = 2 - ih
          occ_def(IEXTR,1,ndef) = 2 - ip
        end do
      end do
      if (semi_r12) then
      ndef = ndef + 1
      occ_def(IHOLE,2,ndef) = 2
      occ_def(IVALE,1,ndef) = 1
      occ_def(IEXTR,1,ndef) = 1
      ndef = ndef + 1
      occ_def(IHOLE,2,ndef) = 1
      occ_def(IVALE,1,ndef) = 1
      occ_def(IVALE,2,ndef) = 1
      occ_def(IEXTR,1,ndef) = 1
      ndef = ndef + 1
      occ_def(IVALE,1,ndef) = 1
      occ_def(IVALE,2,ndef) = 2
      occ_def(IEXTR,1,ndef) = 1
      ndef = ndef + 1
      occ_def(IHOLE,2,ndef) = 1
      occ_def(IEXTR,1,ndef) = 1
      ndef = ndef + 1
      occ_def(IVALE,2,ndef) = 1
      occ_def(IEXTR,1,ndef) = 1
      end if

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

      ! multireference CC-F12 lagrangian
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
        call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &       tgt_info,val_label=(/'Sbar','H   ','S   ','C0  '/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/maxcom/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'MODE',1,
     &     tgt_info,val_str='---')
        call set_arg('F_MRCC_F12_LAG',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &     tgt_info,val_str='ic-MRCC-F12 Lagrangian')
      call set_rule2('F_MRCC_F12_LAG',SELECT_SPECIAL,tgt_info)
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &     tgt_info,val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'OPERATORS',2,
     &     tgt_info,val_label=(/'H','S'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='CHECK_FAC')

      call set_rule2('F_MRCC_F12_LAG',REPLACE,tgt_info)
      call set_dependency('F_MRCC_F12_LAG','T',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','L',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','R12',tgt_info)
      if (semi_r12) then
      call set_dependency('F_MRCC_F12_LAG','sR12-INT',tgt_info)
      call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',REPLACE,'OP_LIST',16,tgt_info,
     &     val_label=(/'S      ','R12si  ','Sbar   ','R12si^+',
     &                 'S      ','R12    ','Sbar   ','R12^+  ',
     &                 'S      ','sR12   ','Sbar   ','sR12^+ ',
     &                 'S      ','T      ','Sbar   ','L      '/))
      if (test_semi) then
        call set_dependency('F_MRCC_F12_LAG','sR12-INT',tgt_info)
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'sR12-INT'/))
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'sR12-INT^+'/))
      end if

      call set_rule2('F_MRCC_F12_LAG',SELECT_SPECIAL,tgt_info)
      call set_dependency('F_MRCC_F12_LAG','Favg',tgt_info)
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &     tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      if (test_semi) then
        call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'OPERATORS',7,
     &     tgt_info,
     &     val_label=(/'L      ','H      ','T      ',
     &                 'R12    ','Favg   ','R12si  ','Rsi-INT'/))
      else
        call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'OPERATORS',7,
     &     tgt_info,
     &     val_label=(/'L    ','H    ','T    ',
     &                 'R12  ','Favg ','R12si','sR12 '/))
      end if
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC_F12')
      descr='MRCC_SI'
      if (semi_red) descr='MRCC_SI_RED'
      call set_arg('F_MRCC_F12_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str=descr)
      else
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
      end if

      ! identify F12 intermediates
      call set_dependency('F_MRCC_F12_LAG','VINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','BINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','BhINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','XINT_R12',tgt_info)
      if(.not.notrunc) then
      call set_dependency('F_MRCC_F12_LAG','CINT_R12',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','C1_formal',tgt_info)
      call set_dependency('F_MRCC_F12_LAG','Vring_formal',tgt_info)
      if (use_u3)
     &call set_dependency('F_MRCC_F12_LAG','U3_formal',tgt_info)
      end if

      call set_rule2('F_MRCC_F12_LAG',FACTOR_OUT,tgt_info)
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_F12_LAG'/))
      if(notrunc) then
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'INTERM',5,tgt_info,
     &     val_label=(/'BINT_R12  ','BhINT_R12 ','XINT_R12  ',
     &                 'VINT_R12  ','VINT_R12^+'/))
      else if (.not.use_u3) then
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'INTERM',10,tgt_info,
     &     val_label=(/'BINT_R12      ','BhINT_R12     ',
     &                 'XINT_R12      ',
     &                 'VINT_R12      ','VINT_R12^+    ',
     &                 'CINT_R12      ','CINT_R12^+    ',
     &                 'Vring_formal  ','Vring_formal^+',
     &                 'C1_formal     '/))
      else
      call set_arg('F_MRCC_F12_LAG',FACTOR_OUT,'INTERM',12,tgt_info,
     &     val_label=(/'BINT_R12      ','BhINT_R12     ',
     &                 'XINT_R12      ',
     &                 'VINT_R12      ','VINT_R12^+    ',
     &                 'CINT_R12      ','CINT_R12^+    ',
     &                 'Vring_formal  ','Vring_formal^+',
     &                 'C1_formal     ',
     &                 'U3_formal     ','U3_formal^+   '/))
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
      if (use_u3) then
        ! re-expand U3
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'U3_formal'/))
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'U3_formal^+'/))
        ! re-expand C as well, otherwise we get into double counting trouble
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'CINT_R12'/))
        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'CINT_R12^+'/))
C        ! C1 intermediate? No -- there is no clash
C        call set_rule2('F_MRCC_F12_LAG',EXPAND,tgt_info)
C        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_RES',1,tgt_info,
C     &       val_label=(/'F_MRCC_F12_LAG'/))
C        call set_arg('F_MRCC_F12_LAG',EXPAND,'LABEL_IN',1,tgt_info,
C     &       val_label=(/'F_MRCC_F12_LAG'/))
C        call set_arg('F_MRCC_F12_LAG',EXPAND,'INTERM',1,tgt_info,
C     &       val_label=(/'C1_formal'/))
        ! and replace R12 by R12_INT
        call set_rule2('F_MRCC_F12_LAG',REPLACE,tgt_info)
        call set_dependency('F_MRCC_F12_LAG','R12-INT',tgt_info)
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'OP_LIST',4,tgt_info,
     &       val_label=(/'R12      ','R12-INT  ',
     &                   'R12^+    ','R12-INT^+'/))
      end if
      if (semi_r12) then
        call set_rule2('F_MRCC_F12_LAG',REPLACE,tgt_info)
        call set_dependency('F_MRCC_F12_LAG','ME_Rsi',tgt_info)
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
        call set_arg('F_MRCC_F12_LAG',REPLACE,'OP_LIST',4,tgt_info,
     &       val_label=(/'R12si    ','Rsi-INT  ',
     &                   'R12si^+  ','Rsi-INT^+'/))
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
        call set_arg('F_prePPint_F12',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_PP'/))
      call set_rule2('F_prePPint_F12',INVARIANT,tgt_info)
      call set_arg('F_prePPint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_prePPint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePPint_F12'/))
      call set_arg('F_prePPint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
        call set_arg('F_prePPint_F12',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
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
        call set_arg('F_preHHint_F12',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_HH'/))
      call set_rule2('F_preHHint_F12',INVARIANT,tgt_info)
      call set_arg('F_preHHint_F12',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_preHHint_F12',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preHHint_F12'/))
      call set_arg('F_preHHint_F12',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
        call set_arg('F_preHHint_F12',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
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

      ! Residual for C0  ---- NOTE: using target F_OMG_C0 instead
      !      (see set_ic_mrcc_targets.f; correct Lagrangian provided there! )
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
c dbg
c      call set_rule2('F_OMG_F12_C0',PRINT_FORMULA,tgt_info)
c      call set_arg('F_OMG_F12_C0',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG_F12_C0'/))
c dbgend
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

      ! Residual
      call add_target2('F_OMG_F12',.false.,tgt_info)
      call set_dependency('F_OMG_F12','F_MRCC_F12_LAG',tgt_info)
      call set_dependency('F_OMG_F12','OMG',tgt_info)
      call set_rule2('F_OMG_F12',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_F12',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_F12'/))
        call set_arg('F_OMG_F12',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
      call set_arg('F_OMG_F12',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_OMG_F12',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))

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
        call set_arg('F_MRCC_F12_E',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'L'/))
      call set_arg('F_MRCC_F12_E',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC-F12 energy expression')
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

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! Residual
      call add_target2('FOPT_OMG_F12',.false.,tgt_info)
      call set_dependency('FOPT_OMG_F12','EVAL-R12-INTER',tgt_info)
      call set_dependency('FOPT_OMG_F12','EVAL_D',tgt_info)
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
      call set_rule2('FOPT_OMG_F12',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG_F12',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG_F12'/))
      if (maxp.ge.2.or.maxh.ge.2) then
        labels(1:20)(1:len_target_name) = ' '
        ndef = 0
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
        call set_arg('FOPT_OMG_F12',OPTIMIZE,'INTERM',ndef,tgt_info,
     &               val_label=labels(1:ndef))
      end if
      labels(1:20)(1:len_target_name) = ' '
      ndef = 0
      if (semi_r12) then
        call set_dependency('FOPT_OMG_F12','make-sR12-INT',tgt_info)
c        call set_dependency('FOPT_OMG_F12','ME_sR12',tgt_info)
c        call set_dependency('FOPT_OMG_F12','F_DENS0',tgt_info)
c        labels(ndef+1) = 'F_DENS0'
c        ndef = ndef + 1
c        labels(ndef+1) = 'sR12-INT'
c        ndef = ndef + 1
      end if
      labels(ndef+1) = 'F_MRCC_F12_E'
      ndef = ndef + 1
      labels(ndef+1) = 'F_OMG_F12'
      ndef = ndef + 1
      call set_arg('FOPT_OMG_F12',OPTIMIZE,'LABELS_IN',ndef,tgt_info,
     &             val_label=labels(1:ndef))

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Solve MR coupled cluster F12 equations
      call add_target2('SOLVE_MRCC_F12',solve,tgt_info)
      call set_dependency('SOLVE_MRCC_F12','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','FOPT_OMG_F12',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('SOLVE_MRCC_F12',trim(dia_label),tgt_info)
      call set_dependency('SOLVE_MRCC_F12','EVAL_D',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','DEF_ME_Dtrdag',tgt_info)
      call set_dependency('SOLVE_MRCC_F12','FOPT_T',tgt_info)
        call set_dependency('SOLVE_MRCC_F12','EVAL_Atr',tgt_info)
      if (restart) ! project out redundant part (if sv_thr. changed)
     &   call set_dependency('SOLVE_MRCC_F12','EVAL_Tproj',tgt_info)
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
      call set_rule2('SOLVE_MRCC_F12',SOLVENLEQ,tgt_info)
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_T'/))
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF')
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &       val_label=(/'ME_OMG'/))
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label)/))
      call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_E',1,tgt_info,
     &     val_label=(/'ME_E(MR)'/))
      if (optref.ne.0.and.(update_prc.gt.0.or.prc_traf)) then
       if (prc_traf) then
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',8,tgt_info,
     &     val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag',
     &                 'ME_Dproj ',
     &                 'ME_D     ','ME_Dinv  ',
     &                 'ME_A     ','ME_Auni  '/))
       else
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag',
     &                 'ME_Dproj ',
     &                 'ME_D     ','ME_Dinv  ',
     &                 'ME_A     '/))
       end if
      else if (optref.ne.0.and.update_prc.le.0) then
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',6,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   '/))
      else
        call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag '/))
      end if
      if (optref.ne.0) then
        if (update_prc.gt.0.or.prc_traf) then
          call set_arg('SOLVE_MRCC_F12',SOLVENLEQ,'FORM_SPC',3,tgt_info,
     &         val_label=(/'FOPT_T  ','FOPT_D  ','FOPT_Atr'/))
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
      if (optref.ne.0) then
!        call form_parameters(-1,parameters,2,
!     &       'final CI coefficients :',0,'LIST')
!        call set_rule('SOLVE_MRCC_F12',ttype_opme,PRINT_MEL,
!     &       'ME_C0',1,0,
!     &       parameters,2,tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call form_parameters(-1,parameters,2,
     &       'Spin expectation value <C0| S^2 |C0> :',0,'SCAL F20.12')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_S(S+1)',1,0,
     &       parameters,2,tgt_info)
      end if

      return
      end
