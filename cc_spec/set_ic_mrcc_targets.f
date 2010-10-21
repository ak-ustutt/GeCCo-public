*----------------------------------------------------------------------*
      subroutine set_ic_mrcc_targets(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set targets for internally contracted MRCC
*
*     matthias 2010
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

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ndef, occ_def(ngastp,2,60),
     &     icnt,
     &     msc, maxexc, ip, ih, ivv, iv, ivv2,
     &     minh, maxh,
     &     minp, maxp, maxv, maxvv, minexc, maxcom, maxcom_en,
     &     n_t_cls, i_cls,
     &     n_tred_cls, len_form, optref, idef, ciroot
      logical ::
     &     pure_vv, update_prc, skip, ci_init
      character(len_target_name) ::
     &     dia_label, dia_label2,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character ::
     &     op_ht*3, f_ht*5, op_ht0to*6, f_ht0to*8, form_str*50,
     &     def_ht*10

      ! first set targets for CASSCF or uncontracted CI wave function
      ! (if not done already)
      if (.not.is_keyword_set('calculate.multiref').gt.0) then
        call set_ic_mrci_targets(tgt_info,orb_info)
      end if

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting multireference targets #3...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('calculate.multiref','minh',
     &     ival=minh)
      call get_argument_value('calculate.multiref','maxh',
     &     ival=maxh)
      call get_argument_value('calculate.multiref','minp',
     &     ival=minp)
      call get_argument_value('calculate.multiref','maxp',
     &     ival=maxp)
      call get_argument_value('calculate.multiref','maxv',
     &     ival=maxv)
      call get_argument_value('calculate.multiref','maxvv',
     &     ival=maxvv)
      call get_argument_value('calculate.multiref','minexc',
     &     ival=minexc)
      call get_argument_value('calculate.multiref','maxexc',
     &     ival=maxexc)
      if (maxh.lt.0) maxh = maxexc
      if (maxp.lt.0) maxp = maxexc
      if (maxv.lt.0) maxv = 2*maxexc
      if (maxvv.lt.0) maxvv = maxexc
      call get_argument_value('calculate.multiref','pure_vv',
     &     lval=pure_vv)
      call get_argument_value('calculate.multiref','ciroot',
     &     ival=ciroot)
      call get_argument_value('calculate.multiref','optref',
     &     ival=optref)
      call get_argument_value('calculate.multiref','update_prc',
     &     lval=update_prc)
      call get_argument_value('calculate.multiref','ci_init',
     &     lval=ci_init)
      if (.not.update_prc.and.optref.ne.-1)
     &      call warn('set_ic_mrcc_targets',
     &                'update_prc=F works only for optref=-1')
      call get_argument_value('method.MRCC','maxcom_res',
     &     ival=maxcom)
      call get_argument_value('method.MRCC','maxcom_en',
     &     ival=maxcom_en)
      
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving deexcitation operator L (for icMRCC)
      call add_target('L',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
cmh         no blocks which can be modeled by a contraction of two operators
c            if (ip.ge.1.and.ih.ge.1) cycle
cmh end
c            ! only hole-particle singles
c            if (max(ip,ih)+ivv.eq.1.and.(ip.ne.1.or.ih.ne.1)) cycle
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,2,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      n_t_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('L',ttype_op,DEF_OP_FROM_OCC,
     &              'L',1,1,
     &              parameters,2,tgt_info)

      ! define particle conserving excitation operator T (for icMRCC)
      call add_target('T',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
cmh         no blocks which can be modeled by a contraction of two operators
c            if (ip.ge.1.and.ih.ge.1) cycle
cmh end
c            ! only hole-particle singles
c            if (max(ip,ih)+ivv.eq.1.and.(ip.ne.1.or.ih.ne.1)) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      n_t_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('T',ttype_op,DEF_OP_FROM_OCC,
     &              'T',1,1,
     &              parameters,2,tgt_info)
cmh      ! T excitation operator
cmh      call add_target2('T',.false.,tgt_info)
      call set_dependency('T','L',tgt_info)
cmh      call set_rule2('T',CLONE_OP,tgt_info)
cmh      call set_arg('T',CLONE_OP,'LABEL',1,tgt_info,
cmh     &     val_label=(/'T'/))
cmh      call set_arg('T',CLONE_OP,'TEMPLATE',1,tgt_info,
cmh     &     val_label=(/'L'/))
cmh      call set_arg('T',CLONE_OP,'ADJOINT',1,tgt_info,
cmh     &     val_log=(/.true./))
c      call set_rule2('T',SET_ORDER,tgt_info)
c      call set_arg('T',SET_ORDER,'LABEL',1,tgt_info,
c     &     val_label=(/'T'/))
c      call set_arg('T',SET_ORDER,'SPECIES',1,tgt_info,
c     &             val_int=(/1/))

      ! TT intermediate
      call add_target('TT',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            if (ip.lt.1.or.ih.lt.1) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('TT',ttype_op,DEF_OP_FROM_OCC,
     &              'TT',1,1,
     &              parameters,2,tgt_info)

      ! transformed T
      call add_target2('Ttr',.false.,tgt_info)
      call set_dependency('Ttr','T',tgt_info)
      call set_rule2('Ttr',CLONE_OP,tgt_info)
      call set_arg('Ttr',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Ttr'/))
      call set_arg('Ttr',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'T'/))

      ! transformed L
      call add_target2('Ltr',.false.,tgt_info)
      call set_dependency('Ltr','L',tgt_info)
      call set_rule2('Ltr',CLONE_OP,tgt_info)
      call set_arg('Ltr',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Ltr'/))
      call set_arg('Ltr',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'L'/))

      ! subset of L for non-redundant valence-only metric
      call add_target('Lred',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = minp, maxp
        do ih = minh, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
c            ! exclude blocks with one or more hole-part. excitations
c            if (min(ih,ip).ge.1) cycle
            ! if same valence structure already exists, skip block
            skip = .false.
            do idef = 1, ndef
              skip = skip.or.(occ_def(IVALE,1,idef).eq.max(ip-ih,0)+ivv
     &                 .and.occ_def(IVALE,2,idef).eq.max(ih-ip,0)+ivv)
            end do
            if (skip) cycle
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,2,ndef) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      n_tred_cls = ndef
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('Lred',ttype_op,DEF_OP_FROM_OCC,
     &              'Lred',1,1,
     &              parameters,2,tgt_info)
      call set_rule2('Lred',SET_ORDER,tgt_info)
      call set_arg('Lred',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'Lred'/))
      call set_arg('Lred',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/1/))

      ! subset of T operator
      call add_target2('Tred',.false.,tgt_info)
      call set_dependency('Tred','Lred',tgt_info)
      call set_rule2('Tred',CLONE_OP,tgt_info)
      call set_arg('Tred',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Tred'/))
      call set_arg('Tred',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'Lred'/))
      call set_arg('Tred',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.true./))

      ! define Residual
      call add_target('OMG',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('OMG',ttype_op,DEF_OP_FROM_OCC,
     &              'OMG',1,1,
     &              parameters,2,tgt_info)

      ! define Jacobian (diagonal blocks only)
      call add_target('A(CC)',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).eq.0.and.(ivv.eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            if (abs(ih-ip).eq.0.and.ivv.eq.0) cycle ! no conv. blocks (for PREC)
            ndef = ndef + 1
            occ_def(IHOLE,1,3*ndef-1) = ih
            occ_def(IHOLE,2,3*ndef-1) = ih
            occ_def(IPART,1,3*ndef-1) = ip
            occ_def(IPART,2,3*ndef-1) = ip
            occ_def(IVALE,1,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-2) = max(ip-ih,0) + ivv
            occ_def(IVALE,1,3*ndef) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('A(CC)',ttype_op,DEF_OP_FROM_OCC,
     &              'A(CC)',1,1,
     &              parameters,2,tgt_info)

      ! Metric operator (for testing, also non-diagonal blocks)
      call add_target('S',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
           do ivv2 = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (abs(ih-ip)+2*max(ivv,ivv2).gt.maxv) cycle
            if (max(ip,ih).eq.0.and.
     &          (min(ivv,ivv2).eq.0.or..not.pure_vv)) cycle
            if (max(ip,ih)+min(ivv,ivv2).lt.minexc) cycle
            if (abs(ih-ip).eq.0.and.min(ivv,ivv2).eq.0) cycle ! no conv. blocks
            ndef = ndef + 1
            occ_def(IHOLE,1,3*ndef-1) = ih
            occ_def(IHOLE,2,3*ndef-1) = ih
            occ_def(IPART,1,3*ndef-1) = ip
            occ_def(IPART,2,3*ndef-1) = ip
            occ_def(IVALE,1,3*ndef-1) = max(ih-ip,0) + ivv2
            occ_def(IVALE,2,3*ndef-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,3*ndef-2) = max(ip-ih,0) + ivv2
            occ_def(IVALE,1,3*ndef) = max(ip-ih,0) + ivv
           end do
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('S',ttype_op,DEF_OP_FROM_OCC,
     &              'S',1,1,
     &              parameters,2,tgt_info)

      ! Diagonal Preconditioner
      call add_target(op_dia//'_T',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_T','T',tgt_info)
      call cloneop_parameters(-1,parameters,'T',.false.)
      call set_rule(op_dia//'_T',ttype_op,CLONE_OP,op_dia//'_T',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner for L
      call add_target(op_dia//'_L',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_L','L',tgt_info)
      call cloneop_parameters(-1,parameters,'L',.false.)
      call set_rule(op_dia//'_L',ttype_op,CLONE_OP,op_dia//'_L',1,1,
     &              parameters,1,tgt_info)

      ! C0 dagger operator (for testing)
      call add_target2('C0dag',.false.,tgt_info)
      call set_dependency('C0dag','C0',tgt_info)
      call set_rule2('C0dag',CLONE_OP,tgt_info)
      call set_arg('C0dag',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'C0dag'/))
      call set_arg('C0dag',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'C0'/))
      call set_arg('C0dag',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.true./))

      ! HT commutator: blocks up to rank 2
      op_ht = 'HT '
      do icnt = 1, 6
        write(op_ht(3:3),'(i1)') icnt
        call add_target2(op_ht,.false.,tgt_info)
        call set_dependency(op_ht,op_ham,tgt_info)
        call set_rule2(op_ht,CLONE_OP,tgt_info)
        call set_arg(op_ht,CLONE_OP,'LABEL',1,tgt_info,
     &       val_label=(/op_ht/))
        call set_arg(op_ht,CLONE_OP,'TEMPLATE',1,tgt_info,
     &       val_label=(/op_ham/))
      end do

      ! sums of HT commutators
      op_ht0to = 'HT0to '
      do icnt = 1, 4
        write(op_ht0to(6:6),'(i1)') icnt
        call add_target2(op_ht0to,.false.,tgt_info)
        call set_dependency(op_ht0to,op_ham,tgt_info)
        call set_rule2(op_ht0to,CLONE_OP,tgt_info)
        call set_arg(op_ht0to,CLONE_OP,'LABEL',1,tgt_info,
     &       val_label=(/op_ht0to/))
        call set_arg(op_ht0to,CLONE_OP,'TEMPLATE',1,tgt_info,
     &       val_label=(/op_ham/))
      end do
*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! multireference CC lagrangian
      ! a) set up
      call add_target2('F_MRCC_LAG',.false.,tgt_info)
      call set_dependency('F_MRCC_LAG','NORM',tgt_info)
      call set_dependency('F_MRCC_LAG','H',tgt_info)
      call set_dependency('F_MRCC_LAG','C0',tgt_info)
      call set_dependency('F_MRCC_LAG','T',tgt_info)
      call set_dependency('F_MRCC_LAG','L',tgt_info)
      call set_rule2('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &     tgt_info,val_label=(/'L','H','T','C0'/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/maxcom/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='---')
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &     tgt_info,val_str='ic-MRCC Lagrangian')
      if (optref.eq.-1) then
        call set_dependency('F_MRCC_LAG','E(MR)',tgt_info)
        call set_rule2('F_MRCC_LAG',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'E(MR)','C0^+','C0'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'IDX_SV',3,
     &       tgt_info,
     &       val_int=(/2,3,4/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &       val_rl8=(/-1d0/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end if
      if (.false..and.maxh.gt.0) then
        ! replace T-T double contractions by TT intermediate
        call set_dependency('F_MRCC_LAG','F_TT',tgt_info)
        call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
      end if
      ! factor out HT intermediates
c      call set_dependency('F_MRCC_LAG','F_HT1',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT2',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT0to1',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT0to2',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT3',tgt_info)
c      call set_dependency('F_MRCC_LAG','F_HT4',tgt_info)
c      call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',4,tgt_info,
c     &     val_label=(/'F_HT2','F_HT1','F_HT0to2','F_HT0to1'/))
      call set_rule2('F_MRCC_LAG',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_LAG',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))

      ! Residual
      call add_target2('F_OMG',.false.,tgt_info)
      call set_dependency('F_OMG','F_MRCC_LAG',tgt_info)
      call set_dependency('F_OMG','OMG',tgt_info)
      call set_rule2('F_OMG',DERIVATIVE,tgt_info)
      call set_arg('F_OMG',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      call set_arg('F_OMG',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_OMG',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_OMG',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))
c dbg
c      call set_rule2('F_OMG',PRINT_FORMULA,tgt_info)
c      call set_arg('F_OMG',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG'/))
c dbgend

      ! Residual for C0
      call add_target2('F_OMG_C0',.false.,tgt_info)
      call set_dependency('F_OMG_C0','F_MRCC_LAG',tgt_info)
      call set_dependency('F_OMG_C0','A_C0',tgt_info)
      call set_rule2('F_OMG_C0',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'C0^+'/))
      call set_rule2('F_OMG_C0',INVARIANT,tgt_info)
      call set_arg('F_OMG_C0',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C0'/))
      call set_arg('F_OMG_C0',INVARIANT,'OPERATORS',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_OMG_C0',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Residual for Reference function')
c dbg
c      call set_rule2('F_OMG_C0',PRINT_FORMULA,tgt_info)
c      call set_arg('F_OMG_C0',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c dbgend

      ! Energy
      call add_target2('F_MRCC_E',.false.,tgt_info)
      call set_dependency('F_MRCC_E','F_MRCC_LAG',tgt_info)
      call set_dependency('F_MRCC_E','E(MR)',tgt_info)
      call set_dependency('F_MRCC_E','L',tgt_info)
      call set_rule2('F_MRCC_E',INVARIANT,tgt_info)
      call set_arg('F_MRCC_E',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))
      call set_arg('F_MRCC_E',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_E',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',2,tgt_info,
     &     val_label=(/'L','E(MR)'/))
      call set_arg('F_MRCC_E',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC energy expression')
      call set_rule2('F_MRCC_E',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_E',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))

      ! multireference CC norm
      ! a) set up using reduced densities
      call add_target2('F_MRCC_NORM',.false.,tgt_info)
      call set_dependency('F_MRCC_NORM','NORM',tgt_info)
      call set_dependency('F_MRCC_NORM','Tred',tgt_info)
      call set_dependency('F_MRCC_NORM','Lred',tgt_info)
      call set_dependency('F_MRCC_NORM','DENS',tgt_info)
      call set_rule2('F_MRCC_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &     tgt_info,
     &     val_label=(/'Lred','Tred'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &     val_int=(/2,3/))
      call set_rule2('F_MRCC_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'DENS','Lred','Tred','DENS'/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,2/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &     val_int=(/1,4/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &     val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
      call set_arg('F_MRCC_NORM',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! b) insert unit operators to allow for differentiation
      call set_dependency('F_MRCC_NORM','1v',tgt_info)
      call set_rule2('F_MRCC_NORM',INSERT,tgt_info)
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1v'/))
      call set_arg('F_MRCC_NORM',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'Lred','Tred'/))
      ! c) replace 1v by 1 (was used because we only needed valence blocks)
      call set_dependency('F_MRCC_NORM','1',tgt_info)
      call set_rule2('F_MRCC_NORM',REPLACE,tgt_info)
      call set_arg('F_MRCC_NORM',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_NORM',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1v','1'/))
c dbg
c      call set_rule2('F_MRCC_NORM',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_NORM'/))
c dbgend

      ! transformed excitation operator
      call add_target2('F_T',.false.,tgt_info)
      call set_dependency('F_T','T',tgt_info)
      call set_dependency('F_T','Dtr',tgt_info)
      call set_dependency('F_T','Ttr',tgt_info)
      do i_cls = 1, n_t_cls
        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'T','Dtr','Ttr','Dtr','T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &       val_int=(/i_cls,1,i_cls,1,i_cls/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',5,tgt_info,
     &       val_int=(/i_cls,-1,i_cls,-1,i_cls/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/2/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'AVOID',4,tgt_info,
     &       val_int=(/3,5,2,4/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_T'/))
      call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'T'/))
      call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'T','Ttr','T'/))
      call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! no active open lines from Ttr
      call set_rule2('F_T',SELECT_LINE,tgt_info)
      call set_arg('F_T',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_T'/))
      call set_arg('F_T',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_T'/))
      call set_arg('F_T',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'T'/))
      call set_arg('F_T',SELECT_LINE,'OP_INCL',1,tgt_info,
     &     val_label=(/'Ttr'/))
      call set_arg('F_T',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_T',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='no_ext')
c dbg
c      call set_rule2('F_T',PRINT_FORMULA,tgt_info)
c      call set_arg('F_T',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_T'/))
c dbgend

      ! transformed deexcitation operator
      call add_target2('F_L',.false.,tgt_info)
      call set_dependency('F_L','L',tgt_info)
      call set_dependency('F_L','Dtr',tgt_info)
      call set_dependency('F_L','Ltr',tgt_info)
      do i_cls = 1, n_t_cls
        call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'L','Dtr^+','Ltr','Dtr^+','L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MIN',5,tgt_info,
     &       val_int=(/i_cls,1,i_cls,1,i_cls/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MAX',5,tgt_info,
     &       val_int=(/i_cls,-1,i_cls,-1,i_cls/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/2/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'AVOID',4,tgt_info,
     &       val_int=(/1,3,2,4/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_L'/))
      call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &     tgt_info,
     &     val_label=(/'L','Ltr','L'/))
      call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &     val_int=(/1,2,1/))
      call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! no active open lines from Ttr
      call set_rule2('F_L',SELECT_LINE,tgt_info)
      call set_arg('F_L',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_L'/))
      call set_arg('F_L',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_L'/))
      call set_arg('F_L',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_L',SELECT_LINE,'OP_INCL',1,tgt_info,
     &     val_label=(/'Ltr'/))
      call set_arg('F_L',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_L',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='no_ext')
c dbg
c      call set_rule2('F_L',PRINT_FORMULA,tgt_info)
c      call set_arg('F_L',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_L'/))
c dbgend

      ! transformed multireference CC energy
      ! a) set up
      call add_target2('F_E(MRCC)tr',.false.,tgt_info)
      call set_dependency('F_E(MRCC)tr','E(MR)',tgt_info)
      call set_dependency('F_E(MRCC)tr','H',tgt_info)
c      call set_dependency('F_E(MRCC)tr','FREF',tgt_info)
      call set_dependency('F_E(MRCC)tr','C0',tgt_info)
      call set_dependency('F_E(MRCC)tr','T',tgt_info)
      call set_dependency('F_E(MRCC)tr','L',tgt_info)
      call set_rule2('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &     tgt_info,val_label=(/'E(MR)'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &     tgt_info,val_label=(/'L','H','T','C0'/))
c     &     tgt_info,val_label=(/'L','FREF','T','C0'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/0/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='NOSCAL')
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='Precursor for linearized Jacobian')
      ! f) insert 1 (particle/hole space) for later differentiation
      call set_dependency('F_E(MRCC)tr','1ph',tgt_info)
      call set_rule2('F_E(MRCC)tr',INSERT,tgt_info)
      call set_arg('F_E(MRCC)tr',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1ph'/))
      call set_arg('F_E(MRCC)tr',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'L','T'/))
      ! replace 1ph by 1
      call set_dependency('F_E(MRCC)tr','1',tgt_info)
      call set_rule2('F_E(MRCC)tr',REPLACE,tgt_info)
      call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1ph','1'/))
      ! g) expand T and L
      call set_dependency('F_E(MRCC)tr','F_T',tgt_info)
      call set_dependency('F_E(MRCC)tr','F_L',tgt_info)
      call set_rule2('F_E(MRCC)tr',EXPAND,tgt_info)
      call set_arg('F_E(MRCC)tr',EXPAND,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',EXPAND,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',EXPAND,'INTERM',2,tgt_info,
     &     val_label=(/'F_T','F_L'/))
c dbg
c      call set_rule2('F_E(MRCC)tr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_E(MRCC)tr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_E(MRCC)tr'/))
c dbgend

      ! (transformed) Jacobian times vector
      call add_target2('F_A_Ttr',.false.,tgt_info)
      call set_dependency('F_A_Ttr','F_E(MRCC)tr',tgt_info)
      call set_dependency('F_A_Ttr','OMG',tgt_info)
      call set_rule2('F_A_Ttr',DERIVATIVE,tgt_info)
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Ltr'/))
c dbg
c      call set_rule2('F_A_Ttr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_A_Ttr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_A_Ttr'/))
c dbgend

      ! (transformed) Jacobian
      call add_target2('F_Atr',.false.,tgt_info)
      call set_dependency('F_Atr','F_A_Ttr',tgt_info)
      call set_dependency('F_Atr','A(CC)',tgt_info)
      call set_rule2('F_Atr',DERIVATIVE,tgt_info)
      call set_arg('F_Atr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Atr'/))
      call set_arg('F_Atr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_Atr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A(CC)'/))
      call set_arg('F_Atr',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Ttr'/))
c dbg
c      call set_rule2('F_Atr',KEEP_TERMS,tgt_info)
c      call set_arg('F_Atr',KEEP_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
c      call set_arg('F_Atr',KEEP_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
c      call set_arg('F_Atr',KEEP_TERMS,'TERMS',2,tgt_info,
c     &     val_int=(/1,7/))
c
c      call set_rule2('F_Atr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Atr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
c dbgend

      ! Metric times amplitudes
      call add_target2('F_MRCC_SC',.false.,tgt_info)
      call set_dependency('F_MRCC_SC','F_MRCC_NORM',tgt_info)
      call set_dependency('F_MRCC_SC','OMG',tgt_info)
      call set_rule2('F_MRCC_SC',DERIVATIVE,tgt_info)
      call set_arg('F_MRCC_SC',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_SC'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_NORM'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_MRCC_SC',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Lred'/))
c dbg
c      call set_rule2('F_MRCC_SC',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_SC',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_SC'/))
c dbgend

      ! Metric (valence-only part)
      call add_target2('F_MRCC_D',.false.,tgt_info)
      call set_dependency('F_MRCC_D','F_MRCC_SC',tgt_info)
      call set_dependency('F_MRCC_D','D',tgt_info)
      call set_rule2('F_MRCC_D',DERIVATIVE,tgt_info)
      call set_arg('F_MRCC_D',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_D'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_SC'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'D'/))
      call set_arg('F_MRCC_D',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'Tred'/))
c dbg
      call set_rule2('F_MRCC_D',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_D',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_D'/))
c dbgend

      ! HT intermediate
      f_ht = 'F_HT '
      op_ht = 'HT '
      do icnt = 1, 2
        write(f_ht(5:5),'(i1)') icnt
        write(op_ht(3:3),'(i1)') icnt
        call add_target2(f_ht,.false.,tgt_info)
        call set_dependency(f_ht,op_ht,tgt_info)
        call set_dependency(f_ht,op_ham,tgt_info)
        call set_dependency(f_ht,'T',tgt_info)
        call set_rule2(f_ht,DEF_MRCC_LAGRANGIAN,tgt_info)
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &       val_label=(/f_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'OP_RES',1,tgt_info,
     &       val_label=(/op_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'OPERATORS',4,tgt_info,
     &       val_label=(/'T  ','H  ','T  ',op_ht/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,tgt_info,
     &       val_int=(/icnt/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,tgt_info,
     &       val_int=(/icnt/))
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &       val_str='HBAR')
        call set_arg(f_ht,DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &       val_str='Commutator intermediate')
        if (icnt.gt.1.and.maxh.gt.0) then
          ! replace T-T double contractions by TT intermediate
          call set_dependency(f_ht,'F_TT',tgt_info)
          call set_rule2(f_ht,FACTOR_OUT,tgt_info)
          call set_arg(f_ht,FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &         val_label=(/f_ht/))
          call set_arg(f_ht,FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &         val_label=(/f_ht/))
          call set_arg(f_ht,FACTOR_OUT,'INTERM',1,tgt_info,
     &         val_label=(/'F_TT'/))
        end if
c dbg
c        if (icnt.eq.2) then
c        ! factor out lower HT intermediates
c        call set_dependency(f_ht,'F_HT1',tgt_info)
c        call set_rule2(f_ht,FACTOR_OUT,tgt_info)
c        call set_arg(f_ht,FACTOR_OUT,'LABEL_RES',1,tgt_info,
c     &       val_label=(/f_ht/))
c        call set_arg(f_ht,FACTOR_OUT,'LABEL_IN',1,tgt_info,
c     &       val_label=(/f_ht/))
c        call set_arg(f_ht,FACTOR_OUT,'INTERM',1,tgt_info,
c     &       val_label=(/'F_HT1'/))
c        end if
c dbgend
        call set_rule2(f_ht,PRINT_FORMULA,tgt_info)
        call set_arg(f_ht,PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht/))
      end do

      ! sums of HT intermediates
      f_ht0to = 'F_HT0to '
      op_ht0to = 'HT0to '
      form_str = 'HT0to =H'
      len_form = 8
      do icnt = 1, 4
        write(f_ht0to(8:8),'(i1)') icnt
        write(op_ht0to(6:6),'(i1)') icnt
        call add_target2(f_ht0to,.false.,tgt_info)
        call set_dependency(f_ht0to,op_ht0to,tgt_info)
        call set_dependency(f_ht0to,op_ham,tgt_info)
        op_ht = 'HT '
        do ih = 1, icnt
          write(op_ht(3:3),'(i1)') ih
          call set_dependency(f_ht0to,op_ht,tgt_info)
        end do
        write(form_str(6:6),'(i1)') icnt
        write(form_str(len_form+1:len_form+4),'(a1,a3)') '+',op_ht
        len_form = len_form + 4
        call set_rule2(f_ht0to,DEF_FORMULA,tgt_info)
        call set_arg(f_ht0to,DEF_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht0to/))
        call set_arg(f_ht0to,DEF_FORMULA,'FORMULA',1,tgt_info,
     &       val_str=form_str(1:len_form))
        call set_rule2(f_ht0to,PRINT_FORMULA,tgt_info)
        call set_arg(f_ht0to,PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/f_ht0to/))
      end do

      call add_target2('F_TT',.false.,tgt_info)
      call set_dependency('F_TT','TT',tgt_info)
      call set_dependency('F_TT','T',tgt_info)
      call set_rule2('F_TT',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'TT','T','T','TT'/))
      call set_arg('F_TT',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_rule2('F_TT',PRINT_FORMULA,tgt_info)
      call set_arg('F_TT',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_TT'/))

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! density matrix
      call add_target2('FOPT_MRCC_D',.false.,tgt_info)
      call set_dependency('FOPT_MRCC_D','F_MRCC_D',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_D',tgt_info)
      call set_dependency('FOPT_MRCC_D','DEF_ME_DENS',tgt_info)
      call set_rule2('FOPT_MRCC_D',OPTIMIZE,tgt_info)
      call set_arg('FOPT_MRCC_D',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_MRCC_D'/))
      call set_arg('FOPT_MRCC_D',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_MRCC_D'/))

      ! transformed Hessian
      call add_target2('FOPT_Atr',.false.,tgt_info)
      call set_dependency('FOPT_Atr','F_Atr',tgt_info)
      call set_dependency('FOPT_Atr',mel_ham,tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_A(CC)',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_Dtr',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_Atr',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Atr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Atr'/))
      call set_arg('FOPT_Atr',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Atr'/))

      ! Residual
      call add_target2('FOPT_OMG',.false.,tgt_info)
      call set_dependency('FOPT_OMG','F_OMG',tgt_info)
      call set_dependency('FOPT_OMG','F_MRCC_E',tgt_info)
      call set_dependency('FOPT_OMG','F_T',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_Dtr',tgt_info)
      if (.false..and.maxh.gt.0)
     &    call set_dependency('FOPT_OMG','DEF_ME_TT',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT1',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT2',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_OMG',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG',mel_ham,tgt_info)
      if (optref.ne.0) then
        call set_dependency('FOPT_OMG','F_DENS0',tgt_info)
        call set_dependency('FOPT_OMG','F_MRCC_D',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_DENS',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_D',tgt_info)
        call set_rule2('FOPT_OMG',ASSIGN_ME2OP,tgt_info)
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_D'/))
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &             val_label=(/'D'/))
        call set_dependency('FOPT_OMG','F_Atr',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_A(CC)',tgt_info)
c        call set_dependency('FOPT_OMG','DEF_ME_FREF',tgt_info)
        if (optref.eq.-1) then
          call set_dependency('FOPT_OMG','F_OMG_C0',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_A_C0',tgt_info)
          call set_rule2('FOPT_OMG',ASSIGN_ME2OP,tgt_info)
          call set_arg('FOPT_OMG',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &               val_label=(/'ME_A_C0'/))
          call set_arg('FOPT_OMG',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &               val_label=(/'A_C0'/))
        end if
      end if
      call set_rule2('FOPT_OMG',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG'/))
      if (optref.eq.-1.and.update_prc) then
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',7,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG','F_OMG_C0','F_T',
     &                         'F_DENS0','F_MRCC_D','F_Atr'/))
      else if (optref.eq.-1.and..not.update_prc) then
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',6,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG','F_OMG_C0','F_T',
     &                         'F_DENS0','F_MRCC_D'/))
      else if (optref.ne.0) then
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',6,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG','F_T',
     &                         'F_DENS0','F_MRCC_D','F_Atr'/))
      else
        call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',3,tgt_info,
     &             val_label=(/'F_MRCC_E','F_OMG','F_T'/))
      end if

      ! Residual for C0
      call add_target2('FOPT_OMG_C0',.false.,tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_T',tgt_info)
      if (.false..and.maxh.gt.0)
     &    call set_dependency('FOPT_OMG_C0','DEF_ME_TT',tgt_info)
c      call set_dependency('FOPT_OMG_C0','DEF_ME_HT1',tgt_info)
c      call set_dependency('FOPT_OMG_C0','DEF_ME_HT2',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG_C0',mel_ham,tgt_info)
      call set_dependency('FOPT_OMG_C0','F_OMG_C0',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_A_C0',tgt_info)
      call set_rule2('FOPT_OMG_C0',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_OMG_C0',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_A_C0'/))
      call set_arg('FOPT_OMG_C0',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'A_C0'/))
      call set_rule2('FOPT_OMG_C0',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG_C0',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG_C0'/))
      call set_arg('FOPT_OMG_C0',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_OMG_C0'/))
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! ME for Hessian
      call add_target2('DEF_ME_A(CC)',.false.,tgt_info)
      call set_dependency('DEF_ME_A(CC)','A(CC)',tgt_info)
      call set_rule2('DEF_ME_A(CC)',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_A(CC)'/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'A(CC)'/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A(CC)',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &     val_int=(/0/))

      ! Diagonal Preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call add_target(trim(dia_label),ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label),'EVAL_FREF',tgt_info)
      call set_dependency(trim(dia_label),
     &                    op_dia//'_'//'T',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)
      labels(2) = op_dia//'_'//'T'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule(trim(dia_label),ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      ! use effective Fock op. (needed only for pure inactive exc.)
      labels(1) = trim(dia_label)
      labels(2) = 'ME_FREF'
      call set_rule(trim(dia_label),ttype_opme,
     &              PRECONDITIONER,
     &              labels,2,1,
     &              parameters,1,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (a):',0,'LIST')
c      call set_rule(trim(dia_label),ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME for T
      call add_target2('DEF_ME_T',.false.,tgt_info)
      call set_dependency('DEF_ME_T','T',tgt_info)
      call set_rule2('DEF_ME_T',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_T',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_T'/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'T'/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for Ttr
      call add_target2('DEF_ME_Ttr',.false.,tgt_info)
      call set_dependency('DEF_ME_Ttr','Ttr',tgt_info)
      call set_rule2('DEF_ME_Ttr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Ttr'/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! ME-List for HT intermediates
      op_ht = 'HT '
      op_ht0to = 'ME_HT '
      def_ht = 'DEF_ME_HT '
      do icnt = 1, 4
        write(op_ht(3:3),'(i1)') icnt
        write(op_ht0to(6:6),'(i1)') icnt
        write(def_ht(10:10),'(i1)') icnt
        call add_target2(def_ht,.false.,tgt_info)
        call set_dependency(def_ht,op_ht,tgt_info)
        call set_rule2(def_ht,DEF_ME_LIST,tgt_info)
        call set_arg(def_ht,DEF_ME_LIST,'LIST',1,tgt_info,
     &               val_label=(/op_ht0to/))
        call set_arg(def_ht,DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &               val_label=(/op_ht/))
        call set_arg(def_ht,DEF_ME_LIST,'MS',1,tgt_info,
     &               val_int=(/0/))
        call set_arg(def_ht,DEF_ME_LIST,'IRREP',1,tgt_info,
     &               val_int=(/1/))
      end do

      ! ME for TT intermediate
      call add_target2('DEF_ME_TT',.false.,tgt_info)
      call set_dependency('DEF_ME_TT','TT',tgt_info)
      call set_rule2('DEF_ME_TT',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_TT'/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'TT'/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! ME for Residual
      call add_target2('DEF_ME_OMG',.false.,tgt_info)
      call set_dependency('DEF_ME_OMG','OMG',tgt_info)
      call set_rule2('DEF_ME_OMG',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMG'/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'OMG'/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! reordered projector matrix ME_Dproj (to eliminate lin. dep.)
      call add_target('DEF_ME_Dproj',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dproj','Dtr',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dproj'
      labels(2) = 'Dtr'
      call me_list_parameters(-1,parameters,
     &     msc,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dproj',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_dependency('DEF_ME_Dproj','DEF_ME_Dinv',tgt_info)
      labels(1) = 'ME_Dproj'
      labels(2) = 'ME_D'
      call form_parameters(-1,parameters,2,
     &     '---',13,'---')
      call set_rule('DEF_ME_Dproj',ttype_opme,
     &              REORDER_MEL,
     &              labels,2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Reordered projector matrix :',0,'LIST')
c      call set_rule('DEF_ME_Dproj',ttype_opme,PRINT_MEL,
c     &     'ME_Dproj',1,0,
c     &     parameters,2,tgt_info)
c dbgend

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! evaluate valence-only metric
      call add_target2('EVAL_MRCC_D',.false.,tgt_info)
      call set_dependency('EVAL_MRCC_D','FOPT_MRCC_D',tgt_info)
      call set_dependency('EVAL_MRCC_D','EVAL_DENS0',tgt_info)
      call set_rule2('EVAL_MRCC_D',EVAL,tgt_info)
      call set_arg('EVAL_MRCC_D',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_MRCC_D'/))
c dbg
c      call set_rule2('EVAL_MRCC_D',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_MRCC_D',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_D'/))
c dbgend

      ! Evaluate diagonal elements of Jacobian
      call add_target('EVAL_Atr',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_Atr','FOPT_Atr',tgt_info)
c      call set_dependency('EVAL_Atr','EVAL_FREF',tgt_info)
      call set_rule('EVAL_Atr',ttype_opme,EVAL,
     &     'FOPT_Atr',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'transformed Jacobian :',0,'LIST')
c      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
c     &     'ME_A(CC)',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      ! put diagonal elements to preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      labels(1) = trim(dia_label)
      labels(2) = 'ME_A(CC)'
      call set_dependency('EVAL_Atr',trim(dia_label),tgt_info)
      call set_rule('EVAL_Atr',ttype_opme,
     &              EXTRACT_DIAG,
     &              labels,2,1,
     &              parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (b) :',0,'LIST')
c      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Solve MR coupled cluster equations
      call add_target2('SOLVE_MRCC',.true.,tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_OMG',tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_REF_S(S+1)',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('SOLVE_MRCC',trim(dia_label),tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_Atr',tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_MRCC_D',tgt_info)
      call set_dependency('SOLVE_MRCC','DEF_ME_Dtrdag',tgt_info)
      if (optref.ne.0) then
        call me_list_label(dia_label2,mel_dia,orb_info%lsym,
     &                     0,0,0,.false.)
        dia_label2 = trim(dia_label2)//'C0'
        call set_dependency('SOLVE_MRCC',trim(dia_label2),tgt_info)
        if (optref.ne.-1) 
     &     call set_dependency('SOLVE_MRCC','FOPT_OMG_C0',tgt_info)
        call set_dependency('SOLVE_MRCC','DEF_ME_Dproj',tgt_info)
      end if
      do icnt = 1, max(1,optref)
      call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
      if (optref.eq.-1) then
        if (.not.ci_init) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &         val_label=(/'ME_T'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &         val_str='TRF')
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &         val_label=(/'ME_OMG'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &         val_label=(/trim(dia_label)/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &       val_label=(/'ME_E(MR)'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &       val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &         val_label=(/'FOPT_OMG'/))
          call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
        end if
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',2,tgt_info,
     &       val_label=(/'ME_T','ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF/NRM')
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',2,tgt_info,
     &       val_label=(/'ME_OMG','ME_A_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',2,tgt_info,
     &       val_label=(/trim(dia_label),trim(dia_label2)/))
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_T'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str='TRF')
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &       val_label=(/'ME_OMG'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label)/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &     val_label=(/'ME_E(MR)'/))
      if (optref.eq.-2) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',8,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag','ME_Dproj',
     &                 'ME_D','ME_Dinv',
     &                 'ME_A(CC)','ME_C0'/))
      else if (optref.ne.0.and.update_prc) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag','ME_Dproj',
     &                 'ME_D','ME_Dinv',
     &                 'ME_A(CC)'/))
      else if (optref.ne.0.and..not.update_prc) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',6,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag','ME_Dproj',
     &                 'ME_D','ME_Dinv'/))
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &     val_label=(/'ME_Ttr','ME_Dtr','ME_Dtrdag'/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &     val_label=(/'FOPT_OMG'/))

      if (optref.gt.0.and.icnt.ne.optref) then !not in last iteration
        call set_rule2('SOLVE_MRCC',SOLVEEVP,tgt_info)
        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'MODE',1,tgt_info,
     &       val_str='DIA')
        call set_arg('SOLVE_MRCC',SOLVEEVP,'N_ROOTS',1,tgt_info,
     &       val_int=(/ciroot/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_MVP',1,tgt_info,
     &       val_label=(/'A_C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label2)/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'OP_SVP',1,tgt_info,
     &     val_label=(/'C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'FORM',1,tgt_info,
     &       val_label=(/'FOPT_OMG_C0'/))
c dbg
        call form_parameters(-1,parameters,2,
     &       'CI coefficients :',0,'LIST')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
c dbgend
      end if
      end do
c dbg
      if (optref.ne.0) then
        call form_parameters(-1,parameters,2,
     &       'final CI coefficients :',0,'LIST')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_C0',1,0,
     &       parameters,2,tgt_info)
      end if
c dbgend

c dbg
c        call form_parameters(-1,parameters,2,
c     &       'final T amplitudes :',0,'LIST')
c        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &       'ME_T',1,0,
c     &       parameters,2,tgt_info)
c dbgend

      return
      end
