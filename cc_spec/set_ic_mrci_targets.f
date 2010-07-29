*----------------------------------------------------------------------*
      subroutine set_ic_mrci_targets(tgt_info,orb_info)
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

      integer ::
     &     min_rank, max_rank, ndef, occ_def(ngastp,2,60),
     &     isim, ncat, nint, icnt,
     &     isym, ms, msc, sym_arr(8), maxexc, ip, ih, ivv, iv,
     &     cminh, cmaxh, cminp, cmaxp, cmaxexc, minh, maxh,
     &     minp, maxp, maxv, maxvv, minexc, cbarc(2),
     &     nlabels, nroots, gno
      logical ::
     &     use_hessian, use_dens, pure_vv, calc
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character*4 ::
     &     op_exc, op_deexc

      ! get maximum excitation rank
      call get_argument_value('calculate.multiref','maxexc',
     &     ival=maxexc)
      call get_argument_value('calculate.multiref','cmaxexc',
     &     ival=cmaxexc)

      ! first set targets for CASSCF or uncontracted CI wave function
      call set_unc_mrci_targets(tgt_info,orb_info,
     &                          .false..or.maxexc.eq.0)

      ! if maxexc = 0: return because call of unc_mrci is sufficient
      if (maxexc.eq.0) then
        return
      else if (cmaxexc.gt.0) then
        call quit(1,'set_ic_mrci_targets',
     &            'Warning: Only tested for CASSCF reference so far')
      end if

      if (iprlvl.gt.0)
     &     write(luout,*) 'setting multireference targets #2...'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1 ! assuming closed shell

      ! which normal ordering is used?
      call get_argument_value('calculate.multiref','GNO',
     &     ival=gno)
      select case(gno)
      case(0)
      case(1)
        write(luout,*) 'Using generalized normal order (GNO)'
c        call set_gno_targets(tgt_info,orb_info,1)
      case default
        call quit(1,'set_ic_mrci_targets','unknown normal order')
      end select
cmh
      call set_gno_targets(tgt_info,orb_info,1)
cmh end
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
      ! minimum excitation: maxexc for GNO=0
      if (gno.eq.0) minexc = maxexc
      if (maxh.lt.0) maxh = maxexc
      if (maxp.lt.0) maxp = maxexc
      if (maxv.lt.0) maxv = 2*maxexc
      if (maxvv.lt.0) maxvv = maxexc
      call get_argument_value('calculate.multiref','pure_vv',
     &     lval=pure_vv)
      call get_argument_value('calculate.multiref','calc',
     &     lval=calc)
      call get_argument_value('calculate.multiref','nroots',
     &     ival=nroots)

      if (ntest.ge.100) then
        print *,'minh    = ',minh
        print *,'maxh    = ',maxh
        print *,'minp    = ',minp
        print *,'maxp    = ',maxp
        print *,'maxv    = ',maxv
        print *,'maxvv   = ',maxvv
        print *,'minexc  = ',minexc
        print *,'maxexc  = ',maxexc
        print *,'pure_vv = ',pure_vv
        print *,'nroots  = ',nroots
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define scalar multireference energy
      call add_target('E(MR)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('E(MR)',ttype_op,DEF_HAMILTONIAN,'E(MR)',
     &              1,1,parameters,1,tgt_info)

      ! define particle conserving excitation operator C (for icMRCI)
      call add_target('C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (.not.pure_vv.and.ivv.gt.0.and.ip.eq.0.and.ih.eq.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if ((max(ip,ih).gt.0.or.pure_vv).and.
     &          max(ip,ih)+ivv.lt.minexc) cycle
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
            do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
              if (.not.pure_vv.and.ivv.gt.0.and.ip.eq.0.and.ih.eq.0)
     &            cycle
              if (abs(ih-ip)+2*ivv.gt.maxv) cycle
              if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
              if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
              if ((max(ip,ih).gt.0.or.pure_vv).and.
     &            max(ip,ih)+ivv.lt.minexc) cycle
              ! exclude blocks with one or more hole-part. excitations
              if (min(ih,ip).ge.1) cycle
              if (max(ih,ip).eq.0.and.ivv.eq.0) cycle
              ndef = ndef + 1
              occ_def(IHOLE,2,ndef) = ih
              occ_def(IPART,1,ndef) = ip
              occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
              occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
            end do
          end do
        end do
        call op_from_occ_parameters(-1,parameters,2,
     &                occ_def,ndef,1,(/0,0/),ndef)
        call set_rule('c',ttype_op,DEF_OP_FROM_OCC,
     &                'c',1,1,
     &                parameters,2,tgt_info)
      end if

      ! define unit operator suited to connect C^+ and C
      call add_target('1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
           do iv = 0, maxexc
            ! not needed for pure inactive excitations
            if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
            ! only pure active or pure non-active occ. classes needed for now
            if (iv.gt.0.and.ih+ip.gt.0) cycle
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,1,ndef) = iv
            occ_def(IVALE,2,ndef) = iv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1',ttype_op,DEF_OP_FROM_OCC,
     &              '1',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1',ttype_op,SET_HERMIT,
     &              '1',1,1,
     &              parameters,1,tgt_info)

      ! define unit operator (scalar and non-valence blocks only)
      call add_target('1ph',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          ! not needed for pure inactive excitations
c          if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
          if (ip.ge.2.and.ih.ge.2) cycle
          ndef = ndef + 1
          occ_def(IHOLE,1,ndef) = ih
          occ_def(IHOLE,2,ndef) = ih
          occ_def(IPART,1,ndef) = ip
          occ_def(IPART,2,ndef) = ip
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1ph',ttype_op,DEF_OP_FROM_OCC,
     &              '1ph',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1ph',ttype_op,SET_HERMIT,
     &              '1ph',1,1,
     &              parameters,1,tgt_info)

      ! define unit operator (valence-only blocks)
      call add_target('1v',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 1, maxexc
        ndef = ndef + 1
        occ_def(IVALE,1,ndef) = iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1v',ttype_op,DEF_OP_FROM_OCC,
     &              '1v',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1v',ttype_op,SET_HERMIT,
     &              '1v',1,1,
     &              parameters,1,tgt_info)

      ! define scalar norm
      call add_target('NORM',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('NORM',ttype_op,DEF_HAMILTONIAN,'NORM',
     &              1,1,parameters,1,tgt_info)

      ! define density matrix
      call add_target('D',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
            if (.not.pure_vv.and.ip.eq.0.and.ih.eq.0.and.ivv.gt.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*3-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*3-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef*3) = max(ip-ih,0) + ivv
            occ_def(IVALE,2,ndef*3-2) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      use_dens = ndef.gt.1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('D',ttype_op,DEF_OP_FROM_OCC,
     &              'D',1,1,
     &              parameters,2,tgt_info)

      ! define daggered density matrix
      call add_target('Ddag',ttype_op,.false.,tgt_info)
      call set_dependency('Ddag','D',tgt_info)
      call cloneop_parameters(-1,parameters,'D',.true.)
      call set_rule('Ddag',ttype_op,CLONE_OP,'Ddag',1,1,
     &              parameters,1,tgt_info)

      ! define transposed density matrix (e.g. for S^(-1/2))
      call add_target('Dtr',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
            if (.not.pure_vv.and.ip.eq.0.and.ih.eq.0.and.ivv.gt.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if (max(ip,ih)+ivv.lt.minexc) cycle
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*2-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,2,ndef*2-1) = max(ih-ip,0) + ivv
            occ_def(IVALE,1,ndef*2) = max(ip-ih,0) + ivv
            occ_def(IVALE,2,ndef*2) = max(ip-ih,0) + ivv
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('Dtr',ttype_op,DEF_OP_FROM_OCC,
     &              'Dtr',1,1,
     &              parameters,2,tgt_info)

      ! define Jacobian
      call add_target('A_C',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (.not.pure_vv.and.ivv.gt.0.and.ip.eq.0.and.ih.eq.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if ((max(ip,ih).gt.0.or.pure_vv)
     &          .and.max(ip,ih)+ivv.lt.minexc) cycle
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
      call set_rule('A_C',ttype_op,DEF_OP_FROM_OCC,
     &              'A_C',1,1,
     &              parameters,2,tgt_info)

      ! define Hessian
      call add_target('A',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do ivv = 0, min(max(max(maxp,maxh),maxexc)-max(ip,ih),maxvv)
            if (.not.pure_vv.and.ivv.gt.0.and.ip.eq.0.and.ih.eq.0) cycle
            if (abs(ih-ip)+2*ivv.gt.maxv) cycle
            if (max(ip,ih).gt.0.and.ip.lt.minp) cycle
            if (max(ip,ih).gt.0.and.ih.lt.minh) cycle
            if ((max(ip,ih).gt.0.or.pure_vv)
     &          .and.max(ip,ih)+ivv.lt.minexc) cycle
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
      use_hessian = ndef.gt.0
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('A',ttype_op,DEF_OP_FROM_OCC,
     &              'A',1,1,
     &              parameters,2,tgt_info)

      ! Metric times icCI vector product
      call add_target('SC',ttype_op,.false.,
     &                tgt_info)
      call set_dependency('SC','A_C',tgt_info)
      call cloneop_parameters(-1,parameters,'A_C',.false.)
      call set_rule('SC',ttype_op,CLONE_OP,'SC',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner
      call add_target(op_dia//'_C',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(op_dia//'_C','C',tgt_info)
      call cloneop_parameters(-1,parameters,'C',.false.)
      call set_rule(op_dia//'_C',ttype_op,CLONE_OP,op_dia//'_C',1,1,
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
     &       val_label=(/'C0^+','Cdag',op_ham,'C','C0'/))
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
     &       val_label=(/'Cdag',op_ham,'C'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/2,3,4/))
        call set_rule2('F_preE(MR)',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_preE(MR)',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'DENS','Cdag',op_ham,'C','DENS'/))
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
     &       val_label=(/op_ham,'C0','DENS'/)) ! C0 is dummy
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
     &       val_label=(/op_ham,'C0','CUM'/)) ! C0 is dummy
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
     &       val_label=(/'Cdag',op_ham/))
        call set_dependency('F_preE(MR)','1',tgt_info)
        call set_rule2('F_preE(MR)',REPLACE,tgt_info)
        call set_arg('F_preE(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1v','1'/))
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
     &       val_label=(/'Cdag','C'/))
        call set_rule2('F_preE(MR)',REPLACE,tgt_info)
        call set_arg('F_preE(MR)',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MR)'/))
        call set_arg('F_preE(MR)',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1v','1'/))
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
     &       val_label=(/'1v','1'/))
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
     &     val_label=(/'Cdag','C^+'/))
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
      call set_dependency('F_E(MR)_diag','1ph',tgt_info)
      call set_rule2('F_E(MR)_diag',INSERT,tgt_info)
      call set_arg('F_E(MR)_diag',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MR)_diag'/))
      call set_arg('F_E(MR)_diag',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MR)'/))
      call set_arg('F_E(MR)_diag',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_E(MR)_diag',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1ph'/))
      call set_arg('F_E(MR)_diag',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/'Cdag','C'/))
      ! replace 1ph by 1
      call set_dependency('F_E(MR)_diag','1',tgt_info)
      call set_rule2('F_E(MR)_diag',REPLACE,tgt_info)
      call set_arg('F_E(MR)_diag',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MR)_diag'/))
      call set_arg('F_E(MR)_diag',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MR)_diag'/))
      call set_arg('F_E(MR)_diag',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1ph','1'/))
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
      labels(3) = 'A_C'
      labels(4) = 'Ctr^+'
      labels(5) = ' '
      call add_target('F_A_C_diag',ttype_frm,.false.,tgt_info)
      call set_dependency('F_A_C_diag','F_E(MR)_diag',tgt_info)
      call set_dependency('F_A_C_diag','A_C',tgt_info)
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

      ! expression for the norm
      ! a) set up norm expression
      op_deexc = 'C^+ '
      op_exc = 'C   '
      if (gno.eq.1) then
        op_deexc = 'c^+ '
        op_exc = 'c    '
      end if
      call add_target2('F_NORM',.false.,tgt_info)
      call set_dependency('F_NORM','NORM',tgt_info)
      call set_dependency('F_NORM',op_exc,tgt_info)
c      if (gno.eq.0) then
c        call set_dependency('F_NORM','C0',tgt_info)
c        call set_rule2('F_NORM',EXPAND_OP_PRODUCT,tgt_info)
c        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &       val_label=(/'F_NORM'/))
c        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &       val_label=(/'NORM'/))
c        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &       tgt_info,
c     &       val_label=(/'C0^+','C^+','C','C0'/))
c        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &       val_int=(/2,3,4,5/))
c      else if (gno.eq.1) then
        call set_dependency('F_NORM','DENS',tgt_info)
        call set_rule2('F_NORM',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &       tgt_info,
     &       val_label=(/op_deexc,op_exc/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &       val_int=(/2,3/))
        call set_rule2('F_NORM',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &       tgt_info,
     &       val_label=(/'DENS',op_deexc,op_exc,'DENS'/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &       val_int=(/2,3,4,2/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/1/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
     &       val_int=(/1,4/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
        if (gno.eq.0)
     &   call set_arg('F_NORM',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &       val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
c      end if
      ! b) insert unit operators to allow for differentiation
      ! and for factoring out of hole densities
      call set_dependency('F_NORM','1v',tgt_info)
      call set_rule2('F_NORM',INSERT,tgt_info)
      call set_arg('F_NORM',INSERT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',INSERT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',INSERT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM',INSERT,'OP_INS',1,tgt_info,
     &     val_label=(/'1v'/))
      call set_arg('F_NORM',INSERT,'OP_INCL',2,tgt_info,
     &     val_label=(/op_deexc,op_exc/))
      ! c) replace 1v by 1 (was used because we only needed valence blocks)
      call set_dependency('F_NORM','1',tgt_info)
      call set_rule2('F_NORM',REPLACE,tgt_info)
      call set_arg('F_NORM',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1v','1'/))
      if (gno.eq.1) then
        ! d) expand reduced densities in terms of cumulants
        call set_dependency('F_NORM','F_DENS',tgt_info)
        call set_rule2('F_NORM',EXPAND,tgt_info)
        call set_arg('F_NORM',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_DENS'/))
c        call set_rule2('F_NORM',PRINT_FORMULA,tgt_info)
c        call set_arg('F_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
c       &     val_label=(/'F_NORM'/))
        ! e) select only terms allowed according to contraction rules
        call set_rule2('F_NORM',SELECT_SPECIAL,tgt_info)
        call set_arg('F_NORM',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',SELECT_SPECIAL,'OPERATORS',3,tgt_info,
     &       val_label=(/op_ham,'C0','CUM'/)) !op_ham and C0 are dummies
        call set_arg('F_NORM',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC')
        ! f) factor out hole density
        call set_dependency('F_NORM','F_HOLE',tgt_info)
        call set_rule2('F_NORM',FACTOR_OUT,tgt_info)
        call set_arg('F_NORM',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_HOLE'/))
      end if
      call set_rule2('F_NORM',PRINT_FORMULA,tgt_info)
      call set_arg('F_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_NORM'/))
c      labels(1:20)(1:len_target_name) = ' '
c      labels(1) = 'F_NORM'
c      labels(2) = 'NORM'
c      labels(3) = 'C0^+'
c      labels(4) = 'C^+'
c      labels(5) = '1' ! pure valence blocks only
c      labels(6) = 'C'
c      labels(7) = 'C0'
c      call add_target('F_NORM',ttype_frm,.false.,tgt_info)
c      call set_dependency('F_NORM','NORM',tgt_info)
c      call set_dependency('F_NORM','C0',tgt_info)
c      call set_dependency('F_NORM','C',tgt_info)
c      call set_dependency('F_NORM','1',tgt_info)
c      call expand_parameters(-1,
c     &     parameters,3,
c     &     'multireference energy expression',5,
c     &     (/2,3,4,5,6/),
c     &     (/-1,-1,-1,-1,-1/),
c     &     (/0,0,maxexc+1,0,0/),
c     &     0,0,
c     &     (/1,3,3,5/),2,
c     &     0,0)
c      call set_rule('F_NORM',ttype_frm,EXPAND_OP_PRODUCT,
c     &              labels,7,1,
c     &              parameters,3,tgt_info)
c      ! delete terms in which C^+ and C are contracted via valence lines
c      labels(2:10)(1:len_target_name) = ' '
c      labels(2) = 'F_NORM'
c      labels(3) = 'NORM'
c      labels(4) = 'C^+'
c      labels(5) = 'C'
c      call form_parameters(-1,
c     &       parameters,2,'norm expression',3,'delete')
c      call set_rule('F_NORM',ttype_frm,SELECT_LINE,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)
c      labels(3:10)(1:len_target_name) = ' '
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_NORM',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! factor out reduced density matrices
c      if (use_dens) then
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
     &     val_label=(/'D','C^+','D','C','D'/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &     val_int=(/2,3,2,4,2/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/5/))
      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'AVOID',10,tgt_info,
     &     val_int=(/1,3,1,4,1,5,2,5,3,5/))
c      call set_rule2('F_NORM_fact',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORM_fact'/))
c      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'OPERATORS',2,
c     &     tgt_info,
c     &     val_label=(/'C^+','C'/))
c      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_arg('F_NORM_fact',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
      ! no active lines between C^+ and C
      call set_rule2('F_NORM_fact',SELECT_LINE,tgt_info)
      call set_arg('F_NORM_fact',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_NORM_fact'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM_fact',SELECT_LINE,'OP_INCL',2,tgt_info,
     &     val_label=(/'C^+','C'/))
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
      if (use_dens) then
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

      ! reduced density matrices
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_SC0'
      labels(2) = 'F_NORM'
      labels(3) = 'SC'
      labels(4) = 'C^+'
      if (gno.eq.1) labels(4) = 'c^+'
      labels(5) = ' '
      call add_target('F_SC0',ttype_frm,.false.,tgt_info)
      call set_dependency('F_SC0','F_NORM',tgt_info)
      call set_dependency('F_SC0','SC',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'precursor for density matrices',1,'---')
      call set_rule('F_SC0',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_SC0',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! reduced density matrices
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_D'
      labels(2) = 'F_SC0'
      labels(3) = 'D'
      labels(4) = 'C'
      if (gno.eq.1) labels(4) = 'c'
      labels(5) = ' '
      call add_target('F_D',ttype_frm,.false.,tgt_info)
      call set_dependency('F_D','F_SC0',tgt_info)
      call set_dependency('F_D','D',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'Reduced density matrices',1,'---')
      call set_rule('F_D',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
      ! delete terms which do not have any open lines (valence space)
      ! => scalar occupation class is excluded in the definition
      labels(2:20)(1:len_target_name) = ' '
      labels(2) = 'F_D'
      labels(3) = 'D'
      labels(4) = '1'
c      labels(5) = 'C0'
      labels(5) = 'DENS'
      if (gno.eq.1) then
        labels(4) = 'HOLE'
        labels(5) = 'CUM'
      end if
      call form_parameters(-1,
     &       parameters,2,'D formula',3,'ext')
      call set_rule('F_D',ttype_frm,SELECT_LINE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c dbg
c      call set_rule2('F_D',KEEP_TERMS,tgt_info)
c      call set_arg('F_D',KEEP_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_D'/))
c      call set_arg('F_D',KEEP_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_D'/))
c      call set_arg('F_D',KEEP_TERMS,'TERMS',2,tgt_info,
c     &     val_int=(/6,9/))
c dbgend
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_D',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)

      ! daggered density
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_Ddag'
      call add_target('F_Ddag',ttype_frm,.false.,tgt_info)
      call set_dependency('F_Ddag','Ddag',tgt_info)
      call set_dependency('F_Ddag','D',tgt_info)
      call def_form_parameters(-1,
     &     parameters,2,'Ddag=D','adjoint of density matrix')
      call set_rule('F_Ddag',ttype_frm,DEF_FORMULA,
     &              labels,1,1,
     &              parameters,2,tgt_info)
      labels(2) = 'F_Ddag'
      labels(3) = 'D'
      labels(4) = 'D^+'
      call form_parameters(-1,
     &     parameters,2,'adjoint of density matrix',1,
     &     '---')
      call set_rule('F_Ddag',ttype_frm,REPLACE,
     &            labels,4,1,
     &            parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_Ddag',ttype_frm,PRINT_FORMULA,
c     &                labels,1,0,parameters,2,tgt_info)

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
      if (use_dens) then
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
      if (use_dens.or..true.)
     &        call set_dependency('FOPT_NORM','DEF_ME_D',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_NORM',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! density matrix
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_D'
      labels(2) = 'F_D'
      call add_target('FOPT_D',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_D','F_D',tgt_info)
c      call set_dependency('FOPT_D','DEF_ME_C0',tgt_info)
      if (gno.eq.0) then
        call set_dependency('FOPT_D','DEF_ME_DENS',tgt_info)
        call set_dependency('FOPT_D','DEF_ME_1',tgt_info)
      else if (gno.eq.1) then
        call set_dependency('FOPT_D','DEF_ME_CUM',tgt_info)
        call set_dependency('FOPT_D','DEF_ME_HOLE',tgt_info)
      end if
      call set_dependency('FOPT_D','DEF_ME_D',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_D',ttype_frm,OPTIMIZE,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! daggered density matrix
      labels(1:20)(1:len_target_name)= ' '
      labels(1) = 'FOPT_Ddag'
      labels(2) = 'F_Ddag'
      call add_target('FOPT_Ddag',ttype_frm,.false.,tgt_info)
      call set_dependency('FOPT_Ddag','F_Ddag',tgt_info)
      call set_dependency('FOPT_Ddag','DEF_ME_D',tgt_info)
      call set_dependency('FOPT_Ddag','DEF_ME_Dinvdag',tgt_info)
      call opt_parameters(-1,parameters,1,0)
      call set_rule('FOPT_Ddag',ttype_frm,OPTIMIZE,
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
      call set_arg('DEF_ME_C',DEF_ME_LIST,'MS',1,tgt_info,
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
      call set_arg('DEF_ME_Ctr',DEF_ME_LIST,'MS',1,tgt_info,
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

      ! ME_A
      call add_target2('DEF_ME_A',.false.,tgt_info)
      call set_dependency('DEF_ME_A','A',tgt_info)
      call set_rule2('DEF_ME_A',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_A',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_A'/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'A'/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &     val_int=(/0/))

      ! ME_1
      call add_target2('DEF_ME_1',.false.,tgt_info)
      call set_dependency('DEF_ME_1','1',tgt_info)
      call set_rule2('DEF_ME_1',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_1',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_1'/))
      call set_arg('DEF_ME_1',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'1'/))
      call set_arg('DEF_ME_1',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_1',DEF_ME_LIST,'MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_1',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call dens_parameters(-1,parameters,0,0,0)
      call set_rule('DEF_ME_1',ttype_opme,UNITY,
     &     'ME_1',1,1,
     &     parameters,1,tgt_info)

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

      ! ME_E(MR)
      call add_target('DEF_ME_E(MR)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_E(MR)','E(MR)',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_E(MR)'
      labels(2) = 'E(MR)'
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DEF_ME_E(MR)',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_NORM
      call add_target('DEF_ME_NORM',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_NORM','NORM',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_NORM'
      labels(2) = 'NORM'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_NORM',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! ME_D
      call add_target('DEF_ME_D',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_D','D',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_D'
      labels(2) = 'D'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_D',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! inverted ME_D
      call add_target('DEF_ME_Dinv',ttype_opme,.false.,tgt_info)
      if (calc) then
        call set_dependency('DEF_ME_Dinv','EVAL_D',tgt_info)
      else
        call set_dependency('DEF_ME_Dinv','EVAL_MRCC_D',tgt_info)
      end if
      call set_dependency('DEF_ME_Dinv','DEF_ME_D',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dinv'
      labels(2) = 'D'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dinv',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      labels(1) = 'ME_D'
      labels(2) = 'D'
      call set_rule('DEF_ME_Dinv',ttype_opme,ASSIGN_ME2OP,
     &     labels,2,1,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     '---',0,'invsqrthalf')
      labels(1) = 'ME_Dinv'   ! output
      labels(2) = 'ME_D'      ! input
      call set_rule('DEF_ME_Dinv',ttype_opme,INVERT,
     &              labels,2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Square root of inverse density matrix :',0,'LIST')
c      labels(1) = 'DEF_ME_Dinv'
c      labels(2) = 'ME_Dinv'
c      call set_rule('DEF_ME_Dinv',ttype_opme,PRINT_MEL,
c     &     'ME_Dinv',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! inverted ME_Ddag
      call add_target('DEF_ME_Dinvdag',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dinvdag','Ddag',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dinvdag'
      labels(2) = 'Ddag'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dinvdag',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! reordered inverted ME_D
      call add_target('DEF_ME_Dtr',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dtr','Dtr',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dtr'
      labels(2) = 'Dtr'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dtr',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_dependency('DEF_ME_Dtr','DEF_ME_Dinv',tgt_info)
      labels(2) = 'ME_Dinv'
      call form_parameters(-1,parameters,2,
     &     '---',13,'---')
      call set_rule('DEF_ME_Dtr',ttype_opme,
     &              REORDER_MEL,
     &              labels,2,1,
     &              parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Reordered inverted Density matrix :',0,'LIST')
c      labels(1) = 'DEF_ME_Dtr'
c      labels(2) = 'ME_Dtr'
c      call set_rule('DEF_ME_Dtr',ttype_opme,PRINT_MEL,
c     &     'ME_Dtr',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! reordered daggered inverted ME_D
      call add_target('DEF_ME_Dtrdag',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dtrdag','Dtr',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dtrdag'
      labels(2) = 'Dtr'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dtrdag',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_dependency('DEF_ME_Dtrdag','EVAL_Dinvdag',tgt_info)
      labels(2) = 'ME_Dinvdag'
      call form_parameters(-1,parameters,2,
     &     '---',13,'---')
      call set_rule('DEF_ME_Dtrdag',ttype_opme,
     &              REORDER_MEL,
     &              labels,2,1,
     &              parameters,2,tgt_info)

      ! Diagonal Preconditioner
      call me_list_label(dia_label,mel_dia,1,
     &     0,0,0,.false.)
      call add_target(trim(dia_label)//'C',ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label)//'C','EVAL_FREF',tgt_info)
      call set_dependency(trim(dia_label)//'C',
     &                    op_dia//'_'//'C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)//'C'
      labels(2) = op_dia//'_'//'C'
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
      call set_rule(trim(dia_label)//'C',ttype_opme,
     &              PRECONDITIONER,
     &              labels,2,1,
     &              parameters,3,tgt_info)

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Evaluate norm
      call add_target('EVAL_NORM',ttype_gen,calc,tgt_info)
      call set_dependency('EVAL_NORM','FOPT_NORM',tgt_info)
      call set_dependency('EVAL_NORM','SOLVE_ICCI',tgt_info)
      call set_rule('EVAL_NORM',ttype_opme,EVAL,
     &     'FOPT_NORM',1,0,
     &     parameters,0,tgt_info)

      ! Evaluate density matrix
      call add_target('EVAL_D',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_D','FOPT_D',tgt_info)
      call set_dependency('EVAL_D','EVAL_REF_S(S+1)',tgt_info)
      if (gno.eq.0) then
        call set_dependency('EVAL_D','EVAL_DENS0',tgt_info)
      else if (gno.eq.1) then
        call set_dependency('EVAL_D','EVAL_HOLE',tgt_info)
        call set_dependency('EVAL_D','EVAL_CUM',tgt_info)
      end if
      call set_rule('EVAL_D',ttype_opme,EVAL,
     &     'FOPT_D',1,0,
     &     parameters,0,tgt_info)
      ! fix: set first element (zero occ) to 1 (had not been defined in F_D)
      labels(1:10)(1:len_target_name) = ' '
      labels(1) = 'ME_D'
      call dens_parameters(-1,parameters,1,1,1)
      call set_rule('EVAL_D',ttype_opme,UNITY,
     &     labels,1,1,
     &     parameters,1,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'Density matrix :',0,'LIST')
c      call set_rule('EVAL_D',ttype_opme,PRINT_MEL,
c     &     'ME_D',1,0,
c     &     parameters,2,tgt_info)

      ! Calculate adjoint of ME_Dinv
      call add_target('EVAL_Dinvdag',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_Dinvdag','FOPT_Ddag',tgt_info)
      call set_dependency('EVAL_Dinvdag','DEF_ME_Dinv',tgt_info)
      labels(1) = 'ME_Dinv'
      labels(2) = 'D'
      call set_rule('EVAL_Dinvdag',ttype_opme,ASSIGN_ME2OP,
     &     labels,2,1,
     &     parameters,0,tgt_info)
      call set_rule('EVAL_Dinvdag',ttype_opme,EVAL,
     &     'FOPT_Ddag',1,0,
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
      labels(1) = trim(dia_label)//'C'
      labels(2) = 'ME_A'
      call me_list_label(dia_label,mel_dia,1,
     &     0,0,0,.false.)
      call set_dependency('EVAL_A_diag',trim(dia_label)//'C',tgt_info)
      call set_rule('EVAL_A_diag',ttype_opme,
     &              EXTRACT_DIAG,
     &              labels,2,1,
     &              parameters,0,tgt_info)
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
      if (use_hessian)
     &       call set_dependency('SOLVE_ICCI','EVAL_A_diag',tgt_info)
      call set_dependency('SOLVE_ICCI',trim(dia_label)//'C',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      if (use_dens) then
        labels(1) = 'ME_D'
        labels(2) = 'D'
        call set_rule('SOLVE_ICCI',ttype_opme,ASSIGN_ME2OP,
     &       labels,2,1,
     &       parameters,0,tgt_info)
      end if
      call solve_parameters(-1,parameters,2,1,nroots,'DIA')
      labels(1) = 'ME_C'
      labels(2) = trim(dia_label)//'C'
      labels(3) = 'A_C'
      labels(4) = 'SC'
      labels(5) = 'FOPT_A_C'
      nlabels = 5
      if (use_dens) then
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
      call add_target('EVAL_E(MR)',ttype_gen,calc,tgt_info)
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
