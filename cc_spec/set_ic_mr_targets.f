*----------------------------------------------------------------------*
      subroutine set_ic_mr_targets(tgt_info,orb_info,
     &                             excrestr,maxh,maxp,use_met,
     &                             nsupD,stndD,nremblk,remblk)
*----------------------------------------------------------------------*
*     sets targets that are common for internally contracted
*     multireference methods
*
*     matthias, march 2011
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
     &     ntest =  00

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     maxh, maxp, excrestr(0:maxh,0:maxp,1:2)
      logical, intent(out) ::
     &     use_met

      integer, intent(out) ::
     &     nsupD, stndD(2,60), nremblk, remblk(60)
      logical ::
     &     first

      integer ::
     &     ndef, occ_def(ngastp,2,124),!60),
     &     msc, maxexc, ip, ih, iv,
     &     gno, idef, iexc, jexc,
     &     version(60), ivers, icnt, prc_type, project,
     &     ciroot, n_states, i_state, optref
      logical ::
     &     l_exist,
     &     l_icci, l_iccc, skip, Op_eqs, svdonly, prc_traf,
     &     multistate
      real(8) ::
     &     prc_shift, densmix
      character(len_target_name) ::
     &     me_label, medef_label, dia_label, mel_dia1,
     &     labels(20), c_st
      character(len_command_par) ::
     &     parameters(3)
      character*4 ::
     &     op_exc, op_deexc

      character(len_target_name), external ::
     &     state_label

      if (iprlvl.gt.0) write(lulog,*)
     &     'setting targets for internally contracted MR methods'

      ! get maximum excitation rank
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)

      ! icMRCI calculation?
      l_icci = is_keyword_set('method.MRCI').gt.0
      ! icMRCC calculation?
      l_iccc = is_keyword_set('method.MRCC').gt.0

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1
      if (orb_info%ims.ne.0) msc = 0

      ! which normal ordering is used?
      call get_argument_value('method.MR','GNO',
     &     ival=gno)
      select case(gno)
      case(0)
      case(1)
        write(lulog,*) 'Using generalized normal order (GNO)'
c        call quit(1,'set_ic_mr_targets','Use of GNO not debugged yet')
      case default
        call quit(1,'set_ic_mr_targets','unknown normal order')
      end select

      call get_argument_value('method.MR','project',
     &     ival=project)
      call get_argument_value('method.MR','svdonly',
     &     lval=svdonly)
      call get_argument_value('method.MR','prc_type',
     &     ival=prc_type)
      call get_argument_value('method.MR','prc_shift',
     &     xval=prc_shift)
      call get_argument_value('method.MR','prc_traf',
     &     lval=prc_traf)
      call get_argument_value('method.MR','densmix',
     &     xval=densmix)
      call get_argument_value('method.MR','multistate',
     &     lval=multistate)
      call get_argument_value('method.MR','ciroot',
     &     ival=ciroot)
      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)

      if(multistate)then
       n_states = ciroot
      else
       n_states = 1
      end if

      if (.not.l_iccc.and.prc_type.ne.0.and.prc_type.ne.3.or.
     &    prc_type.gt.4.or.prc_type.ne.2.and.prc_shift.ne.0d0)
     &  call quit(1,'set_ic_mr_targets','Choose other preconditioner!')
      if (.not.l_iccc.and.prc_traf)
     &   write(lulog,'(a)') 'Ignoring prc_traf option in ic-MRCI'
c     &  call quit(1,'set_ic_mr_targets','prc_traf only for MRCC yet')

      if (ntest.ge.100) then
        print *,'gno     = ',gno
        print *,'sv_thr. = ',sv_thresh
        print *,'Tikhonov= ',tikhonov
        print *,'sv_fix  = ',sv_fix
        print *,'sv_old  = ',sv_old
        print *,'project = ',project
        print *,'prc_type= ',prc_type
        print *,'prc_shift=',prc_shift
      end if

      if (sv_fix) then
        inquire(file='SINGVALS',exist=l_exist)
        if (l_exist.and.sv_old) write(lulog,*)
     &     'Using existing SINGVALS file for singular value selection!'
      end if
      if (l_iccc.and.prc_traf.and.jac_fix) then
        inquire(file='SINGVALS2',exist=l_exist)
        if (l_exist) write(lulog,*)
     &     'Using existing SINGVALS2 file for singular value selection!'
      end if

      if (l_iccc) then
        call get_argument_value('method.MRCC','Op_eqs',
     &       lval=Op_eqs)
        if (project.eq.0.and.Op_eqs) 
     &     project=1 ! no off-diagonal metric blocks
      end if

*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define scalar multireference energy
      call add_target('E(MR)',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('E(MR)',ttype_op,DEF_HAMILTONIAN,'E(MR)',
     &              1,1,parameters,1,tgt_info)
      if (multistate)
     &     call set_multistate_operator(tgt_info,n_states,'E(MR)',1)

      ! define unit operator suited to connect C^+ and C
      call add_target('1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
           do iv = 0, maxexc
            ! no scalar block needed
            if (ip+ih+iv.eq.0) cycle
c            ! not needed for pure inactive excitations
c            if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
            ! only pure active or pure non-active occ. classes needed for now
            if (iv.gt.0.and.ih+ip.gt.0) cycle
            ! actually we only need the one-particle blocks!
            if (ih+ip.gt.1) cycle
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
          !no scalar block needed
          if (ip+ih.eq.0) cycle
c          ! not needed for pure inactive excitations
c          if (ip.eq.maxexc.and.ih.eq.maxexc) cycle
c          if (ip.ge.2.and.ih.ge.2) cycle
          ! actually we only need the one-particle blocks!
          if (ip+ih.gt.1) cycle
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
      do iv = 1, maxexc !max(maxexc,orb_info%nactel)
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

      ! define scalar unit operator
      call add_target('1scal',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('1scal',ttype_op,DEF_OP_FROM_OCC,
     &              '1scal',1,1,
     &              parameters,2,tgt_info)
      call opt_parameters(-1,parameters,+1,0)
      call set_rule('1scal',ttype_op,SET_HERMIT,
     &              '1scal',1,1,
     &              parameters,1,tgt_info)

      ! define scalar norm
      call add_target('NORM',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('NORM',ttype_op,DEF_HAMILTONIAN,'NORM',
     &              1,1,parameters,1,tgt_info)

      ! define density matrix
      call add_target('D',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      nsupD = 0
      nremblk = 0
      do ip = 0, maxp
        do ih = 0, maxh
           first = .true.
          do iexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
           do jexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
            if (project.eq.1.and.gno.eq.0.and.iexc.ne.jexc) cycle
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))) cycle
            if(.not.l_iccc.and.
     &           ip.eq.0.and.ih.eq.0.and.iexc.eq.0.and.jexc.eq.0)
     &           cycle
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*3-1) = iexc - ip
            occ_def(IVALE,2,ndef*3-1) = jexc - ip
            occ_def(IVALE,1,ndef*3) = jexc - ih
            occ_def(IVALE,2,ndef*3-2) = iexc - ih
            if (first) then
               nsupD = nsupD + 1
               stndD(1,nsupD) = ndef
               first = .false.
            end if
            stndD(2,nsupD) = ndef
            if (ip+ih.eq.1) then
              nremblk = nremblk + 1
              remblk(nremblk) = 2*(ndef-1)+1
            end if
           end do
          end do
        end do
      end do
      use_met = ndef.gt.0 !1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('D',ttype_op,DEF_OP_FROM_OCC,
     &              'D',1,1,
     &              parameters,2,tgt_info)

      ! define transposed density matrix (e.g. for S^(-1/2))
      call add_target('Dtr',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
           do jexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
            if (project.eq.1.and.gno.eq.0.and.iexc.ne.jexc) cycle
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))) cycle
            if(.not.l_iccc.and.
     &           ip.eq.0.and.ih.eq.0.and.iexc.eq.0.and.jexc.eq.0)
     &           cycle
            ndef = ndef + 1
            occ_def(IVALE,1,ndef*2-1) = iexc - ip
            occ_def(IVALE,2,ndef*2-1) = jexc - ip
            occ_def(IVALE,1,ndef*2) = jexc - ih
            occ_def(IVALE,2,ndef*2) = iexc - ih
           end do
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('Dtr',ttype_op,DEF_OP_FROM_OCC,
     &              'Dtr',1,1,
     &              parameters,2,tgt_info)

      ! we explicitly define an operator for the projector onto
      ! non-redundant excitation space 
      call add_target3(
     &     (/'target Dproj(CLONE_OPERATOR(label=Dproj,template=Dtr))'/),
     &     tgt_info)

      ! subset of exc. op. for non-redundant valence-only metric
      call add_target('Ex',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          ! more excitations than in related class is problematic
          ! in the derivation of off-diagonal metric blocks later on
          if (ip.ge.1.and.ih.ge.1.and.
     &        .not.(project.eq.1.and.gno.eq.0)) then
            iexc = excrestr(ih-1,ip-1,2)-excrestr(ih-1,ip-1,1)
            if (iexc.ge.0.and.
     &          excrestr(ih,ip,2)-excrestr(ih,ip,1).gt.iexc) then
              write(lulog,'(3(a,i1),a)') 'In the excitation class n_h=',
     &          ih,',n_p=',ip,','
              write(lulog,'(3(a,i1),a)')
     &          'more excitations are defined than in the class n_h=',
     &          ih-1,',n_p=',ip-1,'.'
              write(lulog,'(6(a,i1),a)') 'Try to use MR excrestr=(',
     &          ih,',',ih,',',ip,',',ip,',',excrestr(ih,ip,1),
     &          ',',excrestr(ih,ip,1)+iexc,')'
              call quit(1,'set_ic_mr_targets',
     &                  'Check your excitation class restrictions!')
            end if
          end if
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))
     &          .and..not.(l_icci.and.prc_type.ge.3)) cycle
            ! same valence structure already exists?
            ivers = 1
            do idef = 1, ndef
              if (occ_def(IVALE,2,idef).eq.iexc-ih
     &            .and.occ_def(IVALE,1,idef).eq.iexc-ip)
     &           ivers = ivers + 1
            end do
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
            ! distinguish ops with same valence part by blk_version
            version(ndef) = ivers
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('Ex',ttype_op,DEF_OP_FROM_OCC,
     &              'Ex',1,1,
     &              parameters,2,tgt_info)
      call set_rule2('Ex',SET_ORDER,tgt_info)
      call set_arg('Ex',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'Ex'/))
      call set_arg('Ex',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/1/))
      call set_rule2('Ex',SET_ORDER,tgt_info)
      call set_arg('Ex',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'Ex'/))
      call set_arg('Ex',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('Ex',SET_ORDER,'ORDER',1,tgt_info,
     &             val_int=(/ndef/))
      call set_arg('Ex',SET_ORDER,'IDX_FREQ',ndef,tgt_info,
     &             val_int=version(1:ndef))

      ! subset of Residual
      call add_target('OMGred',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))
     &          .and..not.(l_icci.and.prc_type.ge.3)) cycle
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
      call set_rule('OMGred',ttype_op,DEF_OP_FROM_OCC,
     &              'OMGred',1,1,
     &              parameters,2,tgt_info)
      ! distinguish ops with same valence part by blk_version
      call set_rule2('OMGred',SET_ORDER,tgt_info)
      call set_arg('OMGred',SET_ORDER,'LABEL',1,tgt_info,
     &     val_label=(/'OMGred'/))
      call set_arg('OMGred',SET_ORDER,'SPECIES',1,tgt_info,
     &             val_int=(/-1/))
      call set_arg('OMGred',SET_ORDER,'ORDER',1,tgt_info,
     &             val_int=(/ndef/))
      call set_arg('OMGred',SET_ORDER,'IDX_FREQ',ndef,tgt_info,
     e             val_int=version(1:ndef))

      ! define transformed Residual (for preconditioner)
      call add_target('OMGtr',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      icnt = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))
     &          .and..not.(l_icci.and.prc_type.ge.3)) cycle
            icnt = icnt + 1
            ! for cheap precond.: only valence part needed eventually
            if (prc_type.ge.3.and.version(icnt).ne.1) cycle
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
      call set_rule('OMGtr',ttype_op,DEF_OP_FROM_OCC,
     &              'OMGtr',1,1,
     &              parameters,2,tgt_info)

      ! define Hessian / Jacobian (diagonal blocks only)
      call add_target('A',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      idef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))
     &          .and..not.(l_icci.and.prc_type.ge.3)) cycle
c    strictly, we should better use this: (but: does procedure C work?)
c            if (ip.eq.ih.and.ip.eq.iexc) cycle
c dbg hybrid preconditioner
c           if (ndef.ge.6) cycle
c dbgend
            idef = idef + 1
            ! for cheap precond.: only valence part needed
            if (prc_type.ge.3.and.
     &          (version(idef).ne.1.or.ip.eq.ih.and.ip.eq.iexc
     &                                 .and..not.l_icci)) cycle
            ndef = ndef + 1
            if (prc_type.lt.3) then
              occ_def(IHOLE,1,3*ndef-1) = ih
              occ_def(IHOLE,2,3*ndef-1) = ih
              occ_def(IPART,1,3*ndef-1) = ip
              occ_def(IPART,2,3*ndef-1) = ip
            end if
            occ_def(IVALE,1,3*ndef-1) = iexc - ip
            occ_def(IVALE,2,3*ndef-1) = iexc - ip
            occ_def(IVALE,2,3*ndef-2) = iexc - ih
            occ_def(IVALE,1,3*ndef) = iexc - ih
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
      call set_rule('A',ttype_op,DEF_OP_FROM_OCC,
     &              'A',1,1,
     &              parameters,2,tgt_info)

*----------------------------------------------------------------------*
*     Formulae 
*----------------------------------------------------------------------*

      ! expression for the norm
      ! a) set up norm expression
        op_deexc = 'Ex^+'
        op_exc = 'Ex  '
      call add_target2('F_NORM',.false.,tgt_info)
      call set_dependency('F_NORM','NORM',tgt_info)
      call set_dependency('F_NORM',op_exc,tgt_info)
      call set_rule2('F_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_NORM',EXPAND_OP_PRODUCT,'OPERATORS',2,
     &     tgt_info,
     &     val_label=(/op_deexc,op_exc/))
      call set_arg('F_NORM',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
     &     val_int=(/2,3/))
      if (orb_info%nactel.gt.0) then
        call set_dependency('F_NORM','DENS',tgt_info)
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
c        if (gno.eq.0)
c       & call set_arg('F_NORM',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c       &     val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
        call set_arg('F_NORM',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &       val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
      end if
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
      if (l_iccc) then
        call set_dependency('F_NORM','1scal',tgt_info)
        call set_rule2('F_NORM',INSERT,tgt_info)
        call set_arg('F_NORM',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_NORM'/))
        call set_arg('F_NORM',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_NORM',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1scal'/))
        call set_arg('F_NORM',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/op_deexc,op_exc/))
      end if
      ! c) replace 1v by 1 (was used because we only needed valence blocks)
      call set_dependency('F_NORM','1',tgt_info)
      call set_rule2('F_NORM',REPLACE,tgt_info)
      call set_arg('F_NORM',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_NORM'/))
      call set_arg('F_NORM',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'1v','1 '/))
      ! the following is commented out since our new strategy is to
      ! compute the metric in the standard NO and then transform to GNO
c      if (gno.eq.1) then
c        ! d) expand reduced densities in terms of cumulants
c        call set_dependency('F_NORM','F_DENS',tgt_info)
c        call set_rule2('F_NORM',EXPAND,tgt_info)
c        call set_arg('F_NORM',EXPAND,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_NORM'/))
c        call set_arg('F_NORM',EXPAND,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_NORM'/))
c        call set_arg('F_NORM',EXPAND,'INTERM',1,tgt_info,
c     &       val_label=(/'F_DENS'/))
cc        call set_rule2('F_NORM',PRINT_FORMULA,tgt_info)
cc        call set_arg('F_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
cc       &     val_label=(/'F_NORM'/))
c        ! e) select only terms allowed according to contraction rules
c        call set_rule2('F_NORM',SELECT_SPECIAL,tgt_info)
c        call set_arg('F_NORM',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_NORM'/))
c        call set_arg('F_NORM',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_NORM'/))
c        call set_arg('F_NORM',SELECT_SPECIAL,'OPERATORS',3,tgt_info,
c     &       val_label=(/op_ham,'C0','CUM'/)) !op_ham and C0 are dummies
c        call set_arg('F_NORM',SELECT_SPECIAL,'TYPE',1,tgt_info,
c     &       val_str='MRCC')
cc        ! f) factor out hole density
cc        call set_dependency('F_NORM','F_HOLE',tgt_info)
cc        call set_rule2('F_NORM',FACTOR_OUT,tgt_info)
cc        call set_arg('F_NORM',FACTOR_OUT,'LABEL_RES',1,tgt_info,
cc     &       val_label=(/'F_NORM'/))
cc        call set_arg('F_NORM',FACTOR_OUT,'LABEL_IN',1,tgt_info,
cc     &       val_label=(/'F_NORM'/))
cc        call set_arg('F_NORM',FACTOR_OUT,'INTERM',1,tgt_info,
cc     &       val_label=(/'F_HOLE'/))
c      end if
c dbg
c      call set_rule2('F_NORM',PRINT_FORMULA,tgt_info)
c      call set_arg('F_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORM'/))
c dbgend
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

      ! metric times amplitude vector
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'F_SC0'
      labels(2) = 'F_NORM'
      labels(3) = 'OMGred'
      labels(4) = 'Ex^+'
      labels(5) = ' '
      call add_target('F_SC0',ttype_frm,.false.,tgt_info)
      call set_dependency('F_SC0','F_NORM',tgt_info)
      call set_dependency('F_SC0','OMGred',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'precursor for metric matrix',1,'---')
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
      labels(4) = 'Ex'
      labels(5) = ' '
      call add_target('F_D',ttype_frm,.false.,tgt_info)
      call set_dependency('F_D','F_SC0',tgt_info)
      call set_dependency('F_D','D',tgt_info)
      call form_parameters(-1,
     &     parameters,2,'Metric matrix in the active space',1,'---')
      call set_rule('F_D',ttype_frm,DERIVATIVE,
     &              labels,5,1,
     &              parameters,2,tgt_info)
c      ! delete terms which do not have any open lines (valence space)
c      ! => scalar occupation class is excluded in the definition
c      labels(2:20)(1:len_target_name) = ' '
c      labels(2) = 'F_D'
c      labels(3) = 'D'
c      labels(4) = '1'
cc      labels(5) = 'C0'
c      labels(5) = 'DENS'
c      if (gno.eq.1) then
c        labels(4) = 'HOLE'
c        labels(5) = 'CUM'
c      end if
c      call form_parameters(-1,
c     &       parameters,2,'D formula',3,'ext')
c      call set_rule('F_D',ttype_frm,SELECT_LINE,
c     &              labels,5,1,
c     &              parameters,2,tgt_info)
c dbg
c      call set_rule2('F_D',KEEP_TERMS,tgt_info)
c      call set_arg('F_D',KEEP_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_D'/))
c      call set_arg('F_D',KEEP_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_D'/))
c      call set_arg('F_D',KEEP_TERMS,'TERMS',2,tgt_info,
c     &     val_int=(/6,9/))
c dbgend
c dbg
c      call form_parameters(-1,parameters,2,'stdout',1,'stdout')
c      call set_rule('F_D',ttype_frm,PRINT_FORMULA,
c     &                labels,2,1,parameters,2,tgt_info)
c      call set_rule2('F_D',ABORT,tgt_info)
c dbgend

*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! density matrix
      call add_target2('FOPT_D',.false.,tgt_info)
      call set_dependency('FOPT_D','F_D',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_D','DEF_ME_D',tgt_info)
      if (l_iccc) call set_dependency('FOPT_D','DEF_ME_1scal',tgt_info)
      call set_rule2('FOPT_D',OPTIMIZE,tgt_info)
      call set_arg('FOPT_D',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_D'/))
      if (orb_info%nactel.eq.0) then
        call set_arg('FOPT_D',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &               val_label=(/'F_D'/))
      else if (gno.eq.0) then
        call set_dependency('FOPT_D','DEF_ME_DENS',tgt_info)
       if (densmix.gt.0d0) then
        call set_dependency('FOPT_D','F_DENSmix',tgt_info)
        call set_arg('FOPT_D',OPTIMIZE,'LABELS_IN',2,tgt_info,
     &               val_label=(/'F_DENSmix','F_D      '/))
       else
        call set_dependency('FOPT_D','F_DENS0',tgt_info)
        call set_arg('FOPT_D',OPTIMIZE,'LABELS_IN',2,tgt_info,
     &               val_label=(/'F_DENS0','F_D    '/))
       end if
      else if (gno.eq.1) then
        call set_dependency('FOPT_D','DEF_ME_DENS',tgt_info)
        call set_dependency('FOPT_D','Y_GNO',tgt_info)
       if (densmix.gt.0d0) then
        call set_dependency('FOPT_D','F_DENSmix',tgt_info)
        call set_arg('FOPT_D',OPTIMIZE,'LABELS_IN',3,tgt_info,
     &               val_label=(/'F_DENSmix','F_Y_GNO  ','F_D      '/))
       else
        call set_dependency('FOPT_D','F_DENS0',tgt_info)
        call set_arg('FOPT_D',OPTIMIZE,'LABELS_IN',3,tgt_info,
     &               val_label=(/'F_DENS0','F_Y_GNO','F_D    '/))
       end if
      end if

*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

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
      call set_arg('DEF_ME_A',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &     val_int=(/n_states/))
      call set_arg('DEF_ME_A',DEF_ME_LIST,'REC',1,tgt_info,
     &     val_int=(/1/))
      if (prc_type.lt.3) then
        call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &       val_int=(/1/))
        call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &       val_int=(/1/))
        call set_arg('DEF_ME_A',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &       val_int=(/0/))
      end if

      ! ME_Auni
      call add_target2('DEF_ME_Auni',.false.,tgt_info)
      call set_dependency('DEF_ME_Auni','DEF_ME_A',tgt_info)
      call set_rule2('DEF_ME_Auni',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_Auni'/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'A'/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &     val_int=(/n_states/))
      call set_arg('DEF_ME_Auni',DEF_ME_LIST,'REC',1,tgt_info,
     &     val_int=(/1/))
      if (prc_type.lt.3) then
        call set_arg('DEF_ME_Auni',DEF_ME_LIST,'DIAG_TYPE',1,
     &       tgt_info,val_int=(/1/))
        call set_arg('DEF_ME_Auni',DEF_ME_LIST,'DIAG_IRREP',1,
     &       tgt_info,val_int=(/1/))
        call set_arg('DEF_ME_Auni',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &       val_int=(/0/))
      end if
      ! reassign operator to original list
      labels(1) = 'ME_A'
      labels(2) = 'A'
      call set_rule('DEF_ME_Auni',ttype_opme,ASSIGN_ME2OP,
     &     labels,2,1,
     &     parameters,0,tgt_info)

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
      call set_arg('DEF_ME_1',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_1',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call dens_parameters(-1,parameters,0,0,0)
      call set_rule('DEF_ME_1',ttype_opme,UNITY,
     &     'ME_1',1,1,
     &     parameters,1,tgt_info)

      ! ME_1scal
      call add_target2('DEF_ME_1scal',.false.,tgt_info)
      call set_dependency('DEF_ME_1scal','1scal',tgt_info)
      call set_rule2('DEF_ME_1scal',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_1scal'/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'1scal'/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_1scal',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/1/))
      call dens_parameters(-1,parameters,0,0,0)
      call set_rule('DEF_ME_1scal',ttype_opme,UNITY,
     &     'ME_1scal',1,1,
     &     parameters,1,tgt_info)

      ! ME_E(MR)
      call add_target('DEF_ME_E(MR)',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_E(MR)','E(MR)',tgt_info)
      do i_state = 0,n_states
      if (i_state.gt.0.and..not.multistate) cycle
      c_st = state_label(i_state,.true.)
      if(i_state.eq.0) c_st = ""
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_E(MR)'//trim(c_st)
      labels(2) = 'E(MR)'//trim(c_st)
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DEF_ME_E(MR)',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      end do

      ! ME_D
      call add_target2('DEF_ME_D',.false.,tgt_info)
      call set_dependency('DEF_ME_D','D',tgt_info)
      call set_rule2('DEF_ME_D',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_D',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_D'/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'D'/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &     val_int=(/msc/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &     val_int=(/n_states/))
      call set_arg('DEF_ME_D',DEF_ME_LIST,'REC',1,tgt_info,
     &     val_int=(/1/))
      ! automatic enforcement of S2 only without GNO
      ! In case of GNO, this is done manually, because cumulant-based 
      ! overlap tensors have different spin symmetrization rules
      if (spinadapt.ge.2.and.gno.eq.0)
     &   call set_arg('DEF_ME_D',DEF_ME_LIST,'S2',1,tgt_info,
     &        val_int=(/0/))

      ! inverted ME_D
      call add_target('DEF_ME_Dinv',ttype_opme,svdonly,tgt_info)
      ! invsqrt assumes 'T' as cluster operator: important in print
      if (l_iccc)
     &   call set_dependency('DEF_ME_Dinv','T',tgt_info)
      if (l_icci.and.svdonly)
     &   call set_dependency('DEF_ME_Dinv','C',tgt_info)
      call set_dependency('DEF_ME_Dinv','EVAL_D',tgt_info)
      call set_dependency('DEF_ME_Dinv','DEF_ME_D',tgt_info)
      call set_rule2('DEF_ME_Dinv',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dinv'/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'D'/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'S2',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &             val_log=(/.false./))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_Dinv',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))
      do i_state=1,n_states,1
      call set_rule2('DEF_ME_Dinv',INVERT,tgt_info)
      if (gno.eq.1.and.project.eq.1) then
        call set_dependency('DEF_ME_Dinv','EVAL_GNOSO',tgt_info)
        call set_arg('DEF_ME_Dinv',INVERT,'LIST_INV',2,tgt_info,
     &       val_label=(/'ME_D    ','ME_GNOSO'/))
      else
        call set_arg('DEF_ME_Dinv',INVERT,'LIST_INV',1,tgt_info,
     &       val_label=(/'ME_D'/))
      end if
      if (prc_type.eq.2) then
        call set_dependency('DEF_ME_Dinv','DEF_ME_Dunit',tgt_info)
        call set_arg('DEF_ME_Dinv',INVERT,'LIST',2,tgt_info,
     &       val_label=(/'ME_Dinv ','ME_Dunit'/))
      else
        call set_arg('DEF_ME_Dinv',INVERT,'LIST',1,tgt_info,
     &       val_label=(/'ME_Dinv'/))
      end if
      if (l_icci) then
        call set_arg('DEF_ME_Dinv',INVERT,'MODE',1,tgt_info,
     &       val_str='invsqrthalf')
      else !for MRCC, we also need the projector matrix
        call set_arg('DEF_ME_Dinv',INVERT,'MODE',1,tgt_info,
     &       val_str='invsqrt')
      end if
      if (gno.gt.0.and.l_iccc) then
        ! add terms needed for trafo to GNO basis
        call set_dependency('DEF_ME_Dinv','FOPT_Dinv_GNO',tgt_info)
        call set_rule2('DEF_ME_Dinv',EVAL,tgt_info)
        call set_arg('DEF_ME_Dinv',EVAL,'FORM',1,tgt_info,
     &               val_label=(/'FOPT_Dinv_GNO'/))
        call set_arg('DEF_ME_Dinv',EVAL,'INIT',1,tgt_info,
     &               val_log=(/.false./))
      end if
      labels(1) = 'ME_D'
      labels(2) = 'D'
      call set_rule('DEF_ME_Dinv',ttype_opme,ASSIGN_ME2OP,
     &     labels,2,1,
     &     parameters,0,tgt_info)
      if (gno.gt.0.and.l_iccc) then
        ! add terms needed for Projector using GNO basis
        call set_dependency('DEF_ME_Dinv','FOPT_Dproj_GNO',tgt_info)
        call set_rule2('DEF_ME_Dinv',EVAL,tgt_info)
        call set_arg('DEF_ME_Dinv',EVAL,'FORM',1,tgt_info,
     &               val_label=(/'FOPT_Dproj_GNO'/))
        call set_arg('DEF_ME_Dinv',EVAL,'INIT',1,tgt_info,
     &               val_log=(/.false./))
      end if
c dbg
c      call set_rule2('DEF_ME_Dinv',PRINT_MEL_INFO_,tgt_info)
c      call set_arg('DEF_ME_Dinv',PRINT_MEL_INFO_,'LIST',1,tgt_info,
c     &               val_label=(/'ME_Dinv'/))
c      call form_parameters(-1,parameters,2,
c     &     'Square root of inverse density matrix :',0,'LIST')
c      labels(1) = 'DEF_ME_Dinv'
c      labels(2) = 'ME_Dinv'
c      call set_rule('DEF_ME_Dinv',ttype_opme,PRINT_MEL,
c     &     'ME_Dinv',1,0,
c     &     parameters,2,tgt_info)
c dbgend
       if(multistate) then
        call set_rule2('DEF_ME_Dinv',ADV_STATE,tgt_info)
        call set_arg('DEF_ME_Dinv',ADV_STATE,'LISTS',1,tgt_info,
     &       val_label=['ME_Dinv'])
        call set_arg('DEF_ME_Dinv',ADV_STATE,'N_ROOTS',1,tgt_info,
     &       val_int=[n_states])
        call set_rule2('DEF_ME_Dinv',ADV_STATE,tgt_info)
        call set_arg('DEF_ME_Dinv',ADV_STATE,'LISTS',1,tgt_info,
     &       val_label=['ME_D'])
        call set_arg('DEF_ME_Dinv',ADV_STATE,'N_ROOTS',1,tgt_info,
     &       val_int=[n_states])
       end if
      end do

      ! reordered inverted ME_D
      call add_target('DEF_ME_Dtr',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dtr','Dtr',tgt_info)
      call set_rule2('DEF_ME_Dtr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dtr'/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Dtr'/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'S2',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &             val_log=(/.false./))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_Dtr',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))
      call set_dependency('DEF_ME_Dtr','DEF_ME_Dinv',tgt_info)
      do i_state = 1,n_states
       c_st = state_label(i_state,.false.)
       labels(1:20)(1:len_target_name) = ' '
       labels(1) = 'ME_Dtr'
       labels(2) = 'Dtr'
       labels(2) = 'ME_Dinv'
       call form_parameters(-1,parameters,2,
     &      '---',13,'---')
       call set_rule('DEF_ME_Dtr',ttype_opme,
     &      REORDER_MEL,
     &      labels,2,1,
     &      parameters,2,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Reordered inverted Density matrix :',0,'LIST')
c      call set_rule('DEF_ME_Dtr',ttype_opme,PRINT_MEL,
c     &     'ME_Dtr',1,0,
c     &     parameters,2,tgt_info)
c      if(i_state.eq.n_states)
c     &     call set_rule2('DEF_ME_Dtr',ABORT,tgt_info)
c dbgend
       if(multistate) then
        call set_rule2('DEF_ME_Dtr',ADV_STATE,tgt_info)
        call set_arg('DEF_ME_Dtr',ADV_STATE,'LISTS',2,tgt_info,
     &       val_label=['ME_Dinv','ME_Dtr '])
        call set_arg('DEF_ME_Dtr',ADV_STATE,'N_ROOTS',1,tgt_info,
     &       val_int=[n_states])
       end if
      end do

      ! reordered daggered inverted ME_D
      call add_target('DEF_ME_Dtrdag',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dtrdag','DEF_ME_Dtr',tgt_info)
      if (l_iccc.and.prc_traf)
     &   call set_dependency('DEF_ME_Dtrdag','EVAL_Atr',tgt_info)
      call set_rule2('DEF_ME_Dtrdag',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dtrdag'/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Dtr'/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'S2',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &             val_log=(/.false./))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_Dtrdag',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))
      ! labels for advance state
      labels(1) = 'ME_Dtrdag'
      labels(2) = 'ME_Dtr'
      idef = 2
      if (prc_type.eq.2) then
       idef = idef+1
       labels(idef) = 'ME_Dunit'
      else
       idef = idef+1
       labels(idef) = 'ME_Dinv'
      end if
      call set_dependency('DEF_ME_Dtrdag','DEF_ME_Dinv',tgt_info)
      do i_state = 1, n_states
      c_st = state_label(i_state,.false.)
      call set_rule2('DEF_ME_Dtrdag',REORDER_MEL,tgt_info)
      call set_arg('DEF_ME_Dtrdag',REORDER_MEL,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_Dtrdag'/))
      if (prc_type.eq.2) then
        ! misuse ME_Dtrdag as intermediate in preconditioner calculation
        call set_arg('DEF_ME_Dtrdag',REORDER_MEL,'LIST_IN',1,tgt_info,
     &               val_label=(/'ME_Dunit'/))
      else
        call set_arg('DEF_ME_Dtrdag',REORDER_MEL,'LIST_IN',1,tgt_info,
     &               val_label=(/'ME_Dinv'/))
      end if
      call set_arg('DEF_ME_Dtrdag',REORDER_MEL,'FROMTO',1,tgt_info,
     &             val_int=(/13/))
      call set_arg('DEF_ME_Dtrdag',REORDER_MEL,'ADJOINT',1,tgt_info,
     &             val_log=(/.true./))
      ! reassign operator to original list
      call set_rule2('DEF_ME_Dtrdag',ASSIGN_ME2OP,tgt_info)
      call set_arg('DEF_ME_Dtrdag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('DEF_ME_Dtrdag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Reordered transposed inverted Density matrix :',0,'LIST')
c      call set_rule('DEF_ME_Dtrdag',ttype_opme,PRINT_MEL,
c     &     'ME_Dtrdag',1,0,
c     &     parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'ME_Dtr :',0,'LIST')
c      call set_rule('DEF_ME_Dtrdag',ttype_opme,PRINT_MEL,
c     &     'ME_Dtr',1,0,
c     &     parameters,2,tgt_info)
c      if(i_state.eq.n_states)
c     &      call set_rule2('DEF_ME_Dtrdag',ABORT,tgt_info)
c dbgend
       if(multistate) then
        call set_rule2('DEF_ME_Dtrdag',ADV_STATE,tgt_info)
        call set_arg('DEF_ME_Dtrdag',ADV_STATE,'LISTS',idef,tgt_info,
     &       val_label=labels(1:idef))
        call set_arg('DEF_ME_Dtrdag',ADV_STATE,'N_ROOTS',1,tgt_info,
     &       val_int=[n_states])
       end if
      end do

      ! ME_Dunit
      call add_target('DEF_ME_Dunit',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dunit','D',tgt_info)
      call set_rule2('DEF_ME_Dunit',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dunit'/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'D'/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'S2',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &             val_log=(/.false./))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &             val_int=(/n_states/))
      call set_arg('DEF_ME_Dunit',DEF_ME_LIST,'REC',1,tgt_info,
     &             val_int=(/1/))

      ! reordered daggered inverted ME_D, squared element-wise
      call add_target('DEF_ME_Dudag_2',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_Dudag_2','DEF_ME_Dtrdag',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dudag_2'
      labels(2) = 'Dtr'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule('DEF_ME_Dudag_2',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)
      call set_rule2('DEF_ME_Dudag_2',SCALE_COPY,tgt_info)
      call set_arg('DEF_ME_Dudag_2',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_Dudag_2'/))
      call set_arg('DEF_ME_Dudag_2',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_Dtrdag'/))
      call set_arg('DEF_ME_Dudag_2',SCALE_COPY,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))
      call set_arg('DEF_ME_Dudag_2',SCALE_COPY,'MODE',1,tgt_info,
     &             val_str='square')
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Dtrdag, squared element-wise :',0,'LIST')
c      labels(1) = 'DEF_ME_Dudag_2'
c      labels(2) = 'ME_Dudag_2'
c      call set_rule('DEF_ME_Dudag_2',ttype_opme,PRINT_MEL,
c     &     'ME_Dudag_2',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME_NORM
      call add_target('DEF_ME_NORM',ttype_opme,.false.,tgt_info)
      call set_dependency('DEF_ME_NORM','NORM',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_NORM'
      labels(2) = 'NORM'
      call me_list_parameters(-1,parameters,
     &     0,0,1,0,0,.false.)
      call set_rule('DEF_ME_NORM',ttype_opme,DEF_ME_LIST,
     &              labels,2,1,
     &              parameters,1,tgt_info)

      ! List & form. needed for sequential orthogonalization within GNO
      ! a) set up operator and list
      call add_target2('DEF_ME_GNOSO',.false.,tgt_info)
      call set_rule2('DEF_ME_GNOSO',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('DEF_ME_GNOSO',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'OP_GNOSO'/))
      occ_def = 0
      occ_def(IVALE,1:2,1:2) = 1
      call set_arg('DEF_ME_GNOSO',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('DEF_ME_GNOSO',DEF_OP_FROM_OCC,'OCC',2,tgt_info,
     &             val_occ=occ_def)
      call set_rule2('DEF_ME_GNOSO',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_GNOSO',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'OP_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_GNOSO',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('DEF_ME_GNOSO',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &     val_int=(/msc/))
      ! b) formula: minus gamma_1 for both blocks
      call set_dependency('DEF_ME_GNOSO','C0',tgt_info)
      call set_rule2('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'OP_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'OP_GNOSO','C0^+    ','C0      ','OP_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/1,2,3,1/))
      call set_arg('DEF_ME_GNOSO',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/-1d0/))
c dbg
c      call set_rule2('DEF_ME_GNOSO',PRINT_FORMULA,tgt_info)
c      call set_arg('DEF_ME_GNOSO',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_GNOSO'/))
c dbgend
      ! c) optimize
      call set_dependency('DEF_ME_GNOSO','DEF_ME_C0',tgt_info)
      call set_rule2('DEF_ME_GNOSO',OPTIMIZE,tgt_info)
      call set_arg('DEF_ME_GNOSO',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_GNOSO'/))
      call set_arg('DEF_ME_GNOSO',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_GNOSO'/))

*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Evaluate density matrix
      call add_target('EVAL_D',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_D','FOPT_D',tgt_info)
      call set_dependency('EVAL_D','EVAL_REF_S(S+1)',tgt_info)
      if (gno.gt.0.and.l_iccc)
     &     call set_dependency('EVAL_D','FOPT_D_GNO',tgt_info)
      do i_state = 1,n_states
      c_st = state_label(i_state,.false.)
      call set_rule('EVAL_D',ttype_opme,EVAL,
     &     'FOPT_D',1,0,
     &     parameters,0,tgt_info)
      if (gno.gt.0.and.l_iccc) then
        ! perform spin projection? (not done automatically for GNO)
        if (spinadapt.ge.2) then
          call set_rule2('EVAL_D',SPIN_PROJECT,tgt_info)
          call set_arg('EVAL_D',SPIN_PROJECT,'LIST',1,tgt_info,
     &                 val_label=(/'ME_D'/))
          call set_arg('EVAL_D',SPIN_PROJECT,'S2',1,tgt_info,
     &                 val_int=(/0/))
        end if
        ! transform to GNO basis
        call set_rule2('EVAL_D',EVAL,tgt_info)
        call set_arg('EVAL_D',EVAL,'FORM',1,tgt_info,
     &               val_label=(/'FOPT_D_GNO'/))
        call set_arg('EVAL_D',EVAL,'INIT',1,tgt_info,
     &               val_log=(/.false./))
      end if
c      ! fix: set first element (zero occ) to 1 (had not been defined in F_D)
c dbg
c      labels(1:10)(1:len_target_name) = ' '
c      labels(1) = 'ME_D'
c      call dens_parameters(-1,parameters,1,1,1)
c      call set_rule('EVAL_D',ttype_opme,UNITY,
c     &     labels,1,1,
c     &     parameters,1,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'Density matrix :',0,'LIST')
c      call set_rule('EVAL_D',ttype_opme,PRINT_MEL,
c     &     'ME_D',1,0,
c     &     parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'RDMs :',0,'LIST')
c      call set_rule('EVAL_D',ttype_opme,PRINT_MEL,
c     &     'ME_DENS',1,0,
c     &     parameters,2,tgt_info)
c      call set_rule2('EVAL_D',ABORT,tgt_info)
c dbgend
      if(multistate)then
       call set_rule2('EVAL_D',ADV_STATE,tgt_info)
       call set_arg('EVAL_D',ADV_STATE,'LISTS',1,tgt_info,
     &      val_label=['ME_C0'])
       call set_arg('EVAL_D',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
       call set_rule2('EVAL_D',ADV_STATE,tgt_info)
       call set_arg('EVAL_D',ADV_STATE,'LISTS',1,tgt_info,
     &      val_label=['ME_DENS'])
       call set_arg('EVAL_D',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
       call set_rule2('EVAL_D',ADV_STATE,tgt_info)
       call set_arg('EVAL_D',ADV_STATE,'LISTS',1,tgt_info,
     &      val_label=['ME_D'])
       call set_arg('EVAL_D',ADV_STATE,'N_ROOTS',1,tgt_info,
     &      val_int=[n_states])
      end if
      end do

      ! eval list needed for sequential orthogonalization within GNO
      ! a) evaluate
      call add_target2('EVAL_GNOSO',.false.,tgt_info)
      call set_dependency('EVAL_GNOSO','DEF_ME_GNOSO',tgt_info)
      call set_rule2('EVAL_GNOSO',EVAL,tgt_info)
      call set_arg('EVAL_GNOSO',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_GNOSO'/))
      ! b) substract unity from second block (gives minus hole density)
      call set_rule2('EVAL_GNOSO',UNITY,tgt_info)
      call set_arg('EVAL_GNOSO',UNITY,'LIST',1,tgt_info,
     &             val_label=(/'ME_GNOSO'/))
      call set_arg('EVAL_GNOSO',UNITY,'FAC',1,tgt_info,
     &             val_rl8=(/-1d0/))
      call set_arg('EVAL_GNOSO',UNITY,'MIN_BLK',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('EVAL_GNOSO',UNITY,'MAX_BLK',1,tgt_info,
     &             val_int=(/2/))
      ! c) invert
      call set_rule2('EVAL_GNOSO',INVERT,tgt_info)
      call set_arg('EVAL_GNOSO',INVERT,'LIST_INV',1,tgt_info,
     &             val_label=(/'ME_GNOSO'/))
      call set_arg('EVAL_GNOSO',INVERT,'LIST',1,tgt_info,
     &             val_label=(/'ME_GNOSO'/))
      call set_arg('EVAL_GNOSO',INVERT,'MODE',1,tgt_info,
     &             val_str='pseudoinv')
c dbg
c      call set_rule2('EVAL_GNOSO',PRINT_MEL,tgt_info)
c      call set_arg('EVAL_GNOSO',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_GNOSO'/))
c dbgend

      return
      end
