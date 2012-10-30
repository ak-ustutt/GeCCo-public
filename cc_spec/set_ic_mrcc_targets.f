*----------------------------------------------------------------------*
      subroutine set_ic_mrcc_targets(tgt_info,orb_info,
     &                               excrestr,maxh,maxp,execute)
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

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     maxh, maxp, excrestr(0:maxh,0:maxp,1:2)
      logical, intent(in) ::
     &     execute

      integer ::
     &     ndef, occ_def(ngastp,2,124),!60),
     &     icnt, maxexc, maxv,
     &     msc, ip, ih, ivv, iv, ivv2, jvv,
     &     maxcom, maxcom_en, maxcom_h1bar, h1bar_maxp,
     &     i_cls, len_form, optref, idef, ciroot, maxroot,
     &     stndT(2,60), stndD(2,60), nsupT, nsupD,
     &     G_level, iexc, jexc, maxtt, iblk, jblk, kblk, prc_type,
     &     tred, nremblk, remblk(60), igasreo(3), ngas, lblk, ntrunc,
     &     tfix, maxit, t1ord, maxcum, cum_appr_mode, gno
      logical ::
     &     update_prc, skip, preopt, project, first, Op_eqs,
     &     h1bar, htt, svdonly, fact_tt, ex_t3red, trunc, l_exist,
     &     oldref, solve, use_f12, restart
      character(len_target_name) ::
     &     dia_label, dia_label2,
     &     labels(20)
      character(len_command_par) ::
     &     parameters(3)
      character(len=64) ::
     &     descr, descr_h1, descr_h2pp, descr_h2ppp
      character ::
     &     op_ht*3, f_ht*5, op_ht0to*6, f_ht0to*8, form_str*50,
     &     def_ht*10
      character(len=3) ::
     &     prc_mode_str
      real(8) ::
     &     x_ansatz, prc_shift

      if (iprlvl.gt.0) write(luout,*) 'setting icMRCC targets'

      ! CAVEAT: should be adapted as soon as open-shell version
      !         is up and running
      msc = +1
      if (orb_info%ims.ne.0) msc = 0

      ! get some keywords
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      call get_argument_value('method.MR','ciroot',
     &     ival=ciroot)
      call get_argument_value('method.MR','maxroot',
     &     ival=maxroot)
      if(maxroot.le.0) maxroot=ciroot
      call get_argument_value('method.MR','prc_type',
     &     ival=prc_type)
      call get_argument_value('method.MR','prc_shift',
     &     xval=prc_shift)
      call get_argument_value('method.MR','svdonly',
     &     lval=svdonly)
      call get_argument_value('calculate.solve.non_linear','optref',
     &     ival=optref)
      call get_argument_value('calculate.solve.non_linear','update_prc',
     &     lval=update_prc)
      call get_argument_value('calculate.solve.non_linear','preopt',
     &     lval=preopt)
      call get_argument_value('calculate.solve.non_linear','restart',
     &     lval=restart)
      call get_argument_value('method.MR','project',
     &     lval=project)
      call get_argument_value('method.MRCC','Op_eqs',
     &     lval=Op_eqs)
      call get_argument_value('method.MRCC','maxcom_res',
     &     ival=maxcom)
      call get_argument_value('method.MRCC','maxcom_en',
     &     ival=maxcom_en)
      call get_argument_value('method.MRCC','H1bar',
     &     lval=h1bar)
      call get_argument_value('method.MRCC','maxcom_h1bar',
     &     ival=maxcom_h1bar)
      call get_argument_value('method.MRCC','h1bar_maxp',
     &     ival=h1bar_maxp)
      if (h1bar_maxp.lt.0.and.maxexc.le.2) h1bar_maxp = 2
      if (h1bar_maxp.lt.0.and.maxexc.gt.2) h1bar_maxp = 3
      if (.not.h1bar) h1bar_maxp = 4
      call get_argument_value('method.MRCC','G_level',
     &     ival=G_level)
      call get_argument_value('method.MRCC','HTT',
     &     lval=htt)
      call get_argument_value('method.MRCC','maxtt',
     &     ival=maxtt)
      call get_argument_value('method.MRCC','x_ansatz',
     &     xval=x_ansatz)
      call get_argument_value('method.MRCC','Tred_mode',
     &     ival=tred)
      call get_argument_value('method.MRCC','trunc_order',
     &     ival=ntrunc)
      call get_argument_value('method.MRCC','Tfix',
     &     ival=tfix)
      call get_argument_value('method.MRCC','T1ord',
     &     ival=t1ord)
      call get_argument_value('method.MR','oldref',
     &     lval=oldref)
      call get_argument_value('method.MR','maxcum',
     &     ival=maxcum)
      call get_argument_value('method.MR','cum_appr_mode',
     &     ival=cum_appr_mode)
      call get_argument_value('calculate.solve','maxiter',
     &     ival=maxit)
      if (is_argument_set('calculate.solve.non_linear','maxiter').gt.0)
     &     call get_argument_value('calculate.solve.non_linear',
     &     'maxiter',ival=maxit)
      call get_argument_value('method.MR','maxv',
     &     ival=maxv)
      skip = (is_keyword_set('calculate.skip_E').gt.0)
      if (maxv.lt.0) maxv = 2*maxexc
      trunc = ntrunc.ge.0
      solve = execute.and..not.svdonly.and.(tfix.eq.0.or.maxit.gt.1)
     &               .and..not.skip
      use_f12 = is_keyword_set('method.R12').gt.0
      call get_argument_value('method.MR','GNO',
     &     ival=gno)

      if (ntest.ge.100) then
        write(luout,*) 'maxcom_en    = ', maxcom_en
        write(luout,*) 'maxcom_res   = ', maxcom
        write(luout,*) 'G_level      = ', G_level
        write(luout,*) 'preopt       = ', preopt
        write(luout,*) 'Op_eqs       = ', Op_eqs
        write(luout,*) 'H1bar        = ', h1bar
        if (h1bar) write(luout,*) 'maxcom_h1bar = ', maxcom_h1bar
        if (h1bar) write(luout,*) 'h1bar_maxp   = ', h1bar_maxp
        write(luout,*) 'HTT          = ', htt
        write(luout,*) 'maxtt        = ', maxtt
        write(luout,*) 'x_ansatz     = ', x_ansatz
        write(luout,*) 'Tred_mode    = ', tred
        write(luout,*) 'trunc        = ', trunc
        if (tfix.gt.0) write(luout,*) 'Tfix         = ', tfix
        if (t1ord.ge.0) write(luout,*) 'T1ord        = ', t1ord
      end if

      if (x_ansatz.ne.0.5d0.and.x_ansatz.ne.0d0.and.abs(x_ansatz).ne.1d0
     &    .and.x_ansatz.ne.-2d0.and.x_ansatz.ne.-3d0
     &    .and.(maxcom_en.gt.2.or.maxcom.gt.2
     &          .or.h1bar.and.maxcom_h1bar.gt.2))
     &    call quit(1,'set_ic_mrcc_targets',
     &      'x_ansatz.ne.0/0.5/1 currently works only with Ncom<=2')
      if (tred.gt.0.and.optref.eq.0)
     &    call quit(1,'set_ic_mrcc_targets',
     &     'Tred_mode > 0 not yet available for optref=0')
      if (tfix.gt.0.and.(.not.oldref.or..not.project.or.optref.ne.0))
     &    call quit(1,'set_ic_mrcc_targets',
     &     'Tfix>0 only allowed with oldref=T,project=T,optref=0')
      if (t1ord.ge.0.and.tfix.eq.0)
     &    call quit(1,'set_ic_mrcc_targets',
     &     'Manually setting T1ord only enabled yet for Tfix>0')
      if (h1bar_maxp.eq.0.or.h1bar_maxp.eq.1)
     &    call quit(1,'set_ic_mrcc_targets',
     &     'h1bar_maxp should be either -1 or one of 2,3,4')
      if (gno.gt.0.and.(update_prc.or.tred.gt.0))
     &    call quit(1,'set_ic_mrcc_targets',
     &     'update_prc or tred!=0 not yet available for GNO')
      
*----------------------------------------------------------------------*
*     Operators:
*----------------------------------------------------------------------*

      ! define particle conserving deexcitation operator L (for icMRCC)
      call add_target('L',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ndef = ndef + 1
            occ_def(IHOLE,1,ndef) = ih
            occ_def(IPART,2,ndef) = ip
            occ_def(IVALE,2,ndef) = iexc - ip
            occ_def(IVALE,1,ndef) = iexc - ih
          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('L',ttype_op,DEF_OP_FROM_OCC,
     &              'L',1,1,
     &              parameters,2,tgt_info)

      ! define particle conserving excitation operator T (for icMRCC)
      call add_target('T',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      nsupT=0
      do ip = 0, maxp
        do ih = 0, maxh
          first = .true.
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
            if (first) then
              nsupT = nsupT + 1
              stndT(1,nsupT) = ndef
              first = .false.
            end if
            stndT(2,nsupT) = ndef
          end do
        end do
      end do
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

      ! we need to set "super-block" information of metric
      ! since metric is defined in set_ic_mrci_targets,
      ! we just run over the same loops here:
      ndef = 0
      nsupD = 0
      nremblk = 0
      do ip = 0, maxp
        do ih = 0, maxh
          first = .true.
          do iexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
            do jexc = excrestr(ih,ip,2), excrestr(ih,ip,1),-1
            if ((project.or.Op_eqs).and.iexc.ne.jexc) cycle
            ! not for purely inactive excitation class
            if (ip.eq.ih.and.
     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))) cycle
            ndef = ndef + 1
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
      if (nsupD.gt.nsupT) call quit(1,'set_ic_mrcc_targets',
     &       'more super-blocks for metric than for T?')

      ! TT intermediate
      call add_target('TT',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do iexc = excrestr(0,ip,1), excrestr(0,ip,2)
          do ih = 0, maxh
           do jexc = excrestr(ih,0,1), excrestr(ih,0,2)
            if (iexc.ne.jexc) cycle ! only T1-T1, T2-T2, etc.
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
          end do
         end do
        end do
      end do
      fact_tt = ndef.gt.0.and..false.
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

      ! generic T-shaped output operator for transformations
      call add_target2('Tout',.false.,tgt_info)
      call set_dependency('Tout','T',tgt_info)
      call set_rule2('Tout',CLONE_OP,tgt_info)
      call set_arg('Tout',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Tout'/))
      call set_arg('Tout',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'T'/))

      ! transformed L
      call add_target2('Ltr',.false.,tgt_info)
      call set_dependency('Ltr','L',tgt_info)
      call set_rule2('Ltr',CLONE_OP,tgt_info)
      call set_arg('Ltr',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Ltr'/))
      call set_arg('Ltr',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'L'/))

      ! define Residual
      call add_target('OMG',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = max(excrestr(ih,ip,1),1), excrestr(ih,ip,2)
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
      call set_rule('OMG',ttype_op,DEF_OP_FROM_OCC,
     &              'OMG',1,1,
     &              parameters,2,tgt_info)

c dbg
c      ! Metric operator (for testing, also non-diagonal blocks)
c      call add_target('S',ttype_op,.false.,tgt_info)
c      occ_def = 0
c      ndef = 0
c      do ip = 0, maxp
c        do ih = 0, maxh
c          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
c            do jexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
c            ! not for purely inactive excitation class
c            if (ip.eq.ih.and.
c     &          ip.eq.maxval(excrestr(0:maxh,0:maxp,2))) cycle
c            ndef = ndef + 1
c            occ_def(IHOLE,1,3*ndef-1) = ih
c            occ_def(IHOLE,2,3*ndef-1) = ih
c            occ_def(IPART,1,3*ndef-1) = ip
c            occ_def(IPART,2,3*ndef-1) = ip
c            occ_def(IVALE,1,3*ndef-1) = jexc - ip
c            occ_def(IVALE,2,3*ndef-1) = iexc - ip
c            occ_def(IVALE,2,3*ndef-2) = jexc - ih
c            occ_def(IVALE,1,3*ndef) = iexc - ih
c           end do
c          end do
c        end do
c      end do
c      call op_from_occ_parameters(-1,parameters,2,
c     &              occ_def,ndef,3,(/0,0,0,0,0,0/),ndef)
c      call set_rule('S',ttype_op,DEF_OP_FROM_OCC,
c     &              'S',1,1,
c     &              parameters,2,tgt_info)
c dbgend

      ! Diagonal Preconditioner
      call add_target(trim(op_dia)//'_T',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(trim(op_dia)//'_T','T',tgt_info)
      call cloneop_parameters(-1,parameters,'T',.false.)
      call set_rule(trim(op_dia)//'_T',ttype_op,CLONE_OP,
     &              trim(op_dia)//'_T',1,1,
     &              parameters,1,tgt_info)

      ! Diagonal Preconditioner for L
      call add_target(trim(op_dia)//'_L',ttype_op,.false.,
     &                tgt_info)
      call set_dependency(trim(op_dia)//'_L','L',tgt_info)
      call cloneop_parameters(-1,parameters,'L',.false.)
      call set_rule(trim(op_dia)//'_L',ttype_op,CLONE_OP,
     &              trim(op_dia)//'_L',1,1,
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

      ! define effective hamiltonian
      call add_target('Heff',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 0, orb_info%nactel
        ndef = ndef + 1
        occ_def(IVALE,1,ndef) = iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('Heff',ttype_op,DEF_OP_FROM_OCC,
     &              'Heff',1,1,
     &              parameters,2,tgt_info)

      ! Effective Hamiltonian in the excitation space
      call add_target2('Geff',.false.,tgt_info)
      call set_dependency('Geff','T',tgt_info)
      call set_rule2('Geff',CLONE_OP,tgt_info)
      call set_arg('Geff',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'Geff'/))
      call set_arg('Geff',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'T'/))
c dbg
c      call add_target('Geff',ttype_op,.false.,tgt_info)
c      occ_def = 0
c      ndef = 0
c      do ip = 0, maxexc
c        do ih = 0, min(maxexc,2) !only one hole orbital
c          iv = 2 !CAS(2,2)
c          if (max(ip,ih).eq.0) iv = -1
c          if (ip.eq.ih+1) iv = 1
c          if (ip.eq.ih+2) iv = 0
c          if (ip.gt.ih+2) iv = -1
c          do ivv = 0, iv
c            ndef = ndef + 1
c            occ_def(IHOLE,2,ndef) = ih
c            occ_def(IPART,1,ndef) = ip
c            occ_def(IVALE,1,ndef) = max(ih-ip,0) + ivv
c            occ_def(IVALE,2,ndef) = max(ip-ih,0) + ivv
c          end do
c        end do
c      end do
c      call op_from_occ_parameters(-1,parameters,2,
c     &              occ_def,ndef,1,(/0,0/),ndef)
c      call set_rule('Geff',ttype_op,DEF_OP_FROM_OCC,
c     &              'Geff',1,1,
c     &              parameters,2,tgt_info)
c dbgend

      ! T1 transformed Hamiltonian
      ! define manually using a complicated loop structure
      ! in order to sort all blocks with 3 or more P-lines to the end
      call add_target2('H1bar',.false.,tgt_info)
      call set_rule2('H1bar',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('H1bar',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &     val_label=(/'H1bar'/))
      igasreo(1:3) = (/IVALE,IHOLE,IPART/)
      ngas = 3
      if (orb_info%nactt_hpv(IHOLE).eq.0) then
        igasreo(2) = IPART
        ngas = 2
      end if
      occ_def = 0
      ndef = 1 ! scalar block
      ! one-electron part
      do iblk = 1, ngas
        do jblk = 1, ngas
          ndef = ndef + 1
          occ_def(igasreo(iblk),1,ndef) = 1
          occ_def(igasreo(jblk),2,ndef) = 1
        end do
      end do
      ! two-electron part
      first = .true.
      do iblk = 1, ngas ! C1
        do jblk = 1, ngas ! A1
          do kblk = 1, max(iblk,jblk) ! either C2 or A2
            do lblk = 1, min(kblk,min(iblk,jblk)) ! either A2 or C2
             if (kblk.le.iblk) then 
              ndef = ndef + 1
              occ_def(igasreo(iblk),1,ndef) = 1
              occ_def(igasreo(jblk),2,ndef) = 1
              occ_def(igasreo(kblk),1,ndef) = 
     &            occ_def(igasreo(kblk),1,ndef) + 1
              occ_def(igasreo(lblk),2,ndef) = 
     &            occ_def(igasreo(lblk),2,ndef) + 1
              ! do formal blocks start here?
              if (first) then
                icnt = 0
                if (iblk.eq.ngas) icnt = icnt + 1
                if (jblk.eq.ngas) icnt = icnt + 1
                if (kblk.eq.ngas) icnt = icnt + 1
                if (lblk.eq.ngas) icnt = icnt + 1
                if (icnt.gt.h1bar_maxp) then
                  first = .false.
                  call set_arg('H1bar',DEF_OP_FROM_OCC,'FORMAL',
     &                         1,tgt_info,val_int=(/ndef/))
                end if
              end if
             end if
             if (kblk.le.jblk.and.kblk.ne.lblk) then       
              ndef = ndef + 1
              occ_def(igasreo(iblk),1,ndef) = 1
              occ_def(igasreo(jblk),2,ndef) = 1
              occ_def(igasreo(kblk),2,ndef) =
     &            occ_def(igasreo(kblk),2,ndef) + 1
              occ_def(igasreo(lblk),1,ndef) =
     &            occ_def(igasreo(lblk),1,ndef) + 1
             end if
            end do
          end do
        end do
      end do
      call set_arg('H1bar',DEF_OP_FROM_OCC,'BLOCKS',1,tgt_info,
     &     val_int=(/ndef/))
      call set_arg('H1bar',DEF_OP_FROM_OCC,'OCC',ndef,tgt_info,
     &     val_occ=occ_def(1:ngastp,1:2,1:ndef))
c      call set_dependency('H1bar','H',tgt_info)
c      call set_rule2('H1bar',CLONE_OP,tgt_info)
c      call set_arg('H1bar',CLONE_OP,'LABEL',1,tgt_info,
c     &     val_label=(/'H1bar'/))
c      call set_arg('H1bar',CLONE_OP,'TEMPLATE',1,tgt_info,
c     &     val_label=(/'H'/))

      ! T1
      call add_target('T1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), min(excrestr(ih,ip,2),1)
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
      call set_rule('T1',ttype_op,DEF_OP_FROM_OCC,
     &              'T1',1,1,
     &              parameters,2,tgt_info)

      ! T - T1
      call add_target('T-T1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = max(2,excrestr(ih,ip,1)), excrestr(ih,ip,2)
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
      call set_rule('T-T1',ttype_op,DEF_OP_FROM_OCC,
     &              'T-T1',1,1,
     &              parameters,2,tgt_info)

      ! one-particle density (in a form digestable for preconditioner)
      call add_target('DENS1',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 1
      occ_def(IVALE,1,ndef) = 1
      occ_def(IVALE,2,ndef) = 1
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('DENS1',ttype_op,DEF_OP_FROM_OCC,
     &              'DENS1',1,1,
     &              parameters,2,tgt_info)

      ! a number of operators needed for H0inversion
      ! & Dyall type H0 (as also used for the diagonal preconditioner)
      call add_target3((/
     &     'target OPS_FOR_H0INV(           ',
     &     '    depend (OMG,T)              ',
     &     '    CLONE_OPERATOR(label=OMGprj,template=OMG) ',
     &     '    CLONE_OPERATOR(label=DlT,template=T)      ',
     &     '    CLONE_OPERATOR(label=H0_DlT,template=OMG) ',
     &     '    CLONE_OPERATOR(label=S_DlT,template=OMG) ',
     &     '    DEF_OP_FROM_OCC(label=D00,descr=",|,|,V|V,",join=2)',
     &     '    DEF_OP_FROM_OCC(label=H0Dy,descr="H,H|P,P|V,V|VV,VV")',
     &     ')'/),
     &     tgt_info)
    
      ! copies of selected PP blocks of Hamiltonian
      call add_target('H_PP',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 4-h1bar_maxp, 2
        ndef = ndef + 1
        occ_def(IPART,1,ndef) = 2
        occ_def(IPART,2,ndef) = 2-iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('H_PP',ttype_op,DEF_OP_FROM_OCC,
     &              'H_PP',1,1,
     &              parameters,2,tgt_info)

      call add_target('H_PP0',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 0, 2
        ndef = ndef + 1
        occ_def(IPART,1,ndef) = 2
        occ_def(IPART,2,ndef) = 2-iv
        occ_def(IVALE,2,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('H_PP0',ttype_op,DEF_OP_FROM_OCC,
     &              'H_PP0',1,1,
     &              parameters,2,tgt_info)

      ! copy of selected HH blocks of Hamiltonian
      call add_target('H_HH',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do iv = 0, 2
        ndef = ndef + 1
        occ_def(IHOLE,2,ndef) = 2
        occ_def(IHOLE,1,ndef) = 2-iv
        occ_def(IVALE,1,ndef) = iv
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('H_HH',ttype_op,DEF_OP_FROM_OCC,
     &              'H_HH',1,1,
     &              parameters,2,tgt_info)

      ! define intermediate for combining HH contractions
      call add_target('INT_HH',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
c      does not work for higher excitations like this
c      do ip = 2, maxp !only for blocks with at least two P lines
      do ip = 0, min(2,maxp)
        do ih = 0, maxh
c          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
c            ! for perturbative triples we should not define HH
c            if (trunc.and.iexc.gt.2) cycle
            iexc = 2  ! only doubles for now
            if (ip.eq.0.and.ih.eq.0) cycle ! purely act. not needed
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = iexc - ip
            occ_def(IVALE,2,ndef*2-1) = iexc - ih
c          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('INT_HH',ttype_op,DEF_OP_FROM_OCC,
     &              'INT_HH',1,1,
     &              parameters,2,tgt_info)

      ! define intermediate for combining PP contractions
      call add_target2('INT_PP0',.false.,tgt_info)
      call set_dependency('INT_PP0','INT_HH',tgt_info)
      call set_rule2('INT_PP0',CLONE_OP,tgt_info)
      call set_arg('INT_PP0',CLONE_OP,'LABEL',1,tgt_info,
     &     val_label=(/'INT_PP0'/))
      call set_arg('INT_PP0',CLONE_OP,'TEMPLATE',1,tgt_info,
     &     val_label=(/'INT_HH'/))

      ! define intermediate for combining PP contractions
      call add_target('INT_PP',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
c      does not work for higher excitations like this
c      do ip = 2, maxp !only for blocks with at least two P lines
      do ip = 0, min(min(2,maxp),h1bar_maxp-2)
        do ih = 0, maxh
c          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
c            ! for perturbative triples we should not define PP
c            if (trunc.and.iexc.gt.2) cycle
            iexc = 2  ! only doubles for now
            if (ip.eq.0.and.ih.eq.0) cycle ! purely act. not needed
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef*2) = ih
            occ_def(IPART,1,ndef*2) = ip
            occ_def(IVALE,1,ndef*2) = iexc - ip
            occ_def(IVALE,2,ndef*2-1) = iexc - ih
c          end do
        end do
      end do
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,2,(/0,0,0,0/),ndef)
      call set_rule('INT_PP',ttype_op,DEF_OP_FROM_OCC,
     &              'INT_PP',1,1,
     &              parameters,2,tgt_info)

      descr_h1 = 'H,[HV]|P,[HP]|V,[HPV]'
      descr_h2pp = 'HH,V[HP]|V[HV],[HVP][HVP]|[HV]P,[HV]P'
      descr_h2ppp = 'PP,[HV]P'
      descr = 'HV,PV'
 
      call add_target2('INT_HT2',.false.,tgt_info)
      call set_rule2('INT_HT2',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('INT_HT2',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'INT_HT2'/))
      call set_arg('INT_HT2',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)!_h1//'|'//descr_h2pp//'|'//descr_h2ppp)

      call add_target2('INT_T2H',.false.,tgt_info)
      call set_rule2('INT_T2H',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('INT_T2H',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'INT_T2H'/))
      call set_arg('INT_T2H',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str='PV,PH|PP,H[HPV]|PP,V[PV]')

      descr = ',|,|,V|V,'
      if (orb_info%nactel.ge.4) descr = ',|,|,V|V,|,VV|VV,'      
      if (orb_info%nactel.ge.6) 
     &       descr = trim(descr)//'|,VVV|VVV,|,VVVV|VVVV,|,VVVVV|VVVVV,' 
c dbg
c      print *,'nactel: ',orb_info%nactel
c      print *,'descr = ',trim(descr)
c dbg
      call add_target2('INT_D',.false.,tgt_info)
      call set_rule2('INT_D',DEF_OP_FROM_OCC,tgt_info)
      call set_arg('INT_D',DEF_OP_FROM_OCC,'LABEL',1,tgt_info,
     &             val_label=(/'INT_D'/))
      call set_arg('INT_D',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &             val_int=(/2/))
      call set_arg('INT_D',DEF_OP_FROM_OCC,'DESCR',1,tgt_info,
     &     val_str=descr)
      

      ! "redundant" parts of T amplitudes (in principle arbitrary)
      call add_target('T(2)red',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = max(2-ip,0), min(2-ip,maxh)
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ! purely inactive excitations have no redundant component
            if (max(iexc-ip,iexc-ih).eq.0) cycle
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
      call set_rule('T(2)red',ttype_op,DEF_OP_FROM_OCC,
     &              'T(2)red',1,1,
     &              parameters,2,tgt_info)
      call add_target('T(3)red',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = max(3-ip,0), min(3-ip,maxh)
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            ndef = ndef + 1
            occ_def(IHOLE,2,ndef) = ih
            occ_def(IPART,1,ndef) = ip
            occ_def(IVALE,1,ndef) = iexc - ip
            occ_def(IVALE,2,ndef) = iexc - ih
          end do
        end do
      end do
      ex_t3red = ndef.gt.0
      call op_from_occ_parameters(-1,parameters,2,
     &              occ_def,ndef,1,(/0,0/),ndef)
      call set_rule('T(3)red',ttype_op,DEF_OP_FROM_OCC,
     &              'T(3)red',1,1,
     &              parameters,2,tgt_info)

      ! define operator for fixed, read-in amplitudes
      call add_target('Tfix',ttype_op,.false.,tgt_info)
      occ_def = 0
      ndef = 0
      do ip = 0, maxp
        do ih = 0, maxh
          do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
            if (iexc.gt.tfix) cycle ! tfix is the maximum rank
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
      call set_rule('Tfix',ttype_op,DEF_OP_FROM_OCC,
     &              'Tfix',1,1,
     &              parameters,2,tgt_info)

      ! define fixed energy (e.g. E_icMRCCSD when doing a (T) calc.)
      call add_target('E_fix',ttype_op,.false.,tgt_info)
      call hop_parameters(-1,parameters,0,0,1,.false.)
      call set_rule('E_fix',ttype_op,DEF_HAMILTONIAN,'E_fix',
     &              1,1,parameters,1,tgt_info)
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
      if (.not.Op_eqs) then
        call set_rule2('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,tgt_info)
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'LABEL',1,
     &       tgt_info,val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &       tgt_info,val_label=(/'NORM'/))
        if (h1bar) then
          call set_dependency('F_MRCC_LAG','H1bar',tgt_info)
          call set_dependency('F_MRCC_LAG','T-T1',tgt_info)
          call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &         tgt_info,val_label=(/'L    ','H1bar','T-T1 ','C0   '/))
        else
          call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &         tgt_info,val_label=(/'L ','H ','T ','C0'/))
        end if
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &       tgt_info,val_int=(/maxcom/))
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &       tgt_info,val_int=(/maxcom_en/))
        call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &       val_str='---')
          call set_arg('F_MRCC_LAG',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &       tgt_info,val_str='ic-MRCC Lagrangian')
      else
        call set_dependency('F_MRCC_LAG','Heff',tgt_info)
        call set_dependency('F_MRCC_LAG','Geff',tgt_info)
        call set_rule2('F_MRCC_LAG',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'C0^+','Heff','C0  '/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'IDX_SV',3,
     &       tgt_info,
     &       val_int=(/2,3,4/))
        call set_rule2('F_MRCC_LAG',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'LABEL',1,
     &       tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OP_RES',1,
     &       tgt_info,
     &       val_label=(/'NORM'/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &       tgt_info,
     &       val_label=(/'C0^+','L   ','Geff','C0  '/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'IDX_SV',4,
     &       tgt_info,
     &       val_int=(/2,3,4,5/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end if
      if (optref.eq.-1.or.optref.eq.-2) then
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
     &       val_label=(/'E(MR)','C0^+ ','C0   '/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'IDX_SV',3,
     &       tgt_info,
     &       val_int=(/2,3,4/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &       val_rl8=(/-1d0/))
        call set_arg('F_MRCC_LAG',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/.false./))
      end if
      if (h1bar) then
        call set_rule2('F_MRCC_LAG',REPLACE,tgt_info)
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'T-T1','T   '/))
        if (htt) then
          call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &         tgt_info,val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &         tgt_info,val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',3,
     &         tgt_info,val_label=(/'H1bar','H    ','T    '/))
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &         val_str='HTT')
        end if
c dbg
c        call set_dependency('F_MRCC_LAG','F_H1bar',tgt_info)
c        call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
c        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_MRCC_LAG'/))
c        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_MRCC_LAG'/))
c        call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
c     &       val_label=(/'F_H1bar'/))
c dbgend
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
      if (.not.Op_eqs) then
        call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
c        if (project) then !remove terms that are projected out
c          if (h1bar) then
c            call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',4,
c     &         tgt_info,val_label=(/'H1bar','T','C0','L'/))
c          else
c            call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',4,
c     &         tgt_info,val_label=(/'H','T','C0','L'/))
c          end if
c        else
          if (h1bar) then
            call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,
     &         tgt_info,val_label=(/'H1bar','T    '/))
          else
            call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,
     &         tgt_info,val_label=(/'H','T'/))
          end if
c        end if
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC2')
        if (G_level.lt.0) ! approximation will change factors
     &       call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &                    val_str='CHECK_FAC')
      else
        ! Separate equation for each block of Geff
        call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &       val_label=(/'L   ','Geff'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='SAME')
      end if
      if (tfix.gt.0) then 
        ! remove redundant (zero) part of residuals
        call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',3,tgt_info,
     &       val_label=(/'L   ','T   ','Tfix'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCCrem0res')
      end if
      if (.not.Op_eqs.and.trunc.and.t1ord.ge.0) then
        ! Factor out fixed part of energy (not to be truncated)
        call set_dependency('F_MRCC_LAG','F_Efix',tgt_info)
        call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_Efix'/))
      end if
      if (.not.Op_eqs.and.h1bar) then
        if (trunc) then
          ! prescreening: remove terms with definitely too high order
          call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,
     &         tgt_info,val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,
     &         tgt_info,val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',4,
     &       tgt_info,val_label=(/'H1bar','H1bar','T    ','L    '/)) ! #2 is dummy
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &         val_str='COUNT_L/PRESCREEN')
          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &         val_str='MRCCtrunc')
        end if
        call set_dependency('F_MRCC_LAG','F_H1bar',tgt_info)
        if (h1bar_maxp.lt.4) then
          ! expand formal part of H1bar (store H1bar blks. only up to P^2)
          call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_H1barformal'/))
          call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
     &         val_int=(/2/))
        end if
        if (trunc) then
          ! expand H1bar fully. Factor out again after truncation
          call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_H1bar'/))
          call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
     &         val_int=(/2/))
        end if
      end if
      if (.not.Op_eqs.and.trunc) then
        ! apply perturbative truncation of Lagrangian
        call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
        call set_dependency('F_MRCC_LAG','FREF',tgt_info)
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
c        if (h1bar) then
c          call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',4,
c     &       tgt_info,val_label=(/'H1bar','FREF','T','L'/))
c        else
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',4,
     &     tgt_info,val_label=(/'H   ','FREF','T   ','L   '/))
c        end if
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='COUNT_L')
        call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCCtrunc')
        if (t1ord.ge.0) then
          ! expand fixed part of energy again (not necessary in principle)
          call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_Efix'/))
        end if
        if (h1bar) then
          ! Factor out H1bar again (except the formal part)
          call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
          call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &         val_label=(/'F_MRCC_LAG'/))
          call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &         val_label=(/'F_H1bar'/))
        end if
        ! expand effective Fock operator (if it was inserted)
        call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
c        if (h1bar) then
c          call set_dependency('F_MRCC_LAG','F_FREFbar',tgt_info)
c          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
c     &         val_label=(/'F_FREFbar'/))
c        else
        call set_dependency('F_MRCC_LAG','F_FREF',tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_FREF'/))
c        end if
      end if
c dbg
c      ! h) let only diagonal blocks of Jacobian survive
c      call set_rule2('F_MRCC_LAG',SELECT_SPECIAL,tgt_info)
c      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
c     &     val_label=(/'L','T'/))
c      call set_arg('F_MRCC_LAG',SELECT_SPECIAL,'TYPE',1,tgt_info,
c     &     val_str='SAME')
c dbgend
c dbg
      if (fact_tt) then
        ! factor out TT intermediate
        call set_dependency('F_MRCC_LAG','F_TT',tgt_info)
        call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
        call set_rule2('F_MRCC_LAG',FACTOR_OUT,tgt_info)
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
      end if
c dbgend
      if (tfix.gt.0) then
        ! replace fixed part of T
        call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
        call set_dependency('F_MRCC_LAG','F_Tfix',tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_Tfix'/))
        call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
     &       val_int=(/2/))
        ! do the same for Lambda (somewhat awkwardly in three steps)
        call set_rule2('F_MRCC_LAG',REPLACE,tgt_info)
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'L  ','T^+'/))
        call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_Tfix^+'/))
        call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
     &       val_int=(/2/))
        call set_rule2('F_MRCC_LAG',REPLACE,tgt_info)
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
        call set_arg('F_MRCC_LAG',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'T^+','L  '/))
        if (h1bar) then
          ! we also need to replace fixed part of T in H1bar!
          call set_rule2('F_MRCC_LAG',EXPAND,tgt_info)
          call set_dependency('F_MRCC_LAG','F_Tfix',tgt_info)
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_RES',1,tgt_info,
     &         val_label=(/'F_H1bar'/))
          call set_arg('F_MRCC_LAG',EXPAND,'LABEL_IN',1,tgt_info,
     &         val_label=(/'F_H1bar'/))
          call set_arg('F_MRCC_LAG',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_Tfix'/))
          call set_arg('F_MRCC_LAG',EXPAND,'IMODE',1,tgt_info,
     &         val_int=(/2/))
        end if
      end if
c dbg
c      call set_rule2('F_MRCC_LAG',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_LAG',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_LAG'/))
c dbgend

      ! Residual part of Lagrangian
      call add_target2('F_LAG_L',.false.,tgt_info)
      call set_dependency('F_LAG_L','F_MRCC_LAG',tgt_info)
      call set_dependency('F_LAG_L','NORM',tgt_info)
      call set_dependency('F_LAG_L','L',tgt_info)
      call set_rule2('F_LAG_L',SELECT_TERMS,tgt_info)
      call set_arg('F_LAG_L',SELECT_TERMS,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_LAG_L'/))
      call set_arg('F_LAG_L',SELECT_TERMS,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_LAG_L',SELECT_TERMS,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_LAG_L',SELECT_TERMS,'OP_INCL',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_LAG_L',SELECT_TERMS,'BLK_INCL',1,tgt_info,
     &     val_int=(/0/))
      ! Cumulant approximation?
      if (maxcum.gt.0) then
        ! Factor out reduced density matrices
        call set_rule2('F_LAG_L',FACTOR_OUT,tgt_info)
        call set_dependency('F_LAG_L','F_DENS0',tgt_info)
        call set_arg('F_LAG_L',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_LAG_L'/))
        call set_arg('F_LAG_L',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_LAG_L'/))
        call set_arg('F_LAG_L',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_DENS0'/))
        ! Expand density matrices in terms of cumulants
        call set_rule2('F_LAG_L',EXPAND,tgt_info)
        call set_arg('F_LAG_L',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_LAG_L'/))
        call set_arg('F_LAG_L',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_LAG_L'/))
        if (maxcom.le.2.and.cum_appr_mode.ge.2.and.
     &      maxcum.ge.2.and.maxcum.le.4) then
          call set_dependency('F_LAG_L','F_DENS_appr',tgt_info)
          call set_arg('F_LAG_L',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_DENS_appr'/))
        else
          call set_dependency('F_LAG_L','F_DENS',tgt_info)
          call set_arg('F_LAG_L',EXPAND,'INTERM',1,tgt_info,
     &         val_label=(/'F_DENS'/))
        end if
        call set_arg('F_LAG_L',EXPAND,'IMODE',1,tgt_info,
     &       val_int=(/2/)) ! keep low-rank RDMs
c        ! Delete disconnected terms (cum. only connected to L)
c        call set_rule2('F_LAG_L',SELECT_SPECIAL,tgt_info)
c        call set_arg('F_LAG_L',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        call set_arg('F_LAG_L',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c        call set_arg('F_LAG_L',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
c     &       val_label=(/'CUM','L'/))
c        call set_arg('F_LAG_L',SELECT_SPECIAL,'TYPE',1,tgt_info,
c     &       val_str='MRCC3')
      end if
c dbg
c        call set_rule2('F_LAG_L',PRINT_FORMULA,tgt_info)
c        call set_arg('F_LAG_L',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &       val_label=(/'F_LAG_L'/))
c dbgend

      ! Residual
      call add_target2('F_OMG',.false.,tgt_info)
      call set_dependency('F_OMG','F_MRCC_LAG',tgt_info)
      call set_dependency('F_OMG','OMG',tgt_info)
      call set_rule2('F_OMG',DERIVATIVE,tgt_info)
      call set_arg('F_OMG',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      if (maxcum.gt.0) then
        call set_dependency('F_OMG','F_LAG_L',tgt_info)
        call set_arg('F_OMG',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_LAG_L'/))
      else
        call set_arg('F_OMG',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
      end if
      call set_arg('F_OMG',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      call set_arg('F_OMG',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))
      if (tfix.eq.0) then ! tfix>0: contains both T and Tfix
        call set_rule2('F_OMG',SELECT_SPECIAL,tgt_info)
        call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_OMG'/))
        call set_arg('F_OMG',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_OMG'/))
        call set_arg('F_OMG',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &       val_label=(/'H','T'/))
        call set_arg('F_OMG',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC2')
        call set_arg('F_OMG',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='CHECK    X')
      else if (maxit.eq.1) then ! T=0
        call set_rule2('F_OMG',INVARIANT,tgt_info)
        call set_arg('F_OMG',INVARIANT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_OMG'/))
        call set_arg('F_OMG',INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_OMG'/))
        call set_arg('F_OMG',INVARIANT,'OP_RES',1,tgt_info,
     &       val_label=(/'OMG'/))
        call set_arg('F_OMG',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_OMG',INVARIANT,'TITLE',1,tgt_info,
     &       val_str='Higher-order residual for first iteration')
      end if
      if (tfix.gt.0) then
        call set_rule2('F_OMG',PRINT_FORMULA,tgt_info)
        call set_arg('F_OMG',PRINT_FORMULA,'LABEL',1,tgt_info,
     &       val_label=(/'F_OMG'/))
      end if

      ! Lagrangian without Lambda...
      call add_target2('F_E_C0',.false.,tgt_info)
      call set_dependency('F_E_C0','L',tgt_info)
      call set_rule2('F_E_C0',INVARIANT,tgt_info)
      call set_arg('F_E_C0',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E_C0'/))
      if (use_f12) then
        call set_dependency('F_E_C0','F_MRCC_F12_LAG',tgt_info)
        call set_arg('F_E_C0',INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_F12_LAG'/))
      else
        call set_dependency('F_E_C0','F_MRCC_LAG',tgt_info)
        call set_arg('F_E_C0',INVARIANT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_LAG'/))
      end if
      call set_arg('F_E_C0',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_E_C0',INVARIANT,'OPERATORS',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_E_C0',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC energy expression for C0 equations')
c dbg
c      ! insert unit operators to allow for double differentiation
c      call set_dependency('F_E_C0','1v',tgt_info)
c      call set_rule2('F_E_C0',INSERT,tgt_info)
c      call set_arg('F_E_C0',INSERT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_E_C0'/))
c      call set_arg('F_E_C0',INSERT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_E_C0'/))
c      call set_arg('F_E_C0',INSERT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_E_C0',INSERT,'OP_INS',1,tgt_info,
c     &     val_label=(/'1v'/))
c      call set_arg('F_E_C0',INSERT,'OP_INCL',2,tgt_info,
c     &     val_label=(/'C0^+','C0'/))
c dbgend

      ! Residual for C0
      call add_target2('F_OMG_C0',.false.,tgt_info)
      call set_dependency('F_OMG_C0','F_E_C0',tgt_info)
      call set_dependency('F_OMG_C0','A_C0',tgt_info)
      call set_rule2('F_OMG_C0',DERIVATIVE,tgt_info)
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A_C0'/))
      call set_arg('F_OMG_C0',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'C0^+'/))
c      call set_rule2('F_OMG_C0',INVARIANT,tgt_info)
c      call set_arg('F_OMG_C0',INVARIANT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_OMG_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'OP_RES',1,tgt_info,
c     &     val_label=(/'A_C0'/))
c      call set_arg('F_OMG_C0',INVARIANT,'OPERATORS',1,tgt_info,
c     &     val_label=(/'L'/))
c      call set_arg('F_OMG_C0',INVARIANT,'TITLE',1,tgt_info,
c     &     val_str='Residual for Reference function')
      call set_rule2('F_OMG_C0',SELECT_SPECIAL,tgt_info)
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG_C0'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'H','T'/))
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
      call set_arg('F_OMG_C0',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='CHECK    X')
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
      if (tfix.eq.0) then
        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/'L    ','E(MR)'/))
      else if (maxit.gt.1) then
        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',3,tgt_info,
     &       val_label=(/'L     ','Tfix^+','E(MR) '/))
      else
        call set_arg('F_MRCC_E',INVARIANT,'OPERATORS',4,tgt_info,
     &       val_label=(/'L     ','Tfix^+','E(MR) ','T     '/))
      end if
      call set_arg('F_MRCC_E',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='MRCC energy expression')
      if (tfix.ne.0) then
        call set_rule2('F_MRCC_E',SELECT_SPECIAL,tgt_info)
        call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_MRCC_E'/))
        call set_arg('F_MRCC_E',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_MRCC_E'/))
        call set_arg('F_MRCC_E',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &       val_label=(/'H','T'/))
        call set_arg('F_MRCC_E',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC2')
        call set_arg('F_MRCC_E',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='CHECK    X')
      end if
      call set_rule2('F_MRCC_E',PRINT_FORMULA,tgt_info)
      call set_arg('F_MRCC_E',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_MRCC_E'/))

      ! transformed excitation operator
      call add_target2('F_T',.false.,tgt_info)
      call set_dependency('F_T','T',tgt_info)
      call set_dependency('F_T','Dtr',tgt_info)
      call set_dependency('F_T','Ttr',tgt_info)
      do i_cls = 1, nsupD
        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'T  ','Dtr','Ttr','Dtr','T  '/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndD(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndD(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/4/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &       val_int=(/3,5,2,4,1,4,2,5/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      do i_cls = nsupD+1, nsupT
        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'T  ','Ttr','T  '/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',2,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',2,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
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
c      do i_cls = 1, nsupD
c        call set_rule2('F_T',EXPAND_OP_PRODUCT,tgt_info)
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &       val_label=(/'F_T'/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &       val_label=(/'T'/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'OPERATORS',3,
c     &       tgt_info,
c     &       val_label=(/'T','Ttr','T'/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
c     &       val_int=(/1,2,1/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
c     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls),
c     &                 stndT(1,i_cls)/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
c     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls),
c     &                 stndT(2,i_cls)/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
c     &       val_rl8=(/-1d0/))
c        call set_arg('F_T',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &       val_log=(/.false./))
c      end do
c dbgend
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
      do i_cls = 1, nsupD
        call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'L    ','Dtr^+','Ltr  ','Dtr^+','L    '/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndD(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndD(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/4/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &       val_int=(/1,3,2,4,1,4,2,5/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      do i_cls = nsupD+1, nsupT
        call set_rule2('F_L',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_L',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'L'/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'L  ','Ltr','L  '/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MIN',2,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'BLK_MAX',2,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_L',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
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

      ! generic version of the above
      call add_target2('F_Ttr_GEN',.false.,tgt_info)
      call set_dependency('F_Ttr_GEN','Tout',tgt_info)
      call set_dependency('F_Ttr_GEN','Dtr',tgt_info)
      call set_dependency('F_Ttr_GEN','Ttr',tgt_info)
      do i_cls = 1, nsupD
        call set_rule2('F_Ttr_GEN',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_Ttr_GEN'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'Tout'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'OPERATORS',5,
     &       tgt_info,
     &       val_label=(/'Tout','Dtr ','Ttr ','Dtr ','Tout'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'IDX_SV',5,tgt_info,
     &       val_int=(/1,2,3,2,1/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndD(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndD(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &       val_int=(/4/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'AVOID',8,tgt_info,
     &       val_int=(/3,5,2,4,1,4,2,5/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      do i_cls = nsupD+1, nsupT
        call set_rule2('F_Ttr_GEN',EXPAND_OP_PRODUCT,tgt_info)
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &       val_label=(/'F_Ttr_GEN'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &       val_label=(/'Tout'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &       tgt_info,
     &       val_label=(/'Tout','Ttr ','Tout'/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &       val_int=(/1,2,1/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'BLK_MIN',2,tgt_info,
     &       val_int=(/stndT(1,i_cls),stndT(1,i_cls)/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'BLK_MAX',2,tgt_info,
     &       val_int=(/stndT(2,i_cls),stndT(2,i_cls)/))
        call set_arg('F_Ttr_GEN',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &       val_log=(/i_cls.eq.1/))
      end do
      ! no active open lines from Ttr
      call set_rule2('F_Ttr_GEN',SELECT_LINE,tgt_info)
      call set_arg('F_Ttr_GEN',SELECT_LINE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Ttr_GEN'/))
      call set_arg('F_Ttr_GEN',SELECT_LINE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_Ttr_GEN'/))
      call set_arg('F_Ttr_GEN',SELECT_LINE,'OP_RES',1,tgt_info,
     &     val_label=(/'Tout'/))
      call set_arg('F_Ttr_GEN',SELECT_LINE,'OP_INCL',1,tgt_info,
     &     val_label=(/'Ttr'/))
      call set_arg('F_Ttr_GEN',SELECT_LINE,'IGAST',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('F_Ttr_GEN',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='no_ext')
c dbg
c      call set_rule2('F_Ttr_GEN',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Ttr_GEN',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Ttr_GEN'/))
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
     &     tgt_info,val_label=(/'L ','H ','T ','C0'/))
c     &     tgt_info,val_label=(/'L','FREF','T','C0'/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/1/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/0/))
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='NOSCAL')
      call set_arg('F_E(MRCC)tr',DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='Precursor for linearized Jacobian')
      if (G_level.ge.0) then ! discard disconnected terms
       call set_rule2('F_E(MRCC)tr',SELECT_SPECIAL,tgt_info)
       call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &      val_label=(/'F_E(MRCC)tr'/))
       call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &      val_label=(/'F_E(MRCC)tr'/))
       call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &      val_label=(/'H','T'/))
       call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &      val_str='MRCC2')
      end if
      ! f) insert 1 (particle/hole space) for later differentiation
      if (prc_type.lt.3) then
        call set_dependency('F_E(MRCC)tr','1ph',tgt_info)
        call set_rule2('F_E(MRCC)tr',INSERT,tgt_info)
        call set_arg('F_E(MRCC)tr',INSERT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MRCC)tr'/))
        call set_arg('F_E(MRCC)tr',INSERT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MRCC)tr'/))
        call set_arg('F_E(MRCC)tr',INSERT,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_E(MRCC)tr',INSERT,'OP_INS',1,tgt_info,
     &       val_label=(/'1ph'/))
        call set_arg('F_E(MRCC)tr',INSERT,'OP_INCL',2,tgt_info,
     &       val_label=(/'L','T'/))
        ! replace 1ph by 1
        call set_dependency('F_E(MRCC)tr','1',tgt_info)
        call set_rule2('F_E(MRCC)tr',REPLACE,tgt_info)
        call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_E(MRCC)tr'/))
        call set_arg('F_E(MRCC)tr',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_E(MRCC)tr'/))
        call set_arg('F_E(MRCC)tr',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'1ph','1  '/))
      end if
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
      ! h) let only diagonal blocks of Jacobian survive
      call set_rule2('F_E(MRCC)tr',SELECT_SPECIAL,tgt_info)
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'OPERATORS',2,tgt_info,
     &     val_label=(/'Ltr','Ttr'/))
      call set_arg('F_E(MRCC)tr',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='SAME')
c dbg
      call set_rule2('F_E(MRCC)tr',PRINT_FORMULA,tgt_info)
      call set_arg('F_E(MRCC)tr',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
c dbgend

      ! (transformed) Jacobian times vector
      call add_target2('F_A_Ttr',.false.,tgt_info)
      call set_dependency('F_A_Ttr','F_E(MRCC)tr',tgt_info)
      call set_dependency('F_A_Ttr','OMGtr',tgt_info)
      call set_rule2('F_A_Ttr',DERIVATIVE,tgt_info)
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_E(MRCC)tr'/))
      call set_arg('F_A_Ttr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'OMGtr'/))
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
      call set_dependency('F_Atr','A',tgt_info)
      call set_rule2('F_Atr',DERIVATIVE,tgt_info)
      call set_arg('F_Atr',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Atr'/))
      call set_arg('F_Atr',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_A_Ttr'/))
      call set_arg('F_Atr',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'A'/))
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
c dbg
c      call set_rule2('F_Atr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Atr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Atr'/))
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
      iblk = 0
      kblk = 0
      first = .true.
      do ip = 0, maxp
       do iexc = excrestr(0,ip,1), excrestr(0,ip,2)
        iblk = iblk + 1
        jblk = 0
        do ih = 0, maxh
         do jexc = excrestr(ih,0,1), excrestr(ih,0,2)
          jblk = jblk + 1
          if (iexc.ne.jexc) cycle
          kblk = kblk + 1
          call set_rule2('F_TT',EXPAND_OP_PRODUCT,tgt_info)
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &         val_label=(/'F_TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &         val_label=(/'TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'OPERATORS',3,
     &         tgt_info,
     &         val_label=(/'TT','T ','TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
     &         val_int=(/1,2,1/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'BLK_MIN',3,tgt_info,
     &         val_int=(/kblk,-1,kblk/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'BLK_MAX',3,tgt_info,
     &         val_int=(/kblk,-1,kblk/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &         val_log=(/first/))
          first = .false.
          call set_rule2('F_TT',EXPAND_OP_PRODUCT,tgt_info)
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &         val_label=(/'F_TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &         val_label=(/'TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &         tgt_info,
     &         val_label=(/'TT','T ','T ','TT'/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &         val_int=(/1,2,3,1/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'BLK_MIN',4,tgt_info,
     &         val_int=(/kblk,iblk,jblk,kblk/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
     &         val_int=(/kblk,iblk,jblk,kblk/))
          call set_arg('F_TT',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &         val_log=(/.false./))
         end do
        end do
       end do
       do ih = 1, maxh
        do iexc = excrestr(ih,ip,1), excrestr(ih,ip,2)
          iblk = iblk + 1
        end do
       end do
      end do
      call set_rule2('F_TT',PRINT_FORMULA,tgt_info)
      call set_arg('F_TT',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_TT'/))

c dbg
c      ! transformed multireference CC norm
c      ! a) set up
c      call add_target2('F_NORMtr',.false.,tgt_info)
c      call set_dependency('F_NORMtr','NORM',tgt_info)
c      call set_dependency('F_NORMtr','T',tgt_info)
c      call set_dependency('F_NORMtr','L',tgt_info)
c      call set_dependency('F_NORMtr','DENS',tgt_info)
c      call set_rule2('F_NORMtr',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OPERATORS',2,
c     &     tgt_info,
c     &     val_label=(/'L','T'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'IDX_SV',2,tgt_info,
c     &     val_int=(/2,3/))
c      call set_rule2('F_NORMtr',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'DENS','L','T','DENS'/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/2,3,4,2/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'AVOID',2,tgt_info,
c     &     val_int=(/1,4/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'BLK_MAX',4,tgt_info,
c     &     val_int=(/orb_info%nactel,-1,-1,orb_info%nactel/))
c      call set_arg('F_NORMtr',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &     val_log=(/.false./))
c      ! f) insert 1 (particle/hole space) for later differentiation
c      call set_dependency('F_NORMtr','1ph',tgt_info)
c      call set_rule2('F_NORMtr',INSERT,tgt_info)
c      call set_arg('F_NORMtr',INSERT,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',INSERT,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',INSERT,'OP_RES',1,tgt_info,
c     &     val_label=(/'NORM'/))
c      call set_arg('F_NORMtr',INSERT,'OP_INS',1,tgt_info,
c     &     val_label=(/'1ph'/))
c      call set_arg('F_NORMtr',INSERT,'OP_INCL',2,tgt_info,
c     &     val_label=(/'L','T'/))
c      ! replace 1ph by 1
c      call set_dependency('F_NORMtr','1',tgt_info)
c      call set_rule2('F_NORMtr',REPLACE,tgt_info)
c      call set_arg('F_NORMtr',REPLACE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',REPLACE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',REPLACE,'OP_LIST',2,tgt_info,
c     &     val_label=(/'1ph','1'/))
c      ! g) expand T and L
c      call set_dependency('F_NORMtr','F_T',tgt_info)
c      call set_dependency('F_NORMtr','F_L',tgt_info)
c      call set_rule2('F_NORMtr',EXPAND,tgt_info)
c      call set_arg('F_NORMtr',EXPAND,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_NORMtr',EXPAND,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
cc      call set_arg('F_NORMtr',EXPAND,'INTERM',1,tgt_info,
cc     &     val_label=(/'F_T'/))
c      call set_arg('F_NORMtr',EXPAND,'INTERM',2,tgt_info,
c     &     val_label=(/'F_T','F_L'/))
cc dbg
c      call set_rule2('F_NORMtr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_NORMtr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
cc dbgend
c
c      ! transformed Metric times vector
c      call add_target2('F_S_Ttr',.false.,tgt_info)
c      call set_dependency('F_S_Ttr','F_NORMtr',tgt_info)
cc      call set_dependency('F_S_Ttr','OMG',tgt_info) !tr',tgt_info)
c      call set_dependency('F_S_Ttr','OMGtr',tgt_info)
c      call set_rule2('F_S_Ttr',DERIVATIVE,tgt_info)
c      call set_arg('F_S_Ttr',DERIVATIVE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_NORMtr'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'OP_RES',1,tgt_info,
cc     &     val_label=(/'OMG'/)) !tr'/))
c     &     val_label=(/'OMGtr'/))
c      call set_arg('F_S_Ttr',DERIVATIVE,'OP_DERIV',1,tgt_info,
cc     &     val_label=(/'L'/)) !tr'/))
c     &     val_label=(/'Ltr'/))
c      call set_rule2('F_S_Ttr',PRINT_FORMULA,tgt_info)
c      call set_arg('F_S_Ttr',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c
c      ! transformed Metric (should be unity)
c      call add_target2('F_Str',.false.,tgt_info)
c      call set_dependency('F_Str','F_S_Ttr',tgt_info)
c      call set_dependency('F_Str','S',tgt_info)
c      call set_rule2('F_Str',DERIVATIVE,tgt_info)
c      call set_arg('F_Str',DERIVATIVE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Str'/))
c      call set_arg('F_Str',DERIVATIVE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_S_Ttr'/))
c      call set_arg('F_Str',DERIVATIVE,'OP_RES',1,tgt_info,
c     &     val_label=(/'S'/))
c      call set_arg('F_Str',DERIVATIVE,'OP_DERIV',1,tgt_info,
c     &     val_label=(/'Ttr'/))
c      call set_rule2('F_Str',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Str',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Str'/))
c dbgend

      ! Effective Hamiltonian in the active space
      call add_target2('F_Heff',.false.,tgt_info)
      call set_dependency('F_Heff','Heff',tgt_info)
      call set_dependency('F_Heff','H',tgt_info)
      call set_dependency('F_Heff','T',tgt_info)
      call set_rule2('F_Heff',DEF_MRCC_INTM,tgt_info)
      call set_arg('F_Heff',DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &     val_label=(/'F_Heff'/))
      call set_arg('F_Heff',DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &     val_label=(/'Heff'/))
      call set_arg('F_Heff',DEF_MRCC_INTM,'OPERATORS',2,tgt_info,
     &     val_label=(/'T','H'/))
      call set_arg('F_Heff',DEF_MRCC_INTM,'MAXCOM',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_Heff',DEF_MRCC_INTM,'MODE',1,tgt_info,
     &     val_str='Heff')
      call set_arg('F_Heff',DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &     val_str='Effective Hamiltonian in the active space')
c dbg
      call set_rule2('F_Heff',PRINT_FORMULA,tgt_info)
      call set_arg('F_Heff',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_Heff'/))
c dbgend

      ! Effective Hamiltonian in the excitation space
      call add_target2('F_Geff',.false.,tgt_info)
      call set_dependency('F_Geff','Geff',tgt_info)
      call set_dependency('F_Geff','Heff',tgt_info)
      call set_dependency('F_Geff','H',tgt_info)
      call set_dependency('F_Geff','T',tgt_info)
      call set_rule2('F_Geff',DEF_MRCC_INTM,tgt_info)
      call set_arg('F_Geff',DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &     val_label=(/'F_Geff'/))
      call set_arg('F_Geff',DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &     val_label=(/'Geff'/))
      call set_arg('F_Geff',DEF_MRCC_INTM,'OPERATORS',3,tgt_info,
     &     val_label=(/'T   ','H   ','Heff'/))
      call set_arg('F_Geff',DEF_MRCC_INTM,'MAXCOM',1,
     &     tgt_info,val_int=(/maxcom/))
      call set_arg('F_Geff',DEF_MRCC_INTM,'MODE',1,tgt_info,
     &     val_str='Geff')
      call set_arg('F_Geff',DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &     val_str='Effective Hamiltonian in the excitation space')
c dbg
      call set_rule2('F_Geff',PRINT_FORMULA,tgt_info)
      call set_arg('F_Geff',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_Geff'/))
c dbgend

      ! T1 transformed Hamiltonian
      call add_target2('F_H1barfull',.false.,tgt_info)
      call set_dependency('F_H1barfull','H1bar',tgt_info)
      call set_dependency('F_H1barfull','H',tgt_info)
      call set_dependency('F_H1barfull','T1',tgt_info)
      call set_rule2('F_H1barfull',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'LABEL',1,tgt_info,
     &     val_label=(/'F_H1bar'/))
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &     tgt_info,val_label=(/'H1bar'/))
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &     tgt_info,val_label=(/'T1   ','H    ','T1   ','H1bar'/))
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/0/)) !dummy only
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_h1bar/))
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='HBAR_ALL')
      call set_arg('F_H1barfull',DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='T1 transformed Hamiltonian')
c      if (maxtt.lt.0) then ! do not forbid T-T-contractions in H1bar!
        call set_rule2('F_H1barfull',SELECT_SPECIAL,tgt_info)
        call set_arg('F_H1barfull',SELECT_SPECIAL,'LABEL_RES',1,
     &       tgt_info,val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'OPERATORS',2,
     &       tgt_info,val_label=(/'H ','T1'/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCC2')
        call set_arg('F_H1barfull',SELECT_SPECIAL,'MODE',1,tgt_info,
     &       val_str='CHECK_FACX') ! 'X' means: do not forbid T-T-contr
c      end if
      call set_dependency('F_H1barfull','T',tgt_info)
      call set_rule2('F_H1barfull',REPLACE,tgt_info)
      call set_arg('F_H1barfull',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_H1bar'/))
      call set_arg('F_H1barfull',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_H1bar'/))
      call set_arg('F_H1barfull',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'T1','T '/))
      if (trunc.and.tfix.eq.0) then !tfix>0: no trunc. for fixed sol.
        ! apply perturbative truncation of Lagrangian
        call set_rule2('F_H1barfull',SELECT_SPECIAL,tgt_info)
        call set_dependency('F_H1barfull','FREF',tgt_info)
        call set_arg('F_H1barfull',SELECT_SPECIAL,'LABEL_RES',1,
     &       tgt_info,val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'OPERATORS',3,
     &       tgt_info,val_label=(/'H   ','FREF','T   '/))
        call set_arg('F_H1barfull',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &       val_str='MRCCtrunc')
        ! expand effective Fock operator (if it was inserted)
        call set_rule2('F_H1barfull',EXPAND,tgt_info)
        call set_arg('F_H1barfull',EXPAND,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',EXPAND,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_dependency('F_H1barfull','F_FREF',tgt_info)
        call set_arg('F_H1barfull',EXPAND,'INTERM',1,tgt_info,
     &       val_label=(/'F_FREF'/))
      end if
c dbg
      if (fact_tt) then
        ! factor out TT intermediate
        call set_dependency('F_H1barfull','F_TT',tgt_info)
        call set_rule2('F_H1barfull',FACTOR_OUT,tgt_info)
        call set_arg('F_H1barfull',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
        call set_rule2('F_H1barfull',FACTOR_OUT,tgt_info)
        call set_arg('F_H1barfull',FACTOR_OUT,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',FACTOR_OUT,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_H1bar'/))
        call set_arg('F_H1barfull',FACTOR_OUT,'INTERM',1,tgt_info,
     &       val_label=(/'F_TT'/))
      end if
c dbgend
c      call set_rule2('F_H1barfull',PRINT_FORMULA,tgt_info)
c      call set_arg('F_H1barfull',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_H1bar'/))

      ! T1 transformed Hamiltonian, split into formal part and rest
      call add_target2('F_H1bar',.false.,tgt_info)
      call set_dependency('F_H1bar','F_H1barfull',tgt_info)
      call set_rule2('F_H1bar',SELECT_SPECIAL,tgt_info)
      call set_arg('F_H1bar',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,val_label=(/'F_H1barformal'/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_H1bar'/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'OPERATORS',0,tgt_info,
     &     val_label=(/''/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='FORMAL')
      call set_arg('F_H1bar',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='extract')
      call set_rule2('F_H1bar',SELECT_SPECIAL,tgt_info)
      call set_arg('F_H1bar',SELECT_SPECIAL,'LABEL_RES',1,
     &     tgt_info,val_label=(/'F_H1bar'/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_H1bar'/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'OPERATORS',0,tgt_info,
     &     val_label=(/''/))
      call set_arg('F_H1bar',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='FORMAL')
      call set_arg('F_H1bar',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='delete')
c dbg
c      call set_rule2('F_H1bar',PRINT_FORMULA,tgt_info)
c      call set_arg('F_H1bar',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_H1barformal'/))
c      call set_rule2('F_H1bar',PRINT_FORMULA,tgt_info)
c      call set_arg('F_H1bar',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_H1bar'/))
c dbgend

      ! precursor for combining selected PP contractions into one step
      call add_target2('F_prePPint',.false.,tgt_info)
      call set_dependency('F_prePPint','F_OMG',tgt_info)
      call set_dependency('F_prePPint','H_PP',tgt_info)
      call set_rule2('F_prePPint',REPLACE,tgt_info)
      call set_arg('F_prePPint',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePPint'/))
      call set_arg('F_prePPint',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      if (h1bar) then
        call set_arg('F_prePPint',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H1bar','H_PP '/))
      else
        call set_arg('F_prePPint',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_PP'/))
      end if
      call set_rule2('F_prePPint',INVARIANT,tgt_info)
      call set_arg('F_prePPint',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePPint'/))
      call set_arg('F_prePPint',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePPint'/))
      call set_arg('F_prePPint',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      if (h1bar) then
        call set_arg('F_prePPint',INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/'H1bar','H    '/))
      else
        call set_arg('F_prePPint',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
      end if
      call set_arg('F_prePPint',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Precursor for INT_PP')
c dbg
c      call set_rule2('F_prePPint',PRINT_FORMULA,tgt_info)
c      call set_arg('F_prePPint',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_prePPint'/))
c dbgend

      call add_target2('F_PPint',.false.,tgt_info)
      call set_dependency('F_PPint','F_prePPint',tgt_info)
      call set_dependency('F_PPint','INT_PP',tgt_info)
      call set_dependency('F_PPint','F_OMG',tgt_info)
      call set_rule2('F_PPint',DERIVATIVE,tgt_info)
      call set_arg('F_PPint',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PPint'/))
      call set_arg('F_PPint',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePPint'/))
      call set_arg('F_PPint',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP'/))
      call set_arg('F_PPint',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'H_PP'/))
      call set_rule2('F_PPint',INVARIANT,tgt_info)
      call set_arg('F_PPint',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PPint'/))
      call set_arg('F_PPint',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_PPint'/))
      call set_arg('F_PPint',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP'/))
      call set_arg('F_PPint',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('F_PPint',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_PPint',PRINT_FORMULA,tgt_info)
c      call set_arg('F_PPint',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_PPint'/))
c dbgend

      call add_target2('F_prePP0int',.false.,tgt_info)
      call set_dependency('F_prePP0int','F_OMG',tgt_info)
      call set_dependency('F_prePP0int','H_PP0',tgt_info)
      call set_rule2('F_prePP0int',REPLACE,tgt_info)
      call set_arg('F_prePP0int',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePP0int'/))
      call set_arg('F_prePP0int',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      call set_arg('F_prePP0int',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'H   ','H_PP0'/))
      call set_rule2('F_prePP0int',INVARIANT,tgt_info)
      call set_arg('F_prePP0int',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_prePP0int'/))
      call set_arg('F_prePP0int',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePP0int'/))
      call set_arg('F_prePP0int',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      if (h1bar) then
        call set_arg('F_prePP0int',INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/'H1bar','H    '/))
      else
        call set_arg('F_prePP0int',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
      end if
      call set_arg('F_prePP0int',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Precursor for INT_PP0')
c dbg
c      call set_rule2('F_prePP0int',PRINT_FORMULA,tgt_info)
c      call set_arg('F_prePP0int',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_prePP0int'/))
c dbgend

      call add_target2('F_PP0int',.false.,tgt_info)
      call set_dependency('F_PP0int','F_prePP0int',tgt_info)
      call set_dependency('F_PP0int','INT_PP0',tgt_info)
      call set_dependency('F_PP0int','F_OMG',tgt_info)
      call set_rule2('F_PP0int',DERIVATIVE,tgt_info)
      call set_arg('F_PP0int',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PP0int'/))
      call set_arg('F_PP0int',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_prePP0int'/))
      call set_arg('F_PP0int',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP0'/))
      call set_arg('F_PP0int',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'H_PP0'/))
      call set_rule2('F_PP0int',INVARIANT,tgt_info)
      call set_arg('F_PP0int',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_PP0int'/))
      call set_arg('F_PP0int',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_PP0int'/))
      call set_arg('F_PP0int',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_PP0'/))
      call set_arg('F_PP0int',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('F_PP0int',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_PP0int',PRINT_FORMULA,tgt_info)
c      call set_arg('F_PP0int',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_PP0int'/))
c dbgend

      ! Formulae needed for preconditioning with non-diag. H0
      ! the Dyall Hamiltonian itself:
      !  - core Fock matrix
      !  - two-electron active contribution to effective Fock for HH and PP part
      !  - two-electron VVVV part
      call add_target3((/
     &     'target F_H0Dy(',
     &     '    depend OPS_FOR_H0INV',
     &     '    EXPAND_OP_PRODUCT(new=T,label=F_H0Dy,op_res=H0Dy,',
     &     '       operators=(H0Dy,H,H0Dy),',
     &     '       idx_sv=(1,2,1),',
     &     '       descr="2,,[HPV],[HPV]",avoid=(1,3)',
     &     '    )',
     &     '    EXPAND_OP_PRODUCT(new=F,label=F_H0Dy,op_res=H0Dy,',
     &     '       OPERATORS=(H0Dy,C00^+,H,C00,H0Dy),',
     &     '       IDX_SV=(1,2,3,4,1),',
     &     '       DESCR=("1,,H,","1,,,P","5,,P,","5,,,H",',
     &     '                              "3,,[HPV]V,[HPV]V"),',
     &     '       CONNECT=(2,3,3,4)',
     &     '    )',
     &     '    EXPAND_OP_PRODUCT(new=F,label=F_H0Dy,op_res=H0Dy,',
     &     '       OPERATORS=(H0Dy,H,H0Dy),',
     &     '       IDX_SV=(1,2,1),',
     &     '       DESCR="2,,VV,VV",avoid=(1,3)',
     &     '    )',
     &     '    PRINT_FORMULA(label=F_H0Dy)',
     &     ')'/),
     &     tgt_info)

      ! (a) the right-hand-side: projected OMG
      ! (b) the transformation: <C0|tau [H0,T] |C0>
      call add_target3((/
     &     'target F_FOR_H0INV( ',
     &     '    depend (OPS_FOR_H0INV,OMG,Dproj,C00)',
     &     '    EXPAND_OP_PRODUCT(new=T,label=F_OMGprj,op_res=OMGprj,',
     &     '       FAC=-1d0,',
     &     '       OPERATORS=(OMGprj,OMG,OMGprj,',
     &     '                OMGprj,OMG,OMGprj),',
     &     '       IDX_SV=(1,2,1,1,2,1),AVOID=(2,5)',
     &     '    )',
     &     '    EXPAND_OP_PRODUCT(new=T,label=F_D00,op_res=D00,',
     &     '       FAC=1d0,',
     &     '       OPERATORS=(D00,C00^+,D00,',
     &     '                D00,C00,D00),',
     &     '       IDX_SV=(1,2,1,1,3,1)',
     &     '    )',
c     &     '    EXPAND_OP_PRODUCT(new=T,label=F_OMGprj,op_res=OMGprj,',
c     &     '       FAC=-1d0,',
c     &     '       OPERATORS=(OMGprj,OMG,Dproj^+,OMGprj,',
c     &     '                OMGprj,Dproj^+,OMG,OMGprj),',
c     &     '       IDX_SV=(1,2,3,1,1,3,2,1)',
c     &     '    )',
c     &     '    SELECT_LINE(label_res=F_OMGprj,label_in=F_OMGprj,',
c     &     '        OP_RES=OMGprj,OP_INCL=OMG,IGAST=3,MODE=no_ext)',
     &     '    PRINT_FORMULA(label=F_OMGprj)',
     &     '    EXPAND_OP_PRODUCT(new=T,label=F_S_DlT,op_res=S_DlT,',
     &     '       OPERATORS=(S_DlT,C00^+,S_DlT,       ',
     &     '                  S_DlT,DlT,C00,S_DlT),',
     &     '       IDX_SV=(1,2,1,1,3,4,1)',
     &     '    )',
     &     '    PRINT_FORMULA(label=F_S_DlT)',
     &     '    EXPAND_OP_PRODUCT(new=T,label=F_H0_DlT,op_res=H0_DlT,',
     &     '       OPERATORS=(H0_DlT,C00^+,H0_DlT,       ',
     &     '                  H0_DlT,H0Dy,DlT,C00,H0_DlT),',
     &     '       IDX_SV=(1,2,1,1,3,4,5,1),',
     &     '       CONNECT=(5,6))',
     &     '    EXPAND_OP_PRODUCT(new=F,label=F_H0_DlT,op_res=H0_DlT,',
     &     '       FAC=-1d0,',
     &     '       OPERATORS=(H0_DlT,C00^+,H0_DlT,       ',
     &     '                  H0_DlT,DlT,H0Dy,C00,H0_DlT),',
     &     '       IDX_SV=(1,2,1,1,3,4,5,1),',
     &     '       CONNECT=(5,6))',
     &     '    PRINT_FORMULA(label=F_H0_DlT)',
     &     ')'/),
     &     tgt_info)

      ! precursor for combining selected HH contractions into one step
      call add_target2('F_preHHint',.false.,tgt_info)
      call set_dependency('F_preHHint','F_OMG',tgt_info)
      call set_dependency('F_preHHint','H_HH',tgt_info)
      call set_rule2('F_preHHint',REPLACE,tgt_info)
      call set_arg('F_preHHint',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_preHHint'/))
      call set_arg('F_preHHint',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_OMG'/))
      if (h1bar) then
        call set_arg('F_preHHint',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H1bar','H_HH '/))
      else
        call set_arg('F_preHHint',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'H   ','H_HH'/))
      end if
      call set_rule2('F_preHHint',INVARIANT,tgt_info)
      call set_arg('F_preHHint',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_preHHint'/))
      call set_arg('F_preHHint',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preHHint'/))
      call set_arg('F_preHHint',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'OMG'/))
      if (h1bar) then
        call set_arg('F_preHHint',INVARIANT,'OPERATORS',2,tgt_info,
     &       val_label=(/'H    ','H1bar'/))
      else
        call set_arg('F_preHHint',INVARIANT,'OPERATORS',1,tgt_info,
     &       val_label=(/'H'/))
      end if
      call set_arg('F_preHHint',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Precursor for INT_HH')
c dbg
c      call set_rule2('F_preHHint',PRINT_FORMULA,tgt_info)
c      call set_arg('F_preHHint',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_preHHint'/))
c dbgend

      call add_target2('F_HHint',.false.,tgt_info)
      call set_dependency('F_HHint','F_preHHint',tgt_info)
      call set_dependency('F_HHint','INT_HH',tgt_info)
      call set_dependency('F_HHint','F_OMG',tgt_info)
      call set_rule2('F_HHint',DERIVATIVE,tgt_info)
      call set_arg('F_HHint',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_HHint'/))
      call set_arg('F_HHint',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_preHHint'/))
      call set_arg('F_HHint',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_HH'/))
      call set_arg('F_HHint',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'H_HH'/))
      call set_rule2('F_HHint',INVARIANT,tgt_info)
      call set_arg('F_HHint',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_HHint'/))
      call set_arg('F_HHint',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_HHint'/))
      call set_arg('F_HHint',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'INT_HH'/))
      call set_arg('F_HHint',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('F_HHint',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.true./))
c dbg
c      call set_rule2('F_HHint',PRINT_FORMULA,tgt_info)
c      call set_arg('F_HHint',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_HHint'/))
c dbgend

      ! linear H(1+T2) transform, as suggested intermediate
      call add_target2('F_INT_HT2',.false.,tgt_info)
      call set_dependency('F_INT_HT2','H',tgt_info)
      call set_dependency('F_INT_HT2','T',tgt_info)
      call set_dependency('F_INT_HT2','INT_HT2',tgt_info)
      call set_rule2('F_INT_HT2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_HT2'/))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &             val_log=(/.true./))
C      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
C     &             val_label=(/'INT_HT2'/))
C      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OPERATORS',3,tgt_info,
C     &             val_label=(/'INT_HT2','H','INT_HT2'/))
C      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'IDX_SV',3,tgt_info,
C     &             val_int=(/1,2,1/))
C
C      call set_rule2('F_INT_HT2',EXPAND_OP_PRODUCT,tgt_info)
C      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
C     &             val_label=(/'F_INT_HT2'/))
C      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
C     &             val_log=(/.false./))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &             val_label=(/'INT_HT2'/))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OPERATORS',4,tgt_info,
     &             val_label=(/'INT_HT2','H','T','INT_HT2'/))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &             val_int=(/1,2,3,1/))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'DESCR',1,tgt_info,
     &             val_label=(/'3,,[VP][VP],[HV][HV]'/))
      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'CONNECT',2,tgt_info,
     &             val_int=(/2,3/))

c      call set_rule2('F_INT_HT2',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &             val_label=(/'F_INT_HT2'/))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
c     &             val_log=(/.false./))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &             val_label=(/'INT_HT2'/))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'OPERATORS',4,tgt_info,
c     &             val_label=(/'INT_HT2','T','H','INT_HT2'/))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &             val_int=(/1,2,3,1/))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'DESCR',1,tgt_info,
c     &             val_label=(/'2,,[VP][VP],[HV][HV]'/))
c      call set_arg('F_INT_HT2',EXPAND_OP_PRODUCT,'CONNECT',2,tgt_info,
c     &             val_int=(/2,3/))
      call set_rule2('F_INT_HT2',PRINT_FORMULA,tgt_info)
      call set_arg('F_INT_HT2',PRINT_FORMULA,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_HT2'/))

      call add_target2('F_INT_T2H',.false.,tgt_info)
      call set_dependency('F_INT_T2H','H',tgt_info)
      call set_dependency('F_INT_T2H','T',tgt_info)
      call set_dependency('F_INT_T2H','INT_T2H',tgt_info)
      call set_rule2('F_INT_T2H',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_T2H'/))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &             val_log=(/.true./))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &             val_label=(/'INT_T2H'/))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'OPERATORS',4,tgt_info,
     &             val_label=(/'INT_T2H','T','H','INT_T2H'/))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &             val_int=(/1,2,3,1/))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'DESCR',1,tgt_info,
     &             val_label=(/'2,,[VP][VP],[HV][HV]'/))
      call set_arg('F_INT_T2H',EXPAND_OP_PRODUCT,'CONNECT',2,tgt_info,
     &             val_int=(/2,3/))
      call set_rule2('F_INT_T2H',PRINT_FORMULA,tgt_info)
      call set_arg('F_INT_T2H',PRINT_FORMULA,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_T2H'/))

      ! Norm and 1-density
      call add_target2('F_INT_D',.false.,tgt_info)
      call set_dependency('F_INT_D','C0',tgt_info)
      call set_dependency('F_INT_D','INT_D',tgt_info)
      call set_rule2('F_INT_D',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_INT_D',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_D'/))
      call set_arg('F_INT_D',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &             val_log=(/.true./))
      call set_arg('F_INT_D',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &             val_label=(/'INT_D'/))
      call set_arg('F_INT_D',EXPAND_OP_PRODUCT,'OPERATORS',6,tgt_info,
     &             val_label=(/'INT_D','C0^+','INT_D',
     &                                        'INT_D','C0','INT_D'/))
      call set_arg('F_INT_D',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &             val_int=(/1,2,1,1,3,1/))
      call set_rule2('F_INT_D',PRINT_FORMULA,tgt_info)
      call set_arg('F_INT_D',PRINT_FORMULA,'LABEL',1,tgt_info,
     &             val_label=(/'F_INT_D'/))

      ! "Redundant" part of T operator
      call add_target2('F_T(2)red',.false.,tgt_info)
      call set_dependency('F_T(2)red','T(2)red',tgt_info)
      call set_dependency('F_T(2)red','Dtr',tgt_info)
      call set_dependency('F_T(2)red','T',tgt_info)
      call set_rule2('F_T(2)red',DEF_MRCC_INTM,tgt_info)
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &     val_label=(/'F_T(2)red'/))
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &     val_label=(/'T(2)red'/))
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'OPERATORS',2,tgt_info,
     &     val_label=(/'T  ','Dtr'/))
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'MAXCOM',1,
     &     tgt_info,val_int=(/2/))
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'MODE',1,tgt_info,
     &     val_str='Tred')
      call set_arg('F_T(2)red',DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &     val_str='Redundant part of T(2)')
c dbg
      call set_rule2('F_T(2)red',PRINT_FORMULA,tgt_info)
      call set_arg('F_T(2)red',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_T(2)red'/))
c dbgend
      call add_target2('F_T(3)red',.false.,tgt_info)
      call set_dependency('F_T(3)red','T(3)red',tgt_info)
      call set_dependency('F_T(3)red','Dtr',tgt_info)
      call set_dependency('F_T(3)red','T',tgt_info)
      call set_rule2('F_T(3)red',DEF_MRCC_INTM,tgt_info)
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'LABEL',1,tgt_info,
     &     val_label=(/'F_T(3)red'/))
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'INTERM',1,tgt_info,
     &     val_label=(/'T(3)red'/))
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'OPERATORS',3,tgt_info,
     &     val_label=(/'T  ','Dtr'/))
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'MAXCOM',1,
     &     tgt_info,val_int=(/3/))
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'MODE',1,tgt_info,
     &     val_str='Tred')
      call set_arg('F_T(3)red',DEF_MRCC_INTM,'TITLE',1,tgt_info,
     &     val_str='Redundant part of T(3)')
      ! danger: double counting due to identical blocks in Dtr
      ! remove the duplicates:
      do i_cls = 1, nremblk
        call set_rule2('F_T(3)red',SELECT_TERMS,tgt_info)
        call set_arg('F_T(3)red',SELECT_TERMS,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_T(3)red'/))
        call set_arg('F_T(3)red',SELECT_TERMS,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_T(3)red'/))
        call set_arg('F_T(3)red',SELECT_TERMS,'OP_RES',1,tgt_info,
     &       val_label=(/'T(3)red'/))
        call set_arg('F_T(3)red',SELECT_TERMS,'OP_EXCL',1,tgt_info,
     &       val_label=(/'Dtr'/))
        call set_arg('F_T(3)red',SELECT_TERMS,'BLK_EXCL',1,tgt_info,
     &       val_int=(/remblk(i_cls)/))
      end do
c dbg
      call set_rule2('F_T(3)red',PRINT_FORMULA,tgt_info)
      call set_arg('F_T(3)red',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_T(3)red'/))
c dbgend

      ! Spin "expectation value"
      call add_target2('F_MRCC_S(S+1)',.false.,tgt_info)
      call set_dependency('F_MRCC_S(S+1)','F_MRCC_LAG',tgt_info)
      call set_dependency('F_MRCC_S(S+1)','F_REPL_H_S2',tgt_info)
      call set_dependency('F_MRCC_S(S+1)','S(S+1)',tgt_info)
      call set_rule2('F_MRCC_S(S+1)',DERIVATIVE,tgt_info)
      call set_arg('F_MRCC_S(S+1)',DERIVATIVE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',DERIVATIVE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_MRCC_S(S+1)',DERIVATIVE,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',DERIVATIVE,'OP_DERIV',1,tgt_info,
     &     val_label=(/'L'/))
      call set_arg('F_MRCC_S(S+1)',DERIVATIVE,'OP_MULT',1,tgt_info,
     &     val_label=(/'L'/))
      call set_rule2('F_MRCC_S(S+1)',EXPAND,tgt_info)
      call set_arg('F_MRCC_S(S+1)',EXPAND,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',EXPAND,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',EXPAND,'INTERM',1,tgt_info,
     &     val_label=(/'F_REPL_H_S2'/))
      call set_arg('F_MRCC_S(S+1)',EXPAND,'IMODE',1,tgt_info,
     &     val_int=(/1/))
      call set_rule2('F_MRCC_S(S+1)',REPLACE,tgt_info)
      call set_arg('F_MRCC_S(S+1)',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_MRCC_S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_S(S+1)'/))
      call set_arg('F_MRCC_S(S+1)',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'L  ','T^+'/))
c      call set_rule2('F_MRCC_S(S+1)',SELECT_TERMS,tgt_info)
c      call set_arg('F_MRCC_S(S+1)',SELECT_TERMS,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_MRCC_S(S+1)'/))
c      call set_arg('F_MRCC_S(S+1)',SELECT_TERMS,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_MRCC_S(S+1)'/))
c      call set_arg('F_MRCC_S(S+1)',SELECT_TERMS,'OP_RES',1,tgt_info,
c     &     val_label=(/'S(S+1)'/))
c      call set_arg('F_MRCC_S(S+1)',SELECT_TERMS,'OP_INCL',1,tgt_info,
c     &     val_label=(/'T'/))
c dbg
c      call set_rule2('F_MRCC_S(S+1)',PRINT_FORMULA,tgt_info)
c      call set_arg('F_MRCC_S(S+1)',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_MRCC_S(S+1)'/))
c dbgend

c      ! formula for effective one-el. operator containing H1bar
c      call add_target2('F_FREFbar',.false.,tgt_info)
c      call set_dependency('F_FREFbar','F_FREF',tgt_info)
c      call set_dependency('F_FREFbar','H1bar',tgt_info)
c      call set_rule2('F_FREFbar',REPLACE,tgt_info)
c      call set_arg('F_FREFbar',REPLACE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_FREFbar'/))
c      call set_arg('F_FREFbar',REPLACE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_FREF'/))
c      call set_arg('F_FREFbar',REPLACE,'OP_LIST',2,tgt_info,
c     &     val_label=(/'H','H1bar'/))

      ! fixed part of T operator
      call add_target2('F_Tfix',.false.,tgt_info)
      call set_dependency('F_Tfix','T',tgt_info)
      call set_dependency('F_Tfix','Tfix',tgt_info)
      call set_rule2('F_Tfix',DEF_FORMULA,tgt_info)
      call set_arg('F_Tfix',DEF_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_Tfix'/))
      call set_arg('F_Tfix',DEF_FORMULA,'FORMULA',1,tgt_info,
     &     val_str='T=Tfix')
c dbg
c      call set_rule2('F_Tfix',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Tfix',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Tfix'/))
c dbgend

      ! Lagrangian with L replaced by T^+ (= Energy + correction)
      call add_target2('F_Ecorrected',.false.,tgt_info)
      call set_dependency('F_Ecorrected','F_MRCC_LAG',tgt_info)
      call set_dependency('F_Ecorrected','E(MR)',tgt_info)
      call set_rule2('F_Ecorrected',INVARIANT,tgt_info)
      call set_arg('F_Ecorrected',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Ecorrected'/))
      call set_arg('F_Ecorrected',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_MRCC_LAG'/))
      call set_arg('F_Ecorrected',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E(MR)'/))
      call set_arg('F_Ecorrected',INVARIANT,'OPERATORS',2,tgt_info,
     &     val_label=(/'E(MR)','L    '/))
      call set_arg('F_Ecorrected',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Energy + correction from MRCC Lagrangian')
      ! preliminary: only (T) correction for F12
      ! as our current Lagrangian neglect F12 contributions
      if (orb_info%norb_hpv(IEXTR,1).gt.0) then
        call set_rule2('F_Ecorrected',SELECT_TERMS,tgt_info)
        call set_arg('F_Ecorrected',SELECT_TERMS,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_Ecorrected'/))
        call set_arg('F_Ecorrected',SELECT_TERMS,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_Ecorrected'/))
        call set_arg('F_Ecorrected',SELECT_TERMS,'OP_RES',1,tgt_info,
     &       val_label=(/'E(MR)'/))
        call set_arg('F_Ecorrected',SELECT_TERMS,'OP_INCL',1,tgt_info,
     &       val_label=(/'T'/))
        call set_arg('F_Ecorrected',SELECT_TERMS,'BLK_INCL',1,tgt_info,
     &       val_int=(/0/))
      end if
c      call set_rule2('F_Ecorrected',REPLACE,tgt_info)
c      call set_arg('F_Ecorrected',REPLACE,'LABEL_RES',1,tgt_info,
c     &     val_label=(/'F_Ecorrected'/))
c      call set_arg('F_Ecorrected',REPLACE,'LABEL_IN',1,tgt_info,
c     &     val_label=(/'F_Ecorrected'/))
c      call set_arg('F_Ecorrected',REPLACE,'OP_LIST',2,tgt_info,
c     &     val_label=(/'L  ','T^+'/))
c dbg
      call set_rule2('F_Ecorrected',PRINT_FORMULA,tgt_info)
      call set_arg('F_Ecorrected',PRINT_FORMULA,'LABEL',1,tgt_info,
     &     val_label=(/'F_Ecorrected'/))
c dbgend

      ! Just the fixed part of the energy
      call add_target2('F_Efix',.false.,tgt_info)
      call set_dependency('F_Efix','E_fix',tgt_info)
      call set_dependency('F_Efix','H',tgt_info)
      call set_dependency('F_Efix','C0',tgt_info)
      call set_dependency('F_Efix','T',tgt_info)
      call set_dependency('F_Efix','Tfix',tgt_info)
      call set_dependency('F_Efix','L',tgt_info)
      call set_rule2('F_Efix',DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'LABEL',1,
     &     tgt_info,val_label=(/'F_Efix'/))
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'OP_RES',1,
     &     tgt_info,val_label=(/'E_fix'/))
      if (h1bar) then
        call set_dependency('F_Efix','H1bar',tgt_info)
        call set_dependency('F_Efix','T-T1',tgt_info)
        call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &       tgt_info,val_label=(/'L    ','H1bar','T-T1 ','C0   '/))
      else
        call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'OPERATORS',4,
     &       tgt_info,val_label=(/'L ','H ','T ','C0'/))
      end if
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'MAXCOM_RES',1,
     &     tgt_info,val_int=(/0/)) ! just the energy
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',1,
     &     tgt_info,val_int=(/maxcom_en/))
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'MODE',1,tgt_info,
     &     val_str='---')
      call set_arg('F_Efix',DEF_MRCC_LAGRANGIAN,'TITLE',1,
     &     tgt_info,val_str='Energy equation')
      if (h1bar) then
        call set_rule2('F_Efix',REPLACE,tgt_info)
        call set_arg('F_Efix',REPLACE,'LABEL_RES',1,tgt_info,
     &       val_label=(/'F_Efix'/))
        call set_arg('F_Efix',REPLACE,'LABEL_IN',1,tgt_info,
     &       val_label=(/'F_Efix'/))
        call set_arg('F_Efix',REPLACE,'OP_LIST',2,tgt_info,
     &       val_label=(/'T-T1','T   '/))
      end if
      ! only the fixed part
      call set_rule2('F_Efix',REPLACE,tgt_info)
      call set_arg('F_Efix',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'T   ','Tfix'/))
      call set_rule2('F_Efix',INVARIANT,tgt_info)
      call set_arg('F_Efix',INVARIANT,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',INVARIANT,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',INVARIANT,'OP_RES',1,tgt_info,
     &     val_label=(/'E_fix'/))
      call set_arg('F_Efix',INVARIANT,'OPERATORS',2,tgt_info,
     &     val_label=(/'T','L'/))
      call set_arg('F_Efix',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='Fixed part of energy')
      call set_rule2('F_Efix',REPLACE,tgt_info)
      call set_arg('F_Efix',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'Tfix','T   '/))
      call set_rule2('F_Efix',SELECT_SPECIAL,tgt_info)
      call set_arg('F_Efix',SELECT_SPECIAL,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      call set_arg('F_Efix',SELECT_SPECIAL,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_Efix'/))
      if (h1bar) then
        call set_arg('F_Efix',SELECT_SPECIAL,'OPERATORS',2,
     &     tgt_info,val_label=(/'H1bar','T    '/))
      else
        call set_arg('F_Efix',SELECT_SPECIAL,'OPERATORS',2,
     &     tgt_info,val_label=(/'H','T'/))
      end if
      call set_arg('F_Efix',SELECT_SPECIAL,'TYPE',1,tgt_info,
     &     val_str='MRCC2')
c dbg
c      call set_rule2('F_Efix',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Efix',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Efix'/))
c dbgend

      ! spin expectation value <C0| T^+ S^2 T |C0>
      call add_target2('F_T_S2',.false.,tgt_info)
      call set_dependency('F_T_S2','S(S+1)',tgt_info)
      call set_dependency('F_T_S2','C0',tgt_info)
      call set_dependency('F_T_S2','S+',tgt_info)
      call set_dependency('F_T_S2','S-',tgt_info)
      call set_dependency('F_T_S2','Sz',tgt_info)
      call set_dependency('F_T_S2','Sz_dum',tgt_info)
      call set_dependency('F_T_S2','T',tgt_info)
      ! (a) 1/2*(S+S- + S-S+)
      call set_rule2('F_T_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_T_S2'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &     tgt_info,
     &     val_label=(/'C0^+','T^+ ','S+  ','S-  ','T   ','C0  '/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &     val_int=(/2,3,4,5,6,7/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_rule2('F_T_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_T_S2'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &     tgt_info,
     &     val_label=(/'C0^+','T^+ ','S-  ','S+  ','T   ','C0  '/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &     val_int=(/2,3,4,5,6,7/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/0.5d0/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      ! (b) + Sz^2 (Sz_dum is used to circumvent automatic "BCH" factor)
      call set_rule2('F_T_S2',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_T_S2'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'S(S+1)'/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'OPERATORS',6,
     &     tgt_info,
     &     val_label=(/'C0^+  ','T^+   ','Sz    ','Sz_dum',
     &                 'T     ','C0    '/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'IDX_SV',6,tgt_info,
     &     val_int=(/2,3,4,5,6,7/))
      call set_arg('F_T_S2',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.false./))
      call set_rule2('F_T_S2',REPLACE,tgt_info)
      call set_arg('F_T_S2',REPLACE,'LABEL_RES',1,tgt_info,
     &     val_label=(/'F_T_S2'/))
      call set_arg('F_T_S2',REPLACE,'LABEL_IN',1,tgt_info,
     &     val_label=(/'F_T_S2'/))
      call set_arg('F_T_S2',REPLACE,'OP_LIST',2,tgt_info,
     &     val_label=(/'Sz_dum','Sz    '/))
c dbg
c      call set_rule2('F_T_S2',PRINT_FORMULA,tgt_info)
c      call set_arg('F_T_S2',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_T_S2'/))
c dbgend

      ! norm <C0| T^+ T |C0>
      call add_target2('F_T_NORM',.false.,tgt_info)
      call set_dependency('F_T_NORM','NORM',tgt_info)
      call set_dependency('F_T_NORM','C0',tgt_info)
      call set_dependency('F_T_NORM','T',tgt_info)
      call set_rule2('F_T_NORM',EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('F_T_NORM',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
     &     val_label=(/'F_T_NORM'/))
      call set_arg('F_T_NORM',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
     &     val_label=(/'NORM'/))
      call set_arg('F_T_NORM',EXPAND_OP_PRODUCT,'OPERATORS',4,
     &     tgt_info,
     &     val_label=(/'C0^+','T^+ ','T   ','C0  '/))
      call set_arg('F_T_NORM',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
     &     val_int=(/2,3,4,5/))
c dbg
c      call set_rule2('F_T_NORM',PRINT_FORMULA,tgt_info)
c      call set_arg('F_T_NORM',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_T_NORM'/))
c dbgend

      ! Formula to add transformation to GNO to Trafo matrix
      call add_target2('F_Dinv_GNO',.false.,tgt_info)
      call set_dependency('F_Dinv_GNO','D',tgt_info)
      if (gno.gt.0) then
        call set_dependency('F_Dinv_GNO','CUM',tgt_info)
        do i_cls = 1, nsupD
          call set_rule2('F_Dinv_GNO',EXPAND_OP_PRODUCT,tgt_info)
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'LABEL',1,
     &         tgt_info,val_label=(/'F_Dinv_GNO'/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'OP_RES',1,
     &         tgt_info,val_label=(/'D'/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'OPERATORS',11,
     &         tgt_info,
     &         val_label=(/'D  ','D  ','CUM','D  ','D  ','CUM',
     &                     'D  ','D  ','D  ','D  ','D  '/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'IDX_SV',11,
     &         tgt_info,val_int=(/1,2,3,1,1,3,2,1,1,2,1/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'BLK_MIN',3,
     &         tgt_info,
     &         val_int=(/stndD(1,i_cls),stndD(1,i_cls),1/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'BLK_MAX',3,
     &         tgt_info,
     &         val_int=(/stndD(2,i_cls),stndD(2,i_cls),-1/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &         tgt_info,val_int=(/11/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'AVOID',22,
     &         tgt_info,
     &         val_int=(/2,7,2,8,2,10,5,8,5,10,7,10,
     &                   3,4,5,6,3,6,3,8,3,10/))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'FIX_VTX',1,
     &         tgt_info,val_log=(/.true./))
          call set_arg('F_Dinv_GNO',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &         val_log=(/i_cls.eq.1/))
        end do
      end if
c dbg
c      call set_rule2('F_Dinv_GNO',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Dinv_GNO',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Dinv_GNO'/))
c dbgend

      ! Formula to add transformation to GNO to Projector
      call add_target2('F_Dproj_GNO',.false.,tgt_info)
      call set_dependency('F_Dproj_GNO','D',tgt_info)
      if (gno.gt.0) then
        call set_dependency('F_Dproj_GNO','CUM',tgt_info)
        do i_cls = 1, nsupD
          ! first multiply from left
          call set_rule2('F_Dproj_GNO',EXPAND_OP_PRODUCT,tgt_info)
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'LABEL',1,
     &         tgt_info,val_label=(/'F_Dproj_GNO'/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'OP_RES',1,
     &         tgt_info,val_label=(/'D'/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'OPERATORS',11,
     &         tgt_info,
     &         val_label=(/'D  ','D  ','CUM','D  ','D  ','CUM',
     &                     'D  ','D  ','D  ','D  ','D  '/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'IDX_SV',11,
     &         tgt_info,val_int=(/1,2,3,1,1,3,2,1,1,2,1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'BLK_MIN',3,
     &         tgt_info,
     &         val_int=(/stndD(1,i_cls),stndD(1,i_cls),1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'BLK_MAX',3,
     &         tgt_info,
     &         val_int=(/stndD(2,i_cls),stndD(2,i_cls),-1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &         tgt_info,val_int=(/11/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'AVOID',22,
     &         tgt_info,
     &         val_int=(/2,7,2,8,2,10,5,8,5,10,7,10,
     &                   3,4,5,6,3,6,3,8,3,10/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'FIX_VTX',1,
     &         tgt_info,val_log=(/.true./))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &         val_log=(/i_cls.eq.1/))
          ! now from right: will also lead to contributions from both
          call set_rule2('F_Dproj_GNO',EXPAND_OP_PRODUCT,tgt_info)
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'LABEL',1,
     &         tgt_info,val_label=(/'F_Dproj_GNO'/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'OP_RES',1,
     &         tgt_info,val_label=(/'D'/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'OPERATORS',11,
     &         tgt_info,
     &         val_label=(/'D  ','D  ','D  ','D  ','D  ','CUM',
     &                     'D  ','D  ','CUM','D  ','D  '/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'IDX_SV',11,
     &         tgt_info,val_int=(/1,2,1,1,2,3,1,1,3,2,1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'BLK_MIN',3,
     &         tgt_info,
     &         val_int=(/stndD(1,i_cls),stndD(1,i_cls),1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'BLK_MAX',3,
     &         tgt_info,
     &         val_int=(/stndD(2,i_cls),stndD(2,i_cls),-1/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'N_AVOID',1,
     &         tgt_info,val_int=(/11/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'AVOID',22,
     &         tgt_info,
     &         val_int=(/2,5,2,7,2,10,4,7,4,10,5,10,
     &                   6,10,5,9,6,9,4,9,2,9/))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'FIX_VTX',1,
     &         tgt_info,val_log=(/.true./))
          call set_arg('F_Dproj_GNO',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &         val_log=(/.false./))
        end do
      end if
c dbg
c      call set_rule2('F_Dproj_GNO',PRINT_FORMULA,tgt_info)
c      call set_arg('F_Dproj_GNO',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_Dproj_GNO'/))
c dbgend
*----------------------------------------------------------------------*
*     Opt. Formulae 
*----------------------------------------------------------------------*

      ! transformation
      call add_target2('FOPT_T',.false.,tgt_info)
      call set_dependency('FOPT_T','F_T',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_T','DEF_ME_Dtr',tgt_info)
      call set_rule2('FOPT_T',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T'/))
      call set_arg('FOPT_T',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T'/))

      ! transformation, generic form
      call add_target2('FOPT_Ttr_GEN',.false.,tgt_info)
      call set_dependency('FOPT_Ttr_GEN','F_Ttr_GEN',tgt_info)
      call set_dependency('FOPT_Ttr_GEN','DEF_ME_Tout',tgt_info)
      call set_dependency('FOPT_Ttr_GEN','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_Ttr_GEN','DEF_ME_Dtr',tgt_info)
      call set_rule2('FOPT_Ttr_GEN',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Ttr_GEN',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Ttr_GEN'/))
      call set_arg('FOPT_Ttr_GEN',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Ttr_GEN'/))

      ! transformed Hessian
      call add_target2('FOPT_Atr',.false.,tgt_info)
      call set_dependency('FOPT_Atr','F_Atr',tgt_info)
      call set_dependency('FOPT_Atr',mel_ham,tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_A',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_Dtr',tgt_info)
      call set_dependency('FOPT_Atr','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_Atr',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_Atr',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('FOPT_Atr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      call set_rule2('FOPT_Atr',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Atr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Atr'/))
      call set_arg('FOPT_Atr',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Atr'/))

      ! transformed Hessian times vector
      call add_target2('FOPT_A_Ttr',.false.,tgt_info)
      call set_dependency('FOPT_A_Ttr','F_A_Ttr',tgt_info)
      call set_dependency('FOPT_A_Ttr',mel_ham,tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_OMGtr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_Ttr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_1',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_Dtr',tgt_info)
      call set_dependency('FOPT_A_Ttr','DEF_ME_C0',tgt_info)
      call set_rule2('FOPT_A_Ttr',ASSIGN_ME2OP,tgt_info)
      call set_arg('FOPT_A_Ttr',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('FOPT_A_Ttr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      call set_rule2('FOPT_A_Ttr',OPTIMIZE,tgt_info)
      call set_arg('FOPT_A_Ttr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_A_Ttr'/))
      call set_arg('FOPT_A_Ttr',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_A_Ttr'/))

      call add_target3((/
     &   'target FOPT_H0Dy(',
     &   '  depend (MELS_FOR_H0INV,F_H0Dy,DEF_ME_C00)',
     &   '  OPTIMIZE(LABEL_OPT=FOPT_H0Dy,',
     &   '           LABELS_IN=F_H0Dy,INTERM=F_D00))'
     &     /),tgt_info)

      ! extended preconditioner
      call add_target3((/
     &   'target FOPT_H0INV(',
     &   '  depend (MELS_FOR_H0INV,DEF_ME_Dproj,',
     &   '          DEF_ME_C0,F_FOR_H0INV,DEF_ME_C00)',
     &   '  OPTIMIZE(LABEL_OPT=FOPT_H0INV,',
     &   '      LABELS_IN=(F_OMGprj,F_H0_DlT,F_S_DlT),INTERM=F_D00))'
     &     /),tgt_info)

      ! Residual
      call add_target2('FOPT_OMG',.false.,tgt_info)
      call set_dependency('FOPT_OMG','F_OMG',tgt_info)
      call set_dependency('FOPT_OMG','F_MRCC_E',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_T',tgt_info)
      if (Op_eqs) then
        call set_dependency('FOPT_OMG','F_Heff',tgt_info)
        call set_dependency('FOPT_OMG','F_Geff',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_Heff',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_Geff',tgt_info)
      end if
      if (.false..and.maxh.gt.0)
     &    call set_dependency('FOPT_OMG','DEF_ME_TT',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT1',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_HT2',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_OMG',tgt_info)
      call set_dependency('FOPT_OMG','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_OMG',mel_ham,tgt_info)
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_dependency('FOPT_OMG','F_OMG_C0',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_A_C0',tgt_info)
c      call set_dependency('FOPT_OMG','DEF_ME_1v',tgt_info)
        call set_rule2('FOPT_OMG',ASSIGN_ME2OP,tgt_info)
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_A_C0'/))
        call set_arg('FOPT_OMG',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &             val_label=(/'A_C0'/))
      end if
      if (tfix.gt.0)
     &    call set_dependency('FOPT_OMG','DEF_ME_Tfix',tgt_info)
      call set_rule2('FOPT_OMG',OPTIMIZE,tgt_info)
      call set_arg('FOPT_OMG',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_OMG'/))
      if ((maxp.ge.2.or.maxh.ge.2).and.tfix.eq.0) then
        labels(1:20)(1:len_target_name) = ' '
        ndef = 0
        if (maxp.ge.2.and.h1bar_maxp.lt.4) then
          call set_dependency('FOPT_OMG','F_PP0int',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_INT_PP0',tgt_info)
          ndef = ndef + 1
          labels(ndef) = 'F_PP0int'
        end if
        if (maxp.ge.2.and.(maxh.gt.0.or.h1bar_maxp.gt.2)) then
          call set_dependency('FOPT_OMG','F_PPint',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_INT_PP',tgt_info)
          ndef = ndef + 1
          labels(ndef) = 'F_PPint'
        end if
        if (maxh.ge.2) then
          call set_dependency('FOPT_OMG','F_HHint',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_INT_HH',tgt_info)
          ndef = ndef + 1
          labels(ndef) = 'F_HHint'
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
        call set_arg('FOPT_OMG',OPTIMIZE,'INTERM',ndef,tgt_info,
     &               val_label=labels(1:ndef))
      end if

      labels(1:20)(1:len_target_name) = ' '
      ndef = 0
      if (maxcum.gt.0) then
        call set_dependency('FOPT_OMG','F_DENS0',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_DENS',tgt_info)
        labels(ndef+1) = 'F_DENS0'
        ndef = ndef + 1
        if (cum_appr_mode.eq.0) then
          call set_dependency('FOPT_OMG','F_CENT',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_CENT',tgt_info)
          call set_dependency('FOPT_OMG','F_CUM',tgt_info)
          call set_dependency('FOPT_OMG','DEF_ME_CUM',tgt_info)
          labels(ndef+1) = 'F_CENT'
          labels(ndef+2) = 'F_CUM'
          ndef = ndef + 2
        end if
        if (maxcom.le.2.and.cum_appr_mode.ge.2) then
          select case(maxcum)
          case(2)
            call set_dependency('FOPT_OMG','F_D_INT08',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT09',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT10',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT11',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT12',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT08',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT09',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT10',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT11',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT12',tgt_info)
            labels(ndef+1) = 'F_D_INT08'
            labels(ndef+2) = 'F_D_INT09'
            labels(ndef+3) = 'F_D_INT10'
            labels(ndef+4) = 'F_D_INT11'
            labels(ndef+5) = 'F_D_INT12'
            ndef = ndef + 5
          case(3)
            call set_dependency('FOPT_OMG','F_D_INT04',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT05',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT06',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT07',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT04',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT05',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT06',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT07',tgt_info)
            labels(ndef+1) = 'F_D_INT04'
            labels(ndef+2) = 'F_D_INT05'
            labels(ndef+3) = 'F_D_INT06'
            labels(ndef+4) = 'F_D_INT07'
            ndef = ndef + 4
          case(4)
            call set_dependency('FOPT_OMG','F_D_INT01',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT02',tgt_info)
            call set_dependency('FOPT_OMG','F_D_INT03',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT01',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT02',tgt_info)
            call set_dependency('FOPT_OMG','DEF_ME_D_INT03',tgt_info)
            labels(ndef+1) = 'F_D_INT01'
            labels(ndef+2) = 'F_D_INT02'
            labels(ndef+3) = 'F_D_INT03'
            ndef = ndef + 3
            if (cum_appr_mode.eq.3) then
              call set_dependency('FOPT_OMG','F_D_INT13',tgt_info)
              call set_dependency('FOPT_OMG','DEF_ME_D_INT13',tgt_info)
              labels(ndef+1) = 'F_D_INT13'
              ndef = ndef + 1
            end if
          case default
          end select
        end if
      end if
      if (Op_eqs) then
        labels(ndef+1) = 'F_Heff'
        ndef = ndef + 1
      end if
      if (h1bar) then
        call set_dependency('FOPT_OMG','F_H1bar',tgt_info)
        call set_dependency('FOPT_OMG','DEF_ME_H1bar',tgt_info)
        labels(ndef+1) = 'F_H1bar'
        ndef = ndef + 1
      end if
      labels(ndef+1) = 'F_MRCC_E'
      ndef = ndef + 1
      if (Op_eqs) then
        labels(ndef+1) = 'F_Geff'
        ndef = ndef + 1
      end if
      labels(ndef+1) = 'F_OMG'
      ndef = ndef + 1
      if (optref.eq.-1.or.optref.eq.-2) then
        labels(ndef+1) = 'F_OMG_C0'
        ndef = ndef + 1
      end if
      call set_arg('FOPT_OMG',OPTIMIZE,'LABELS_IN',ndef,tgt_info,
     &             val_label=labels(1:ndef))

      ! Residual for C0
      call add_target2('FOPT_OMG_C0',.false.,tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_OMG_C0','DEF_ME_T',tgt_info)
c      call set_dependency('FOPT_OMG_C0','DEF_ME_1v',tgt_info)
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

c dbg
c      ! transformed metric
c      call add_target2('FOPT_Str',.false.,tgt_info)
c      call set_dependency('FOPT_Str','F_Str',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_S',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_1',tgt_info)
c      call set_dependency('FOPT_Str','DEF_ME_Dtr',tgt_info)
c      call set_rule2('FOPT_Str',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_Str',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_Str'/))
c      call set_arg('FOPT_Str',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_Str'/))
c
c      ! transformed metric times vector
c      call add_target2('FOPT_S_Ttr',.false.,tgt_info)
c      call set_dependency('FOPT_S_Ttr','F_S_Ttr',tgt_info)
c      call set_dependency('FOPT_S_Ttr','DEF_ME_OMG',tgt_info) !tr',tgt_info)
c      call set_dependency('FOPT_S_Ttr','DEF_ME_Ttr',tgt_info)
c      call set_dependency('FOPT_S_Ttr','DEF_ME_1',tgt_info)
c      call set_dependency('FOPT_S_Ttr','DEF_ME_Dtr',tgt_info)
c      call set_dependency('FOPT_S_Ttr','DEF_ME_C0',tgt_info)
c      call set_rule2('FOPT_S_Ttr',ASSIGN_ME2OP,tgt_info)
c      call set_arg('FOPT_S_Ttr',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_Dtr'/))
c      call set_arg('FOPT_S_Ttr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'Dtr'/))
c      call set_rule2('FOPT_S_Ttr',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_S_Ttr',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_S_Ttr'/))
c      call set_arg('FOPT_S_Ttr',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_S_Ttr'/))
c dbgend

c dbg
c      ! effective Hamiltonian
c      call add_target2('FOPT_Heff',.false.,tgt_info)
c      call set_dependency('FOPT_Heff','F_Heff',tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_1v',tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_Heff',tgt_info)
c      call set_dependency('FOPT_Heff',mel_ham,tgt_info)
c      call set_dependency('FOPT_Heff','DEF_ME_T',tgt_info)
c      call set_rule2('FOPT_Heff',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_Heff',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_Heff'/))
c      call set_arg('FOPT_Heff',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_Heff'/))
c dbgend

      ! "redundant" T components
      call add_target2('FOPT_T(2)red',.false.,tgt_info)
      call set_dependency('FOPT_T(2)red','F_T(2)red',tgt_info)
      call set_dependency('FOPT_T(2)red','DEF_ME_T(2)red',tgt_info)
      call set_dependency('FOPT_T(2)red','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T(2)red','DEF_ME_Dtr',tgt_info)
      call set_rule2('FOPT_T(2)red',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T(2)red',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T(2)red'/))
      call set_arg('FOPT_T(2)red',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T(2)red'/))
      call add_target2('FOPT_T(3)red',.false.,tgt_info)
      call set_dependency('FOPT_T(3)red','F_T(3)red',tgt_info)
      call set_dependency('FOPT_T(3)red','DEF_ME_T(3)red',tgt_info)
      call set_dependency('FOPT_T(3)red','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T(3)red','DEF_ME_Dtr',tgt_info)
      call set_rule2('FOPT_T(3)red',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T(3)red',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T(3)red'/))
      call set_arg('FOPT_T(3)red',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T(3)red'/))

      ! Spin "expectation value"
      call add_target2('FOPT_MRCC_S(S+1)',.false.,tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','F_MRCC_S(S+1)',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_S(S+1)',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_S+',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_S-',tgt_info)
      call set_dependency('FOPT_MRCC_S(S+1)','DEF_ME_Sz',tgt_info)
      call set_rule2('FOPT_MRCC_S(S+1)',OPTIMIZE,tgt_info)
      call set_arg('FOPT_MRCC_S(S+1)',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_MRCC_S(S+1)'/))
      call set_arg('FOPT_MRCC_S(S+1)',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_MRCC_S(S+1)'/))

      ! Spin "expectation value" <C0| T^+ S^2 T |C0>
      call add_target2('FOPT_T_S2',.false.,tgt_info)
      call set_dependency('FOPT_T_S2','F_T_S2',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_S(S+1)',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_S+',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_S-',tgt_info)
      call set_dependency('FOPT_T_S2','DEF_ME_Sz',tgt_info)
      call set_rule2('FOPT_T_S2',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T_S2',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T_S2'/))
      call set_arg('FOPT_T_S2',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T_S2'/))

      ! Norm <C0| T^+ T |C0>
      call add_target2('FOPT_T_NORM',.false.,tgt_info)
      call set_dependency('FOPT_T_NORM','F_T_NORM',tgt_info)
      call set_dependency('FOPT_T_NORM','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_T_NORM','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_T_NORM','DEF_ME_NORM',tgt_info)
      call set_rule2('FOPT_T_NORM',OPTIMIZE,tgt_info)
      call set_arg('FOPT_T_NORM',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_T_NORM'/))
      call set_arg('FOPT_T_NORM',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_T_NORM'/))

      ! Energy with Lagrangian based corrections
      call add_target2('FOPT_Ecorrected',.false.,tgt_info)
      call set_dependency('FOPT_Ecorrected','F_Ecorrected',tgt_info)
      call set_dependency('FOPT_Ecorrected','DEF_ME_C0',tgt_info)
      call set_dependency('FOPT_Ecorrected','DEF_ME_T',tgt_info)
      call set_dependency('FOPT_Ecorrected','DEF_ME_E(MR)',tgt_info)
      call set_dependency('FOPT_Ecorrected',mel_ham,tgt_info)
      if ((maxp.ge.2.or.maxh.ge.2).and.tfix.eq.0) then
        if (h1bar_maxp.lt.4)
     &     call set_dependency('FOPT_Ecorrected',
     &                         'DEF_ME_INT_PP0',tgt_info)
        if (maxp.ge.2.and.(maxh.gt.0.or.h1bar_maxp.gt.2)) 
     &     call set_dependency('FOPT_Ecorrected',
     &                         'DEF_ME_INT_PP',tgt_info)
        if (maxh.ge.2)
     &     call set_dependency('FOPT_Ecorrected',
     &                         'DEF_ME_INT_HH',tgt_info)
      end if
      if (h1bar)
     &   call set_dependency('FOPT_Ecorrected','DEF_ME_H1bar',tgt_info)
      if (tfix.gt.0) then
        call set_dependency('FOPT_Ecorrected','DEF_ME_Tfix',tgt_info)
        inquire(file='ME_Tfix_list.da',exist=l_exist)
        if (.not.l_exist) call quit(1,'set_ic_mrcc_targets',
     &           'File for fixed T amplitudes not found!')
      end if
      call set_rule2('FOPT_Ecorrected',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Ecorrected',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Ecorrected'/))
      call set_arg('FOPT_Ecorrected',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Ecorrected'/))

      ! formula for modification of trafo matrix (add. trafo to GNO)
      call add_target2('FOPT_Dinv_GNO',.false.,tgt_info)
      call set_dependency('FOPT_Dinv_GNO','F_Dinv_GNO',tgt_info)
      call set_dependency('FOPT_Dinv_GNO','DEF_ME_D',tgt_info)
      if (gno.gt.0) call set_dependency('FOPT_Dinv_GNO','DEF_ME_CUM',
     &                                  tgt_info)
      call set_rule2('FOPT_Dinv_GNO',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Dinv_GNO',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Dinv_GNO'/))
      call set_arg('FOPT_Dinv_GNO',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Dinv_GNO'/))

      ! formula for modification of Projector (add. trafo to GNO)
      call add_target2('FOPT_Dproj_GNO',.false.,tgt_info)
      call set_dependency('FOPT_Dproj_GNO','F_Dproj_GNO',tgt_info)
      call set_dependency('FOPT_Dproj_GNO','DEF_ME_D',tgt_info)
      if (gno.gt.0) call set_dependency('FOPT_Dproj_GNO','DEF_ME_CUM',
     &                                  tgt_info)
      call set_rule2('FOPT_Dproj_GNO',OPTIMIZE,tgt_info)
      call set_arg('FOPT_Dproj_GNO',OPTIMIZE,'LABEL_OPT',1,tgt_info,
     &             val_label=(/'FOPT_Dproj_GNO'/))
      call set_arg('FOPT_Dproj_GNO',OPTIMIZE,'LABELS_IN',1,tgt_info,
     &             val_label=(/'F_Dproj_GNO'/))
*----------------------------------------------------------------------*
*     ME-lists
*----------------------------------------------------------------------*

      ! Diagonal Preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call add_target(trim(dia_label),ttype_opme,.false.,tgt_info)
      call set_dependency(trim(dia_label),'EVAL_FREF',tgt_info)
      call set_dependency(trim(dia_label),
     &                    trim(op_dia)//'_'//'T',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = trim(dia_label)
      labels(2) = trim(op_dia)//'_'//'T'
      call me_list_parameters(-1,parameters,
     &     0,0,1,
     &     0,0,.false.)
      call set_rule(trim(dia_label),ttype_opme,
     &              DEF_ME_LIST,
     &     labels,2,1,
     &     parameters,1,tgt_info)
      if (prc_type.eq.2) then
        call set_dependency(trim(dia_label),'DEF_ME_Ttr',tgt_info)
        labels(1) = 'ME_Ttr'
        call set_dependency(trim(dia_label),'DEF_ME_DENS1',tgt_info)
        call set_rule2(trim(dia_label),PRECONDITIONER,tgt_info)
        call set_arg(trim(dia_label),PRECONDITIONER,'LIST_PRC',1,
     &       tgt_info,val_label=(/'ME_Ttr'/))
        call set_arg(trim(dia_label),PRECONDITIONER,'LIST_INP',2,
     &       tgt_info,val_label=(/'ME_FREF ','ME_DENS1'/))
        call set_arg(trim(dia_label),PRECONDITIONER,'MODE',1,tgt_info,
     &       val_str='dia-Fshift')
        call set_arg(trim(dia_label),PRECONDITIONER,'SHIFT',1,tgt_info,
     &       val_rl8=(/prc_shift/))
      else if (prc_type.ge.0) then
        call set_rule2(trim(dia_label),PRECONDITIONER,tgt_info)
        call set_arg(trim(dia_label),PRECONDITIONER,'LIST_PRC',1,
     &       tgt_info,val_label=(/trim(dia_label)/))
        call set_arg(trim(dia_label),PRECONDITIONER,'LIST_INP',1,
     &       tgt_info,val_label=(/'ME_FREF'/))
        call set_arg(trim(dia_label),PRECONDITIONER,'MODE',1,tgt_info,
     &       val_str='dia-F')
      end if
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (a):',0,'LIST')
c      call set_rule(trim(dia_label),ttype_opme,PRINT_MEL,
c     &     trim(labels(1)),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! ME for T
      call add_target2('DEF_ME_T',.false.,tgt_info)
      call set_dependency('DEF_ME_T','T',tgt_info)
      call set_rule2('DEF_ME_T',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_T',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_T'/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'T   '/))
      call set_arg('DEF_ME_T',DEF_ME_LIST,'2MS',1,tgt_info,
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
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Ttr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! ME for Tout
      call add_target2('DEF_ME_Tout',.false.,tgt_info)
      call set_dependency('DEF_ME_Tout','Tout',tgt_info)
      call set_rule2('DEF_ME_Tout',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Tout',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Tout'/))
      call set_arg('DEF_ME_Tout',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Tout'/))
      call set_arg('DEF_ME_Tout',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Tout',DEF_ME_LIST,'IRREP',1,tgt_info,
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
        call set_arg(def_ht,DEF_ME_LIST,'2MS',1,tgt_info,
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
      call set_arg('DEF_ME_TT',DEF_ME_LIST,'2MS',1,tgt_info,
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
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_OMG',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for transformed Residual
      call add_target2('DEF_ME_OMGtr',.false.,tgt_info)
      call set_dependency('DEF_ME_OMGtr','OMGtr',tgt_info)
      call set_rule2('DEF_ME_OMGtr',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMGtr'/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'OMGtr'/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_OMGtr',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))

      ! reordered projector matrix ME_Dproj (to eliminate lin. dep.)
      call add_target('DEF_ME_Dproj',ttype_opme,.false.,tgt_info)
!      call set_dependency('DEF_ME_Dproj','Dtr',tgt_info)
      call set_dependency('DEF_ME_Dproj','Dproj',tgt_info)
      labels(1:20)(1:len_target_name) = ' '
      labels(1) = 'ME_Dproj'
!      labels(2) = 'Dtr'
      labels(2) = 'Dproj'
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

c dbg
      ! ME_Heff
      call add_target2('DEF_ME_Heff',.false.,tgt_info)
      call set_dependency('DEF_ME_Heff','Heff',tgt_info)
      call set_rule2('DEF_ME_Heff',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'LIST',1,tgt_info,
     &     val_label=(/'ME_Heff'/))
      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &     val_label=(/'Heff'/))
      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'IRREP',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'2MS',1,tgt_info,
     &     val_int=(/0/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
c     &     val_int=(/orb_info%lsym/))
c      call set_arg('DEF_ME_Heff',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
c     &     val_int=(/mod(orb_info%imult-1,2)/))
c
c      ! ME_1v
c      call add_target2('DEF_ME_1v',.false.,tgt_info)
c      call set_dependency('DEF_ME_1v','1v',tgt_info)
c      call set_rule2('DEF_ME_1v',DEF_ME_LIST,tgt_info)
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'LIST',1,tgt_info,
c     &     val_label=(/'ME_1v'/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'OPERATOR',1,tgt_info,
c     &     val_label=(/'1v'/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'IRREP',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'2MS',1,tgt_info,
c     &     val_int=(/0/))
c      call set_arg('DEF_ME_1v',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
c     &     val_int=(/1/))
c      call dens_parameters(-1,parameters,0,0,0)
c      call set_rule('DEF_ME_1v',ttype_opme,UNITY,
c     &     'ME_1v',1,1,
c     &     parameters,1,tgt_info)

      ! ME for Geff
      call add_target2('DEF_ME_Geff',.false.,tgt_info)
      call set_dependency('DEF_ME_Geff','Geff',tgt_info)
      call set_rule2('DEF_ME_Geff',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Geff',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Geff'/))
      call set_arg('DEF_ME_Geff',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Geff'/))
      call set_arg('DEF_ME_Geff',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Geff',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Geff',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
c dbgend

      ! ME for H1bar
      call add_target2('DEF_ME_H1bar',.false.,tgt_info)
      call set_dependency('DEF_ME_H1bar','H1bar',tgt_info)
      call set_rule2('DEF_ME_H1bar',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_H1bar',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_H1bar'/))
      call set_arg('DEF_ME_H1bar',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'H1bar'/))
      call set_arg('DEF_ME_H1bar',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_H1bar',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_H1bar',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for one-particle density (for precond.)
      call add_target2('DEF_ME_DENS1',.false.,tgt_info)
      call set_dependency('DEF_ME_DENS1','DENS1',tgt_info)
      call set_dependency('DEF_ME_DENS1','EVAL_D',tgt_info)
      call set_rule2('DEF_ME_DENS1',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_DENS1',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_DENS1'/))
      call set_arg('DEF_ME_DENS1',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'DENS1'/))
      call set_arg('DEF_ME_DENS1',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_DENS1',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_rule2('DEF_ME_DENS1',SCALE_COPY,tgt_info)
      call set_arg('DEF_ME_DENS1',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_DENS1'/))
      call set_arg('DEF_ME_DENS1',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_DENS'/))
      call set_arg('DEF_ME_DENS1',SCALE_COPY,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))

      ! MEs for H0 Inversion
      call add_target3((/
     &     'target MAKE_H0Dy(',
     &     '    depend (MELS_FOR_H0INV,FOPT_H0Dy)',
     &     '    EVALUATE(form=FOPT_H0Dy)',
     &     '    PRINT_MEL(list=ME_H0Dy)',
     &     ')'
     &     /),tgt_info)
      call add_target3((/
     &     'target MELS_FOR_H0INV(',
     &     '    depend OPS_FOR_H0INV',
     &     '    DEF_ME_LIST(LIST=ME_OMGprj,OPERATOR=OMGprj,',
     &     '                  2MS=0,IRREP=1)',          
     &     '    DEF_ME_LIST(LIST=ME_DlT,OPERATOR=DlT,',
     &     '                  2MS=0,IRREP=1)',          
     &     '    DEF_ME_LIST(LIST=ME_H0_DlT,OPERATOR=H0_DlT,',
     &     '                  2MS=0,IRREP=1)',
     &     '    DEF_ME_LIST(LIST=ME_S_DlT,OPERATOR=S_DlT,',
     &     '                  2MS=0,IRREP=1)',
     &     '    DEF_ME_LIST(LIST=ME_D00,OPERATOR=D00,',
     &     '                  2MS=0,IRREP=1)',
     &     '    DEF_ME_LIST(LIST=ME_H0Dy,OPERATOR=H0Dy,',
     &     '                  2MS=0,IRREP=1)',
     &     ')'/),
     &     tgt_info)

      ! ME for Intermediate(s)
      call add_target2('DEF_ME_INT_PP',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_PP','INT_PP',tgt_info)
      call set_rule2('DEF_ME_INT_PP',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_PP',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_PP'/))
      call set_arg('DEF_ME_INT_PP',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_PP'/))
      call set_arg('DEF_ME_INT_PP',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_PP',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_PP',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call add_target2('DEF_ME_INT_PP0',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_PP0','INT_PP0',tgt_info)
      call set_rule2('DEF_ME_INT_PP0',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_PP0',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_PP0'/))
      call set_arg('DEF_ME_INT_PP0',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_PP0'/))
      call set_arg('DEF_ME_INT_PP0',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_PP0',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_PP0',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      call add_target2('DEF_ME_INT_HH',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_HH','INT_HH',tgt_info)
      call set_rule2('DEF_ME_INT_HH',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_HH',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_HH'/))
      call set_arg('DEF_ME_INT_HH',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_HH'/))
      call set_arg('DEF_ME_INT_HH',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_HH',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_HH',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for INT_HT2
      call add_target2('DEF_ME_INT_HT2',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_HT2','INT_HT2',tgt_info)
      call set_rule2('DEF_ME_INT_HT2',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_HT2',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_HT2'/))
      call set_arg('DEF_ME_INT_HT2',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_HT2'/))
      call set_arg('DEF_ME_INT_HT2',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_HT2',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_HT2',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      ! ME for INT_T2H
      call add_target2('DEF_ME_INT_T2H',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_T2H','INT_T2H',tgt_info)
      call set_rule2('DEF_ME_INT_T2H',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_T2H',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_T2H'/))
      call set_arg('DEF_ME_INT_T2H',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_T2H'/))
      call set_arg('DEF_ME_INT_T2H',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_T2H',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_T2H',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      ! ME for INT_D
      call add_target2('DEF_ME_INT_D',.false.,tgt_info)
      call set_dependency('DEF_ME_INT_D','INT_D',tgt_info)
      call set_rule2('DEF_ME_INT_D',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_INT_D',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_INT_D'/))
      call set_arg('DEF_ME_INT_D',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'INT_D'/))
      call set_arg('DEF_ME_INT_D',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_INT_D',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_INT_D',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for T(2)red
      call add_target2('DEF_ME_T(2)red',.false.,tgt_info)
      call set_dependency('DEF_ME_T(2)red','T',tgt_info)
      call set_rule2('DEF_ME_T(2)red',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_T(2)red',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_T(2)red'/))
      call set_arg('DEF_ME_T(2)red',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'T(2)red'/))
      call set_arg('DEF_ME_T(2)red',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_T(2)red',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_T(2)red',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))
      ! ME for T(3)red
      call add_target2('DEF_ME_T(3)red',.false.,tgt_info)
      call set_dependency('DEF_ME_T(3)red','T',tgt_info)
      call set_rule2('DEF_ME_T(3)red',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_T(3)red',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_T(3)red'/))
      call set_arg('DEF_ME_T(3)red',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'T(3)red'/))
      call set_arg('DEF_ME_T(3)red',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_T(3)red',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_T(3)red',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

      ! ME for Tfix
      call add_target2('DEF_ME_Tfix',.false.,tgt_info)
      call set_dependency('DEF_ME_Tfix','Tfix',tgt_info)
      call set_rule2('DEF_ME_Tfix',DEF_ME_LIST,tgt_info)
      call set_arg('DEF_ME_Tfix',DEF_ME_LIST,'LIST',1,tgt_info,
     &             val_label=(/'ME_Tfix'/))
      call set_arg('DEF_ME_Tfix',DEF_ME_LIST,'OPERATOR',1,tgt_info,
     &             val_label=(/'Tfix'/))
      call set_arg('DEF_ME_Tfix',DEF_ME_LIST,'2MS',1,tgt_info,
     &             val_int=(/0/))
      call set_arg('DEF_ME_Tfix',DEF_ME_LIST,'IRREP',1,tgt_info,
     &             val_int=(/1/))
      call set_arg('DEF_ME_Tfix',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &             val_int=(/msc/))

c dbg
c      ! ME_S
c      call add_target2('DEF_ME_S',.false.,tgt_info)
c      call set_dependency('DEF_ME_S','S',tgt_info)
c      call set_rule2('DEF_ME_S',DEF_ME_LIST,tgt_info)
c      call set_arg('DEF_ME_S',DEF_ME_LIST,'LIST',1,tgt_info,
c     &     val_label=(/'ME_S'/))
c      call set_arg('DEF_ME_S',DEF_ME_LIST,'OPERATOR',1,tgt_info,
c     &     val_label=(/'S'/))
c      call set_arg('DEF_ME_S',DEF_ME_LIST,'IRREP',1,tgt_info,
c     &     val_int=(/1/))
c      call set_arg('DEF_ME_S',DEF_ME_LIST,'MS',1,tgt_info,
c     &     val_int=(/0/))
c dbgend
*----------------------------------------------------------------------*
*     "phony" targets: solve equations, evaluate expressions
*----------------------------------------------------------------------*

      ! Evaluate diagonal elements of Jacobian
      call add_target('EVAL_Atr',ttype_gen,.false.,tgt_info)
c dbg hybrid preconditioner
c      call set_dependency('EVAL_Atr','EVAL_A_Ttr',tgt_info)
c dbgend
      call set_dependency('EVAL_Atr','FOPT_Atr',tgt_info)
c      call set_dependency('EVAL_Atr','EVAL_FREF',tgt_info)
      call set_rule('EVAL_Atr',ttype_opme,EVAL,
     &     'FOPT_Atr',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'transformed Jacobian :',0,'LIST')
c      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
c     &     'ME_A',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      ! put diagonal elements to preconditioner
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('EVAL_Atr',trim(dia_label),tgt_info)
      call set_rule2('EVAL_Atr',EXTRACT_DIAG,tgt_info)
      call set_arg('EVAL_Atr',EXTRACT_DIAG,'LIST_RES',1,tgt_info,
     &             val_label=(/trim(dia_label)/))
      call set_arg('EVAL_Atr',EXTRACT_DIAG,'LIST_IN',1,tgt_info,
     &             val_label=(/'ME_A'/))
      if (prc_type.ge.3)
     &  call set_arg('EVAL_Atr',EXTRACT_DIAG,'EXTEND',1,tgt_info,
     &               val_log=(/.true./))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (b) :',0,'LIST')
c      call set_rule('EVAL_Atr',ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Evaluate approximation of diagonal elements of Jacobian
      call add_target('EVAL_A_Ttr',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_A_Ttr','FOPT_A_Ttr',tgt_info)
      ! (a) set all elements of right hand vector to one
      call set_rule2('EVAL_A_Ttr',SET_MEL,tgt_info)
      call set_arg('EVAL_A_Ttr',SET_MEL,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
      call set_arg('EVAL_A_Ttr',SET_MEL,'VAL_LIST',2,tgt_info,
     &             val_rl8=(/1d0,+1d0/))
      call set_arg('EVAL_A_Ttr',SET_MEL,'IDX_LIST',2,tgt_info,
     &             val_int=(/-1,-2/))
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Ttr set to :',0,'LIST')
c      call set_rule('EVAL_A_Ttr',ttype_opme,PRINT_MEL,
c     &     'ME_Ttr',1,0,
c     &     parameters,2,tgt_info)
c dbgend
      ! (b) evaluate
      call set_rule('EVAL_A_Ttr',ttype_opme,EVAL,
     &     'FOPT_A_Ttr',1,0,
     &     parameters,0,tgt_info)
      ! (c) scale back and copy onto preconditioner list
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('EVAL_A_Ttr',trim(dia_label),tgt_info)
      call set_rule2('EVAL_A_Ttr',SCALE_COPY,tgt_info)
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/trim(dia_label)/))
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_OMGtr'/))
      call set_arg('EVAL_A_Ttr',SCALE_COPY,'FAC',2,tgt_info,
     &             val_rl8=(/1d0,+1d0/))
      ! (d) free up Ttr list
      call set_rule('EVAL_A_Ttr',ttype_opme,RES_ME_LIST,
     &     'ME_Ttr',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (b) :',0,'LIST')
c      call set_rule('EVAL_A_Ttr',ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Preconditioner by transforming the diagonal pre-preconditioner
      call add_target('PREC_diag',ttype_gen,.false.,tgt_info)
      call set_dependency('PREC_diag','FOPT_T',tgt_info)
      call set_dependency('PREC_diag','DEF_ME_Dudag_2',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('PREC_diag',trim(dia_label),tgt_info)
      call set_rule2('PREC_diag',ASSIGN_ME2OP,tgt_info)
      call set_arg('PREC_diag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/trim(dia_label)/))
      call set_arg('PREC_diag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'T'/))
      call set_rule2('PREC_diag',ASSIGN_ME2OP,tgt_info)
      call set_arg('PREC_diag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dudag_2'/))
      call set_arg('PREC_diag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      call set_rule('PREC_diag',ttype_opme,EVAL,
     &     'FOPT_T',1,0,
     &     parameters,0,tgt_info)
      call set_rule2('PREC_diag',ASSIGN_ME2OP,tgt_info)
      call set_arg('PREC_diag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_T'/))
      call set_arg('PREC_diag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'T'/))
      call set_rule2('PREC_diag',ASSIGN_ME2OP,tgt_info)
      call set_arg('PREC_diag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dtr'/))
      call set_arg('PREC_diag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      labels(1) = trim(op_dia)//'_'//'T'
      call set_rule2('PREC_diag',ASSIGN_ME2OP,tgt_info)
      call set_arg('PREC_diag',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/trim(dia_label)/))
      call set_arg('PREC_diag',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/trim(labels(1))/))
      call set_rule('PREC_diag',ttype_opme,RES_ME_LIST,
     &     'ME_Ttr',1,0,
     &     parameters,0,tgt_info)
      call set_rule('PREC_diag',ttype_opme,RES_ME_LIST,
     &     'ME_Dudag_2',1,0,
     &     parameters,0,tgt_info)
c dbg
c      call form_parameters(-1,parameters,2,
c     &     'Preconditioner (b) :',0,'LIST')
c      call set_rule('PREC_diag',ttype_opme,PRINT_MEL,
c     &     trim(dia_label),1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Solve MR coupled cluster equations
      call add_target2('SOLVE_MRCC',solve,tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_REF_S(S+1)',tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_OMG',tgt_info)
      call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
      dia_label = trim(dia_label)//'_T'
      call set_dependency('SOLVE_MRCC',trim(dia_label),tgt_info)
      call set_dependency('SOLVE_MRCC','EVAL_D',tgt_info)
      call set_dependency('SOLVE_MRCC','DEF_ME_Dtrdag',tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_T',tgt_info)
      if (restart) ! project out redundant part (if sv_thr. changed)
     &   call set_dependency('SOLVE_MRCC','EVAL_Tproj',tgt_info)
      prc_mode_str = 'TRF'  ! use diagonal in orth. basis (for T part)
      select case(prc_type)
      case(-1) !do nothing: use old preconditioner file!
        call warn('set_ic_mrcc_targets','Using old preconditioner file')
      case(0,3)
c dbg hybrid preconditioner
c        call warn('set_ic_mrcc_targets','Using hybrid preconditioner')
c dbgend
        if (maxv.gt.0.or.prc_type.eq.0)
     &     call set_dependency('SOLVE_MRCC','EVAL_Atr',tgt_info)
      case(1)
        call set_dependency('SOLVE_MRCC','EVAL_A_Ttr',tgt_info)
      case(2)
        call set_dependency('SOLVE_MRCC','PREC_diag',tgt_info)
        ! we misused ME_Dtrdag. Now get correct one:
        call set_rule2('SOLVE_MRCC',REORDER_MEL,tgt_info)
        call set_arg('SOLVE_MRCC',REORDER_MEL,'LIST_RES',1,tgt_info,
     &               val_label=(/'ME_Dtrdag'/))
        call set_arg('SOLVE_MRCC',REORDER_MEL,'LIST_IN',1,tgt_info,
     &               val_label=(/'ME_Dinv'/))
        call set_arg('SOLVE_MRCC',REORDER_MEL,'FROMTO',1,tgt_info,
     &               val_int=(/13/))
        call set_arg('SOLVE_MRCC',REORDER_MEL,'ADJOINT',1,tgt_info,
     &               val_log=(/.true./))
      case(4)
        if (maxv.gt.0.or.prc_type.eq.0)
     &     call set_dependency('SOLVE_MRCC','EVAL_Atr',tgt_info)
        call set_dependency('SOLVE_MRCC','FOPT_H0INV',tgt_info)
        call set_dependency('SOLVE_MRCC','MAKE_H0Dy',tgt_info)
        call set_dependency('SOLVE_MRCC','FOPT_Ttr_GEN',tgt_info)
        prc_mode_str='IH0'
      case default
        call quit(1,'set_ic_mrcc_targets','unknown prc_type')
      end select
      if (optref.ne.0) then
        call me_list_label(dia_label2,mel_dia,orb_info%lsym,
     &                     0,0,0,.false.)
        dia_label2 = trim(dia_label2)//'C0'
        call set_dependency('SOLVE_MRCC',trim(dia_label2),tgt_info)
        if (optref.ne.-1.and.optref.ne.-2) 
     &     call set_dependency('SOLVE_MRCC','FOPT_OMG_C0',tgt_info)
        call set_dependency('SOLVE_MRCC','DEF_ME_Dproj',tgt_info)
        call set_dependency('SOLVE_MRCC','DEF_ME_C00',tgt_info)
      end if
      do icnt = 1, max(1,optref)
      call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
      if (optref.lt.0) then
        if (preopt) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &         val_label=(/'ME_T'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &         val_str=prc_mode_str)
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &         val_label=(/'ME_OMG'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &         val_label=(/trim(dia_label)/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &       val_label=(/'ME_E(MR)'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &       val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',1,tgt_info,
     &       val_label=(/'FOPT_T'/))
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &         val_label=(/'FOPT_OMG'/))
          call set_rule2('SOLVE_MRCC',SOLVENLEQ,tgt_info)
        end if
      end if
      if (optref.eq.-1.or.optref.eq.-2) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',2,tgt_info,
     &       val_label=(/'ME_T ','ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str=prc_mode_str//'/NRM')
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',2,tgt_info,
     &       val_label=(/'ME_OMG ','ME_A_C0'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',2,tgt_info,
     &       val_label=(/trim(dia_label),trim(dia_label2)/))
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_T'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'MODE',1,tgt_info,
     &       val_str=prc_mode_str)
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_RESID',1,tgt_info,
     &       val_label=(/'ME_OMG'/))
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_PRC',1,tgt_info,
     &       val_label=(/trim(dia_label)/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_E',1,tgt_info,
     &     val_label=(/'ME_E(MR)'/))
      if (optref.ne.0.and.update_prc) then
        if (tred.eq.0) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr   ','ME_Dtr   ','ME_Dtrdag',
     &                 'ME_Dproj ',
     &                 'ME_D     ','ME_Dinv  ',
     &                 'ME_A     '/))
        else if (ex_t3red) then
        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
        call set_dependency('SOLVE_MRCC','FOPT_T(3)red',tgt_info)
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',9,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   ',
     &                 'ME_A      ','ME_T(2)red','ME_T(3)red'/))
        else
        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',8,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   ',
     &                 'ME_A      ','ME_T(2)red'/))
        end if
      else if (optref.ne.0.and..not.update_prc) then
        if (tred.eq.0) then
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',6,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   '/))
        else if (ex_t3red) then
        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
        call set_dependency('SOLVE_MRCC','FOPT_T(3)red',tgt_info)
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',8,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   ','ME_T(2)red',
     &                 'ME_T(3)red'/))
        else
        call set_dependency('SOLVE_MRCC','FOPT_T(2)red',tgt_info)
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',7,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag ',
     &                 'ME_Dproj  ',
     &                 'ME_D      ','ME_Dinv   ','ME_T(2)red'/))
        end if
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'LIST_SPC',3,tgt_info,
     &     val_label=(/'ME_Ttr    ','ME_Dtr    ','ME_Dtrdag '/))
      end if
      if (optref.ne.0) then
        if (update_prc) then
          if (tred.eq.0) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',3,tgt_info,
     &         val_label=(/'FOPT_T  ','FOPT_D  ','FOPT_Atr'/))
          else if (ex_t3red) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',5,tgt_info,
     &         val_label=(/'FOPT_T      ','FOPT_D      ',
     &                     'FOPT_Atr    ',
     &                     'FOPT_T(2)red','FOPT_T(3)red'/))
          else
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',4,tgt_info,
     &         val_label=(/'FOPT_T      ','FOPT_D      ',
     &                     'FOPT_Atr    ',
     &                     'FOPT_T(2)red'/))
          end if
        else
          if (tred.eq.0) then
            if (gno.gt.0) then
              call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',4,tgt_info,
     &             val_label=(/'FOPT_T        ','FOPT_D        ',
     &                         'FOPT_Dinv_GNO ','FOPT_Dproj_GNO'/))
            else
              call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',2,tgt_info,
     &             val_label=(/'FOPT_T','FOPT_D'/))
            end if
          else if (ex_t3red) then
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',4,tgt_info,
     &         val_label=(/'FOPT_T      ','FOPT_D      ',
     &                     'FOPT_T(2)red','FOPT_T(3)red'/))
          else
          call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',3,tgt_info,
     &         val_label=(/'FOPT_T      ','FOPT_D      ',
     &                     'FOPT_T(2)red'/))
          end if
        end if
      else
        call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM_SPC',1,tgt_info,
     &       val_label=(/'FOPT_T'/))
      end if
      call set_arg('SOLVE_MRCC',SOLVENLEQ,'FORM',1,tgt_info,
     &     val_label=(/'FOPT_OMG'/))
c dbg
c      ! print effective Hamiltonian
c      ! Evaluate transformed metric
c      call set_dependency('SOLVE_MRCC','FOPT_Heff',tgt_info)
c      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
c     &     'FOPT_Heff',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'effective Hamiltonian :',0,'LIST')
c      call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &     'ME_Heff',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      if (optref.gt.0.and.icnt.ne.optref) then !not in last iteration
        call set_rule2('SOLVE_MRCC',SOLVEEVP,tgt_info)
        call set_arg('SOLVE_MRCC',SOLVEEVP,'LIST_OPT',1,tgt_info,
     &       val_label=(/'ME_C0'/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'MODE',1,tgt_info,
     &       val_str='DIA')
        call set_arg('SOLVE_MRCC',SOLVEEVP,'N_ROOTS',1,tgt_info,
     &       val_int=(/maxroot/))
        call set_arg('SOLVE_MRCC',SOLVEEVP,'TARG_ROOTS',1,tgt_info,
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
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call set_rule2('SOLVE_MRCC',PRINT_MEL,tgt_info)
        call set_arg('SOLVE_MRCC',PRINT_MEL,'LIST',1,tgt_info,
     &       val_label=(/'ME_S(S+1)'/))
        call set_arg('SOLVE_MRCC',PRINT_MEL,'COMMENT',1,tgt_info,
     &       val_str='Spin expectation value <C0| S^2 |C0> :')
        call set_arg('SOLVE_MRCC',PRINT_MEL,'FORMAT',1,tgt_info,
     &       val_str='SCAL F20.12')
        call set_arg('SOLVE_MRCC',PRINT_MEL,'CHECK_THRESH',1,tgt_info,
     &       val_rl8=(/1d-2/))
        call set_arg('SOLVE_MRCC',PRINT_MEL,'EXPECTED',1,tgt_info,
     &       val_rl8=(/(dble(orb_info%imult**2)-1d0)/4d0/))
c        call form_parameters(-1,parameters,2,
c     &       'Spin expectation value <C0| S^2 |C0> :',0,'SCAL F20.12')
c        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
c     &       'S(S+1)',1,0,
c     &       parameters,2,tgt_info)
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
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_REF_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call form_parameters(-1,parameters,2,
     &       'Spin expectation value <C0| S^2 |C0> :',0,'SCAL F20.12')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_S(S+1)',1,0,
     &       parameters,2,tgt_info)
      end if
      call set_dependency('SOLVE_MRCC','FOPT_T_S2',tgt_info)
      call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
     &     'ME_S(S+1)',1,0,
     &     parameters,0,tgt_info)
      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &     'FOPT_T_S2',1,0,
     &     parameters,0,tgt_info)
      call set_dependency('SOLVE_MRCC','FOPT_T_NORM',tgt_info)
      call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
     &     'ME_NORM',1,0,
     &     parameters,0,tgt_info)
      call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &     'FOPT_T_NORM',1,0,
     &     parameters,0,tgt_info)
      call set_rule2('SOLVE_MRCC',SCALE_COPY,tgt_info)
      call set_arg('SOLVE_MRCC',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_S(S+1)'/))
      call set_arg('SOLVE_MRCC',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_NORM'/))
      call set_arg('SOLVE_MRCC',SCALE_COPY,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))
      call set_arg('SOLVE_MRCC',SCALE_COPY,'MODE',1,tgt_info,
     &             val_str='precond')
      call form_parameters(-1,parameters,2,
     &     'Spin expectation value <C0|T^+ S^2 T|C0>/<C0|T^+ T|C0> :',
     &     0,'SCAL F20.12')
      call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &     'ME_S(S+1)',1,0,
     &     parameters,2,tgt_info)
c dbg
      if (.false.) then
c dbg      if (.not.h1bar.and.tfix.eq.0.and.maxcum.le.0) then
        call set_dependency('SOLVE_MRCC','FOPT_MRCC_S(S+1)',tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,RES_ME_LIST,
     &       'ME_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call set_rule('SOLVE_MRCC',ttype_opme,EVAL,
     &       'FOPT_MRCC_S(S+1)',1,0,
     &       parameters,0,tgt_info)
        call form_parameters(-1,parameters,2,
     &       'Spin expectation value <C0| T^+ e^-T S^2 e^T |C0> :',
     &       0,'SCAL F20.12')
        call set_rule('SOLVE_MRCC',ttype_opme,PRINT_MEL,
     &       'ME_S(S+1)',1,0,
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

      ! Non-iterative higher-order correction
      call add_target2('EVAL_PERT_CORR',.not.svdonly.and.tfix.gt.0,
     &                 tgt_info)
      call set_dependency('EVAL_PERT_CORR','FOPT_Ecorrected',tgt_info)
      if (maxit.gt.1) then
        ! Use nonlinear solver
        call set_dependency('EVAL_PERT_CORR','SOLVE_MRCC',tgt_info)
      else
        ! Do first iteration without solver (saves virtual memory)
        call set_dependency('EVAL_PERT_CORR','FOPT_OMG',tgt_info)
        call me_list_label(dia_label,mel_dia,1,0,0,0,.false.)
        dia_label = trim(dia_label)//'_T'
        call set_dependency('EVAL_PERT_CORR',trim(dia_label),tgt_info)
        call set_dependency('EVAL_PERT_CORR','EVAL_D',tgt_info)
        call set_dependency('EVAL_PERT_CORR','DEF_ME_Dtrdag',tgt_info)
        call set_dependency('EVAL_PERT_CORR','FOPT_T',tgt_info)
        select case(prc_type)
        case(-1) !use old preconditioner file, but warn!
          call warn('set_ic_mrcc_targets',
     &       'Be sure you used prc_type=3 when creating precond. file!')
        case(3)
          call set_dependency('EVAL_PERT_CORR','EVAL_Atr',tgt_info)
        case default
          call quit(1,'set_ic_mrcc_targets',
     &         'Non-iterative higher-order corr. should use prc_type=3')
        end select
        ! (a) evaluate residual
        call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
     &       'FOPT_OMG',1,0,
     &       parameters,0,tgt_info)
        ! (b) transform residual
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dtrdag'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'Dtr'/))
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'T'/))
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMG'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'Ttr'/))
        call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
     &       'FOPT_T',1,0,
     &       parameters,0,tgt_info)
        ! (c) preconditioning step
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_OMG'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'OMG'/))
        call set_rule2('EVAL_PERT_CORR',SCALE_COPY,tgt_info)
        call set_arg('EVAL_PERT_CORR',SCALE_COPY,'LIST_RES',1,tgt_info,
     &               val_label=(/'ME_Ttr'/))
        call set_arg('EVAL_PERT_CORR',SCALE_COPY,'LIST_INP',1,tgt_info,
     &               val_label=(/trim(dia_label)/))
        call set_arg('EVAL_PERT_CORR',SCALE_COPY,'LIST_SHAPE',1,
     &               tgt_info,val_label=(/'ME_OMG'/))
        call set_arg('EVAL_PERT_CORR',SCALE_COPY,'FAC',1,tgt_info,
     &               val_rl8=(/-1d0/))
        call set_arg('EVAL_PERT_CORR',SCALE_COPY,'MODE',1,tgt_info,
     &               val_str='precond')
        ! (d) transform triples vector
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_Dtr'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'Dtr'/))
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_T'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'T'/))
        call set_rule2('EVAL_PERT_CORR',ASSIGN_ME2OP,tgt_info)
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
        call set_arg('EVAL_PERT_CORR',ASSIGN_ME2OP,'OPERATOR',1,
     &             tgt_info,val_label=(/'Ttr'/))
        call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
     &       'FOPT_T',1,0,
     &       parameters,0,tgt_info)
      end if
      call set_rule('EVAL_PERT_CORR',ttype_opme,RES_ME_LIST,
     &     'ME_E(MR)',1,0,
     &     parameters,0,tgt_info)
      call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
     &     'FOPT_Ecorrected',1,0,
     &     parameters,0,tgt_info)
      call form_parameters(-1,parameters,2,
     &     '>>> Total energy :',0,'SCAL F20.12')
      call set_rule('EVAL_PERT_CORR',ttype_opme,PRINT_MEL,
     &     'ME_E(MR)',1,0,
     &     parameters,2,tgt_info)
c dbg
c      ! Calculate and print <C0|T^+ S^2 T|C0>/<C0|S^2|C0>
c      call set_dependency('EVAL_PERT_CORR','FOPT_T_S2',tgt_info)
c      call set_rule('EVAL_PERT_CORR',ttype_opme,RES_ME_LIST,
c     &     'ME_S(S+1)',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
c     &     'FOPT_T_S2',1,0,
c     &     parameters,0,tgt_info)
c      call set_dependency('EVAL_PERT_CORR','FOPT_T_NORM',tgt_info)
c      call set_rule('EVAL_PERT_CORR',ttype_opme,RES_ME_LIST,
c     &     'ME_NORM',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('EVAL_PERT_CORR',ttype_opme,EVAL,
c     &     'FOPT_T_NORM',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule2('EVAL_PERT_CORR',SCALE_COPY,tgt_info)
c      call set_arg('EVAL_PERT_CORR',SCALE_COPY,'LIST_RES',1,tgt_info,
c     &             val_label=(/'ME_S(S+1)'/))
c      call set_arg('EVAL_PERT_CORR',SCALE_COPY,'LIST_INP',1,tgt_info,
c     &             val_label=(/'ME_NORM'/))
c      call set_arg('EVAL_PERT_CORR',SCALE_COPY,'FAC',1,tgt_info,
c     &             val_rl8=(/1d0/))
c      call set_arg('EVAL_PERT_CORR',SCALE_COPY,'MODE',1,tgt_info,
c     &             val_str='precond')
c      call form_parameters(-1,parameters,2,
c     &     'Spin expectation value '//
c     &                      '<C0|PT^+ S^2 PT|C0>/<C0|PT^+ PT|C0> :',
c     &     0,'SCAL F20.12')
c      call set_rule('EVAL_PERT_CORR',ttype_opme,PRINT_MEL,
c     &     'ME_S(S+1)',1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      ! Evaluate transformed metric
c      call add_target('EVAL_Str',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_Str','FOPT_Str',tgt_info)
c      call set_dependency('EVAL_Str','SOLVE_REF',tgt_info)
cc      call set_rule('EVAL_Str',ttype_opme,RES_ME_LIST,
cc     &     'ME_A',1,0,
cc     &     parameters,0,tgt_info)
c      call set_rule('EVAL_Str',ttype_opme,EVAL,
c     &     'FOPT_Str',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'transformed metric :',0,'LIST')
c      call set_rule('EVAL_Str',ttype_opme,PRINT_MEL,
c     &     'ME_S',1,0,
c     &     parameters,2,tgt_info)
c
c      ! Evaluate metric times vector
c      call add_target('EVAL_S_Ttr',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_S_Ttr','FOPT_S_Ttr',tgt_info)
c      ! (a) set all elements of right hand vector to one
c      call set_rule2('EVAL_S_Ttr',SET_MEL,tgt_info)
c      call set_arg('EVAL_S_Ttr',SET_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_Ttr'/))
c      call set_arg('EVAL_S_Ttr',SET_MEL,'VAL_LIST',2,tgt_info,
c     &             val_rl8=(/1d0,1d0/))
c      call set_arg('EVAL_S_Ttr',SET_MEL,'IDX_LIST',2,tgt_info,
c     &             val_int=(/-1,-2/))
cc dbg
cc      call form_parameters(-1,parameters,2,
cc     &     'Ttr set to :',0,'LIST')
cc      call set_rule('EVAL_S_Ttr',ttype_opme,PRINT_MEL,
cc     &     'ME_Ttr',1,0,
cc     &     parameters,2,tgt_info)
cc dbgend
c      ! (b) evaluate
c      call set_rule2('EVAL_S_Ttr',ASSIGN_ME2OP,tgt_info)
c      call set_arg('EVAL_S_Ttr',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_Dtr'/))
c      call set_arg('EVAL_S_Ttr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'Dtr'/))
c      call set_rule('EVAL_S_Ttr',ttype_opme,EVAL,
c     &     'FOPT_S_Ttr',1,0,
c     &     parameters,0,tgt_info)
c      ! (c) scale back and copy onto vector list
c      call set_dependency('EVAL_S_Ttr','DEF_ME_T',tgt_info)
c      call set_rule2('EVAL_S_Ttr',SCALE_COPY,tgt_info)
c      call set_arg('EVAL_S_Ttr',SCALE_COPY,'LIST_RES',1,tgt_info,
c     &             val_label=(/'ME_Ttr'/)) !tr'/))
c      call set_arg('EVAL_S_Ttr',SCALE_COPY,'LIST_INP',1,tgt_info,
c     &             val_label=(/'ME_OMG'/)) !tr'/))
c      call set_arg('EVAL_S_Ttr',SCALE_COPY,'FAC',2,tgt_info,
c     &             val_rl8=(/1d0,1d0/))
cc dbg
c      call form_parameters(-1,parameters,2,
c     &     'Result of metric times vector :',0,'LIST')
c      call set_rule('EVAL_S_Ttr',ttype_opme,PRINT_MEL,
c     &     'ME_Ttr',1,0, !tr',1,0,
c     &     parameters,2,tgt_info)
cc dbgend
c      ! (d) transformation
c      call set_dependency('EVAL_S_Ttr','DEF_ME_Dtrdag',tgt_info)
c      call set_dependency('EVAL_S_Ttr','FOPT_T',tgt_info)
c      call set_rule2('EVAL_S_Ttr',ASSIGN_ME2OP,tgt_info)
c      call set_arg('EVAL_S_Ttr',ASSIGN_ME2OP,'LIST',1,tgt_info,
c     &           val_label=(/'ME_Dtrdag'/))
c      call set_arg('EVAL_S_Ttr',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
c     &           val_label=(/'Dtr'/))
c      call set_rule('EVAL_S_Ttr',ttype_opme,EVAL,
c     &     'FOPT_T',1,0,
c     &     parameters,0,tgt_info)
cc dbg
c      call form_parameters(-1,parameters,2,
c     &     'After transformation :',0,'LIST')
c      call set_rule('EVAL_S_Ttr',ttype_opme,PRINT_MEL,
c     &     'ME_T',1,0,
c     &     parameters,2,tgt_info)
cc dbgend
c
c      ! (d) free up Ttr list
c      call set_rule('EVAL_S_Ttr',ttype_opme,RES_ME_LIST,
c     &     'ME_Ttr',1,0,
c     &     parameters,0,tgt_info)
c dbgend

c dbg
c      call add_target2('PRINT_Heff',.false.,tgt_info)
c      call set_dependency('PRINT_Heff','FOPT_Heff',tgt_info)
c      call set_dependency('PRINT_Heff','SOLVE_REF',tgt_info)
c      call set_rule2('PRINT_Heff',PRINT_MEL,tgt_info)
c      call set_arg('PRINT_Heff',PRINT_MEL,'LIST',1,tgt_info,
c     &             val_label=(/'ME_T'/))
c      call set_rule('PRINT_Heff',ttype_opme,EVAL,
c     &     'FOPT_Heff',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'effective Hamiltonian :',0,'LIST')
c      call set_rule('PRINT_Heff',ttype_opme,PRINT_MEL,
c     &     'ME_Heff',1,0,
c     &     parameters,2,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'bare Hamiltonian :',0,'LIST')
c      call set_rule('PRINT_Heff',ttype_opme,PRINT_MEL,
c     &     mel_ham,1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      ! Evaluate Residual
c      call add_target('EVAL_OMG',ttype_gen,.false.,tgt_info)
c      call set_dependency('EVAL_OMG','FOPT_OMG',tgt_info)
c      call set_dependency('EVAL_OMG','SOLVE_REF',tgt_info)
cc      call set_dependency('EVAL_OMG','PRINT_Heff',tgt_info)
c      call set_rule('EVAL_OMG',ttype_opme,RES_ME_LIST,
c     &     'ME_OMG',1,0,
c     &     parameters,0,tgt_info)
c      call set_rule('EVAL_OMG',ttype_opme,EVAL,
c     &     'FOPT_OMG',1,0,
c     &     parameters,0,tgt_info)
c      call form_parameters(-1,parameters,2,
c     &     'Residual :',0,'LIST')
c      call set_rule('EVAL_OMG',ttype_opme,PRINT_MEL,
c     &     'ME_OMG',1,0,
c     &     parameters,2,tgt_info)
c dbgend

      ! Evaluate projected T (needed e.g. for restart)
      call add_target('EVAL_Tproj',ttype_gen,.false.,tgt_info)
      call set_dependency('EVAL_Tproj','DEF_ME_Dproj',tgt_info)
      call set_dependency('EVAL_Tproj','FOPT_T',tgt_info)
      call set_dependency('EVAL_Tproj','DEF_ME_T',tgt_info)
      call set_dependency('EVAL_Tproj','DEF_ME_Ttr',tgt_info)
      ! (a) first copy T list to Ttr
      call set_rule2('EVAL_Tproj',SCALE_COPY,tgt_info)
      call set_arg('EVAL_Tproj',SCALE_COPY,'LIST_RES',1,tgt_info,
     &             val_label=(/'ME_Ttr'/))
      call set_arg('EVAL_Tproj',SCALE_COPY,'LIST_INP',1,tgt_info,
     &             val_label=(/'ME_T'/))
      call set_arg('EVAL_Tproj',SCALE_COPY,'FAC',1,tgt_info,
     &             val_rl8=(/1d0/))
      ! (just for safety: should be already assigned this way)
      call set_rule2('EVAL_Tproj',ASSIGN_ME2OP,tgt_info)
      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'LIST',1,tgt_info,
     &           val_label=(/'ME_Dproj'/))
      call set_arg('EVAL_Tproj',ASSIGN_ME2OP,'OPERATOR',1,tgt_info,
     &           val_label=(/'Dtr'/))
      ! (b) evaluate projection
      call set_rule2('EVAL_Tproj',EVAL,tgt_info)
      call set_arg('EVAL_Tproj',EVAL,'FORM',1,tgt_info,
     &             val_label=(/'FOPT_T'/))
c dbg
c      call set_rule('EVAL_Tproj',ttype_opme,PRINT_MEL,
c     &     'ME_Ttr',1,0,
c     &     parameters,2,tgt_info)
c dbgend

c dbg
c      ! check of orthogonality: <C0|T^+ T|C0>
c      call add_target2('F_CTTC',.false.,tgt_info)
c      call set_dependency('F_CTTC','E(MR)',tgt_info)
c      call set_dependency('F_CTTC','T',tgt_info)
c      call set_dependency('F_CTTC','C0',tgt_info)
c      call set_rule2('F_CTTC',EXPAND_OP_PRODUCT,tgt_info)
c      call set_arg('F_CTTC',EXPAND_OP_PRODUCT,'LABEL',1,tgt_info,
c     &     val_label=(/'F_CTTC'/))
c      call set_arg('F_CTTC',EXPAND_OP_PRODUCT,'OP_RES',1,tgt_info,
c     &     val_label=(/'E(MR)'/))
c      call set_arg('F_CTTC',EXPAND_OP_PRODUCT,'OPERATORS',4,
c     &     tgt_info,
c     &     val_label=(/'C0^+','T^+','T','C0'/))
c      call set_arg('F_CTTC',EXPAND_OP_PRODUCT,'IDX_SV',4,tgt_info,
c     &     val_int=(/2,3,4,5/))
c      call set_rule2('F_CTTC',PRINT_FORMULA,tgt_info)
c      call set_arg('F_CTTC',PRINT_FORMULA,'LABEL',1,tgt_info,
c     &     val_label=(/'F_CTTC'/))
c      ! optimized CTTC
c      call add_target2('FOPT_CTTC',.false.,tgt_info)
c      call set_dependency('FOPT_CTTC','F_CTTC',tgt_info)
c      call set_dependency('FOPT_CTTC','DEF_ME_T',tgt_info)
c      call set_dependency('FOPT_CTTC','DEF_ME_C0',tgt_info)
c      call set_dependency('FOPT_CTTC','DEF_ME_E(MR)',tgt_info)
c      call set_rule2('FOPT_CTTC',OPTIMIZE,tgt_info)
c      call set_arg('FOPT_CTTC',OPTIMIZE,'LABEL_OPT',1,tgt_info,
c     &             val_label=(/'FOPT_CTTC'/))
c      call set_arg('FOPT_CTTC',OPTIMIZE,'LABELS_IN',1,tgt_info,
c     &             val_label=(/'F_CTTC'/))
c      ! evaluate CTTC
c      call add_target2('EVAL_CTTC',.true.,tgt_info)
c      call set_dependency('EVAL_CTTC','SOLVE_MRCC',tgt_info)
c      call set_dependency('EVAL_CTTC','FOPT_CTTC',tgt_info)
c      call set_rule2('EVAL_CTTC',EVAL,tgt_info)
c      call set_arg('EVAL_CTTC',EVAL,'FORM',1,tgt_info,
c     &             val_label=(/'FOPT_CTTC'/))
c dbgend

      return
      end
