*----------------------------------------------------------------------*
      subroutine set_mr_targets(tgt_info,orb_info,env_type,
     &     name_infile,fforbinf)
*----------------------------------------------------------------------*
*     calls target generators for multireference methods
*
*     matthias, march 2011
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest =  00

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(inout) ::
     &     orb_info
      character(len=*), intent(in) ::
     &     env_type
      type(filinf), intent(in) ::
     &     fforbinf
      character(*), intent(in) ::
     &     name_infile

      integer ::
     &     maxexc, cmaxexc, maxh, maxp, mult, ms, sym, maxtop, maxcum
      logical ::
     &     l_icci, l_iccc, use_met, fixed, use_f12,response,multistate,
     &     set_up_T_corr,skip
      integer, allocatable ::
     &     excrestr(:,:,:)
      real(8) ::
     &     densmix

      integer ::
     &     stndD(2,60), nsupD, nremblk, remblk(60), len
      character(len=256) ::
     &     gecco_path
      

      skip = (is_keyword_set('calculate.skip_E').gt.0)
      ! redefine spin/symmetry of reference state (if requested)
      call get_argument_value('method.MR','mult',
     &     ival=mult)
      call get_argument_value('method.MR','ms',
     &     ival=ms)
      call get_argument_value('method.MR','sym',
     &     ival=sym)
      call get_argument_value('method.MR','densmix',
     &     xval=densmix)
      if (mult.gt.0.and.mult.ne.orb_info%imult) then
        orb_info%imult = mult
        if (ntest.ge.100) write(lulog,*) 'spin mult. = ', mult
      end if
      if (ms.le.orb_info%imult.and.ms.ne.orb_info%ims) then
        orb_info%ims = ms
        if (ntest.ge.100) write(lulog,*) '2Ms        = ', ms
      end if
      ! Ms possible?
      if (orb_info%ims.lt.1-orb_info%imult
     &    .or.orb_info%ims.gt.orb_info%imult-1
     &    .or.mod(orb_info%imult-orb_info%ims,2).eq.0)
     &   call quit(1,'set_mr_targets','impossible Ms')
      if (sym.gt.0.and.sym.ne.orb_info%lsym) then
        orb_info%lsym = sym
        if (ntest.ge.100) write(lulog,*) 'symmetry   = ', sym
        if (sym.gt.orb_info%nsym) call quit(1,'set_mr_targets',
     &           'impossible symmetry')
      end if

      ! print orbital informations (after get it changed from input)

      call put_orbinfo(orb_info, fforbinf)

      call get_argument_value('method.MR','multistate',
     &     lval=multistate)

      call get_argument_value('method.MRCC','set_up_T_corr',
     &     lval=set_up_T_corr)

      ! get maximum excitation rank
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      call get_argument_value('method.MR','cmaxexc',
     &     ival=cmaxexc)

      call get_environment_variable( "GECCO_DIR", value=gecco_path,
     &     length = len)

      if (len.EQ.0)
     &     call quit(1,'set_mr_targets',
     &     "Please, set the GECCO_DIR environment variable.")

      ! icMRCI calculation?
      l_icci = is_keyword_set('method.MRCI').gt.0
      ! icMRCC calculation?
      l_iccc = is_keyword_set('method.MRCC').gt.0
      ! both at same time is not allowed
      if (l_icci.and.l_iccc) call quit(1,'set_mr_targets',
     &   'Now don''t be greedy, choose either icMRCI or icMRCC!')

      if ((l_icci.or.l_iccc).and.cmaxexc.gt.0)
     &            call quit(1,'set_mr_targets',
     &            'Warning: Only tested for CASSCF reference so far')

      ! set R12 intermediates, if necessary
      if (is_keyword_set('method.R12').gt.0) then
        use_f12 = .true.
        call get_argument_value('method.R12','fixed',lval=fixed)
        if (.not.fixed) 
     &      call warn('set_mr_targets', 
     &           'MR + R12: fixed amplitudes are the only choice!')
        call set_r12f_general_targets(tgt_info,orb_info,env_type)
      else
        use_f12 = .false.
      end if
      ! set response targets, if necessary
      if (is_keyword_set('calculate.excitation').gt.0) then
        response=.true.
      else
        response=.false.
      endif

      ! first set targets for CASSCF or uncontracted CI wave function
      call set_unc_mrci_targets(tgt_info,orb_info,
     &                          .not.((l_icci.or.l_iccc).or.skip),
     &                          name_infile,fforbinf%name)
c dbg for calculating cumulants
c      call set_gno_targets(tgt_info,orb_info,1)
c dbgend

      ! if maxexc = 0: return because call of unc_mrci is sufficient
      if (maxexc.eq.0.or..not.(l_icci.or.l_iccc)) then
        call get_argument_value('method.MR','maxcum',
     &     ival=maxcum)
        if (maxcum.gt.0) call set_gno_targets(tgt_info,orb_info,0)
        return
      end if

      ! set targets associated with generalized normal order
      ! (includes reduced density matrices)
      maxtop = 1
      if (l_iccc) then
        call get_argument_value('method.MRCC','maxcom_res',
     &       ival=maxtop)
      end if
      call set_gno_targets(tgt_info,orb_info,maxtop)

      ! get restrictions on excitation classes
      call get_argument_value('method.MR','maxh',
     &     ival=maxh)
      call get_argument_value('method.MR','maxp',
     &     ival=maxp)
      if (max(maxh,maxp).gt.maxexc)
     &     call quit(1,'set_mr_targets',
     &               'maxh and maxp must not exceed maxexc')
      if (maxh.lt.0) maxh = maxexc
      if (maxp.lt.0) maxp = maxexc
      ! additional restrictions from a small number of occupied
      ! or virtual orbitals?
      maxh = min(maxh,2*orb_info%nactt_hpv(IHOLE))
      maxp = min(maxp,2*orb_info%nactt_hpv(IPART))
      allocate(excrestr(0:maxh,0:maxp,1:2))
      call get_exc_restr(excrestr,maxh,maxp,
     &                   orb_info%nactel,orb_info%nactorb)

      ! set targets common for internally contracted methods
      call set_ic_mr_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp,use_met,
     &                       nsupD,stndD,nremblk,remblk)

      ! set icMRCI or icMRCC targets
      if (l_icci) call set_ic_mrci_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp,use_met)
      if (l_iccc) call set_ic_mrcc_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp,.not.use_f12,
     &                       nsupD,stndD,nremblk,remblk,
     &                       name_infile,fforbinf%name)
      if (use_f12) call set_ic_mrcc_f12_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp)
      if (response) call set_python_targets(tgt_info,
     &     trim(gecco_path)//"/python_spec/icmrcc_ee_targets.py",
     &     name_infile,fforbinf%name)
      deallocate(excrestr)

      if (multistate) call set_python_targets(tgt_info,
     &     trim(gecco_path)//"/python_spec/multistate_eff_ham.py",
     &     name_infile,fforbinf%name)

      if (densmix.gt.0d0.and.maxexc.lt.3) 
     &     call set_python_targets(tgt_info,
     &     trim(gecco_path)//"/python_spec/printe.py",
     &     name_infile,fforbinf%name)

      if (set_up_T_corr) call set_python_targets(tgt_info,
     &     trim(gecco_path)//"/python_spec/set_up_T_corr.py",
     &     name_infile,fforbinf%name)



      return
      end
