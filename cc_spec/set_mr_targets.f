*----------------------------------------------------------------------*
      subroutine set_mr_targets(tgt_info,orb_info)
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

      integer, parameter ::
     &     ntest = 100

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     maxexc, cmaxexc, maxh, maxp, mult, ms, sym, maxtop, maxcum
      logical ::
     &     l_icci, l_iccc, use_met
      integer, allocatable ::
     &     excrestr(:,:,:)

      ! redefine spin/symmetry of reference state (if requested)
      call get_argument_value('method.MR','mult',
     &     ival=mult)
      call get_argument_value('method.MR','ms',
     &     ival=ms)
      call get_argument_value('method.MR','sym',
     &     ival=sym)
      if (mult.gt.0.and.mult.ne.orb_info%imult) then
        orb_info%imult = mult
        if (ntest.ge.100) write(luout,*) 'spin mult. = ', mult
      end if
      if (sym.gt.0.and.sym.ne.orb_info%lsym) then
        orb_info%lsym = sym
        if (ntest.ge.100) write(luout,*) 'symmetry   = ', sym
        if (sym.gt.orb_info%nsym) call quit(1,'set_mr_targets',
     &           'impossible symmetry')
      end if
      orb_info%ims = ms
      if (ntest.ge.100) write(luout,*) 'Ms         = ', ms
      ! Ms possible?
      if (ms.lt.1-orb_info%imult.or.ms.gt.orb_info%imult-1
     &    .or.mod(orb_info%imult-ms,2).eq.0)
     &   call quit(1,'set_mr_targets','impossible Ms')

      ! get maximum excitation rank
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      call get_argument_value('method.MR','cmaxexc',
     &     ival=cmaxexc)

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

      ! first set targets for CASSCF or uncontracted CI wave function
      call set_unc_mrci_targets(tgt_info,orb_info,
     &                          .not.(l_icci.or.l_iccc))
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
      allocate(excrestr(0:maxh,0:maxp,1:2))
      call get_exc_restr(excrestr,maxh,maxp,
     &                   orb_info%nactel,orb_info%nactorb)

      ! set targets common for internally contracted methods
      call set_ic_mr_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp,use_met)

      ! set icMRCI or icMRCC targets
      if (l_icci) call set_ic_mrci_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp,use_met)
      if (l_iccc) call set_ic_mrcc_targets(tgt_info,orb_info,
     &                       excrestr,maxh,maxp)
      deallocate(excrestr)

      return
      end
