
*---------------------------------------------------------------------*
!>     processes the input up to the next calculate block 
!!     sets several common blocks 
*----------------------------------------------------------------------*
      subroutine process_input(one_more,orb_info)
*----------------------------------------------------------------------*

      use keyword_trees, only:inp_postprocess
      use parse_input, only: inp_show,reg_show
      implicit none
      include 'stdunit.h'
      include 'ioparam.h'
      include 'ifc_input.h'
      include 'def_orbinf.h'
      include 'routes.h'

      
      logical, intent(out) ::
     &     one_more
      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     icnt, len, nfreeze, ncnt, ncnt2, nactel, iread
      integer, allocatable ::
     &     iscr(:)
      character ::
     &     str*256
      logical ::
     &     allowed(6)

      one_more=.false.
      call inp_postprocess(one_more)
      if (.not. one_more) return


      ! check input -- start version
      icnt = is_keyword_set('method')

c      if (icnt.le.0) then
c        call quit(0,'process_input','no "method" block specified')
c      end if
      call get_argument_value('general','print',ival=iprlvl)
      write(lulog,*) 'printlevel is set to ',iprlvl

      if (iprlvl.ge.10)
     &   call show_keywords(lulog)

      ! set file block-length
      call get_argument_value('general','da_block',ival=iread)
      lblk_da = iread*1024/nrecfc

      ncnt = is_keyword_set('orb_space.shell')

      allowed(1:6) = .true.
      do icnt = 1, ncnt
        ncnt2 = is_argument_set('orb_space.shell','type',keycount=icnt)
        if (ncnt2.ne.1)
     &       call quit(0,'process_input','single shell? frozen?')
        str(1:256) = ' '
        call get_argument_value('orb_space.shell','type',keycount=icnt,
     &                          str=str)
        select case(trim(str))
        case('frozen')
          if (.not.allowed(1)) cycle
          allowed(1) = .false.
          if (is_argument_set('orb_space.shell','def',
     &                        keycount=icnt).gt.0) then
            call get_argument_dimension(len,'orb_space.shell','def',
     &                                  keycount=icnt)
            allocate(iscr(len))
            call get_argument_value('orb_space.shell','def',
     &                              keycount=icnt,iarr=iscr)
            nfreeze = sum(iscr(1:len))
          else if (is_argument_set('orb_space.shell',
     &                             'nfreeze',keycount=icnt).gt.0) then
            call get_argument_value('orb_space.shell',
     &                             'nfreeze',keycount=icnt,ival=nfreeze)
            len = orb_info%nsym
            allocate(iscr(len))
            call auto_freeze(iscr,nfreeze,orb_info)
          else
            nfreeze = -1
            len = orb_info%nsym
            allocate(iscr(len))
            call auto_freeze(iscr,nfreeze,orb_info)
          end if


          if (nfreeze.gt.0) then 
             call add_frozen_shell(iscr,len,'frz',orb_info)
          end if 
          deallocate(iscr)

        case('occ')
          if (.not.allowed(5)) cycle
          allowed(5) = .false.
          if (.not.allowed(1).or..not.allowed(6))
     &      call quit(0,'process_input',
     &                'type=occ must preceed closed and frozen')

          if (is_argument_set('orb_space.shell','def',
     &                        keycount=icnt).gt.0) then
            call get_argument_dimension(len,'orb_space.shell','def',
     &                                  keycount=icnt)
            allocate(iscr(len))
            call get_argument_value('orb_space.shell','def',
     &                              keycount=icnt,iarr=iscr)
          else 
            call quit(0,'process_input',
     &                'type=occ must be accompanied by a definition')
          end if
          call add_frozen_shell(iscr,len,'occ',orb_info)
          deallocate(iscr)
        case('closed')
          if (.not.allowed(6)) cycle
          allowed(6) = .false.
          if (.not.allowed(1))
     &      call quit(0,'process_input',
     &                'type=closed must preceed frozen')

          if (is_argument_set('orb_space.shell','def',
     &                        keycount=icnt).gt.0) then
            call get_argument_dimension(len,'orb_space.shell','def',
     &                                  keycount=icnt)
            allocate(iscr(len))
            call get_argument_value('orb_space.shell','def',
     &                              keycount=icnt,iarr=iscr)
          else
            call quit(0,'process_input',
     &                'type=occ must be accompanied by a definition')
          end if
          call add_frozen_shell(iscr,len,'cls',orb_info)
          deallocate(iscr)
        case('occorb')
          if (.not.allowed(2)) cycle
          allowed(2) = .false.
cmh       Change of inactive orbitals currently leads to wrong Fock Op.
          call quit(0,'process_input','core Fock op. would be wrong!')
          call get_argument_dimension(len,'orb_space.shell','def',
     &                                keycount=icnt)
          allocate(iscr(len))
          call get_argument_value('orb_space.shell','def',
     &                            keycount=icnt,iarr=iscr)
          call modify_actspc(iscr,len,-1,orb_info,1)
          deallocate(iscr)
        case('actorb')
          if (.not.allowed(3)) cycle
          allowed(3) = .false.
          call get_argument_dimension(len,'orb_space.shell','def',
     &                                keycount=icnt)
          allocate(iscr(len))
          call get_argument_value('orb_space.shell','def',
     &                            keycount=icnt,iarr=iscr)
          call get_argument_value('orb_space.shell','nactel',
     &                            keycount=icnt,ival=nactel)
          call modify_actspc(iscr,len,nactel,orb_info,2)
          deallocate(iscr)
        case('deleted')
          if (.not.allowed(4)) cycle
          allowed(4) = .false.

          if (is_argument_set('orb_space.shell','def',
     &                        keycount=icnt).gt.0) then
            call get_argument_dimension(len,'orb_space.shell','def',
     &                                  keycount=icnt)
            allocate(iscr(len))
            call get_argument_value('orb_space.shell','def',
     &                              keycount=icnt,iarr=iscr)
            nfreeze = sum(iscr(1:len))
          else if (is_argument_set('orb_space.shell',
     &                             'nfreeze',keycount=icnt).gt.0) then
            call get_argument_value('orb_space.shell',
     &                             'nfreeze',keycount=icnt,ival=nfreeze)
            call quit(1,'process_input','no automatic delete yet')
            !len = orb_info%nsym
            !allocate(iscr(len))
            !call auto_freeze(iscr,nfreeze,orb_info)
          else
            call quit(1,'process_input','no automatic delete yet')
            !nfreeze = -1
            !len = orb_info%nsym
            !allocate(iscr(len))
            !call auto_freeze(iscr,nfreeze,orb_info)
          end if

          if (nfreeze.gt.0) 
     &         call add_deleted_shell(iscr,len,orb_info)
          deallocate(iscr)

        case default
          call quit(0,'process_input','unexpected shell type: "'//
     &         trim(str)//'"')
        end select
      end do

      ! set routes for core routines
      call get_argument_value('calculate.routes','schedule',
     &     ival=irt_sched)
      call get_argument_value('calculate.routes','contract',
     &     ival=irt_contr)
      call get_argument_value('calculate.routes','str_block',
     &     ival=len_str_block)
      call get_argument_value('calculate.routes','cnt_block',
     &     ival=len_cnt_block)
      call get_argument_value('calculate.routes','force_batching',
     &     ival=force_batching)
      if (force_batching.gt.2)
     &       call quit(0,'process_input',
     &       'illegal value for force_batching (>2)')          
      call get_argument_value('calculate.routes','force_ooc_sort',
     &     ival=force_ooc_sort)
      call get_argument_value('calculate.routes','maxbranch',
     &     ival=maxbranch)
      call get_argument_value('calculate.routes','auto_opt',
     &     lval=use_auto_opt)
      call get_argument_value('calculate.routes','reo_factor',
     &     xval=reo_factor)
      call get_argument_value('calculate.routes','use_tr',
     &     lval=use_tr)
      call get_argument_value('calculate.routes','sv_fix',
     &     lval=sv_fix)
      call get_argument_value('calculate.routes','sv_old',
     &     lval=sv_old)
      call get_argument_value('calculate.routes','sv_thresh',
     &     xval=sv_thresh)
      call get_argument_value('calculate.routes','jac_fix',
     &     lval=jac_fix)
      call get_argument_value('calculate.routes','jac_thresh',
     &     xval=jac_thresh)
      call get_argument_value('calculate.routes','Tikhonov',
     &     xval=tikhonov)
      call get_argument_value('calculate.routes','spinadapt',
     &     ival=spinadapt)
      if (spinadapt.gt.3)
     &   call quit(0,'process_input','illegal value for spinadapt (>3')
      if (spinadapt.lt.0) then ! set up by default
        ! full spin adaptation for non-singlet MRCC, else none
        spinadapt = 0
        if (is_keyword_set('method.MRCC').gt.0.and.orb_info%imult.ne.1)
     &     spinadapt = 3
      end if
      return

      end subroutine process_input
