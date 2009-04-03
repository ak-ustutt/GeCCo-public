*----------------------------------------------------------------------*
      subroutine process_input(one_more,orb_info)
*----------------------------------------------------------------------*
*     post-process input up to next "calculate" block
*----------------------------------------------------------------------*

      use parse_input

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

      type(keyword), pointer ::
     &     current

      integer ::
     &     icnt, len, nfreeze
      integer, allocatable ::
     &     iscr(:)
      character ::
     &     str*256

      if (.not.associated(history_pointer)) then
        ! advance to first keyword
        history_pointer => keyword_history
        if (.not.associated(history_pointer%down_h)) then
          call quit(0,'process_input','not a single keyword given?')
        end if
        history_pointer => history_pointer%down_h
      else if (associated(history_pointer%next)) then
        ! next keyword
        history_pointer => history_pointer%next
      else
        one_more = .false.
        return
      end if
      one_more = .true.

      current => history_pointer 

      ! advance until next "calculate" block
      ! set present block "active" and all preceeding
      ! blocks with same label "inactive"
      do
        call set_keyword_status(current,+1)
        call unset_previous_keyword(current)

        if (trim(current%key).eq.'calculate') exit

        if (associated(current%next)) then
          current => current%next
        else
          exit
        end if

      end do

      ! save history pointer
      history_pointer => current

      if (iprlvl.ge.10)
     &   call keyword_list(luout,keyword_history)

      ! check input -- start version
      icnt = is_keyword_set('method')

c      if (icnt.le.0) then
c        call quit(0,'process_input','no "method" block specified')
c      end if

      ! set file block-length
      call get_argument_value('general','da_block',ival=lblk_da)
      lblk_da = lblk_da*1024/nrecfc

      icnt = is_keyword_set('orb_space.shell')

      if (icnt.eq.1) then
        icnt = is_argument_set('orb_space.shell','type')
        if (icnt.ne.1)
     &       call quit(0,'process_input','single shell? frozen?')
        call get_argument_value('orb_space.shell','type',str=str)

        if (str(1:6).ne.'frozen')
     &       call quit(0,'process_input','unexpected shell type: "'//
     &         trim(str)//'"')

        if (is_argument_set('orb_space.shell','def').gt.0) then
          call get_argument_dimension(len,'orb_space.shell','def')
          allocate(iscr(len))
          call get_argument_value('orb_space.shell','def',iarr=iscr)
        else if (is_argument_set('orb_space.shell','nfreeze').gt.0) then
          call get_argument_value('orb_space.shell',
     &                            'nfreeze',ival=nfreeze)
          len = orb_info%nsym
          allocate(iscr(len))
          call auto_freeze(iscr,nfreeze,orb_info)
        else
          nfreeze = -1
          len = orb_info%nsym
          allocate(iscr(len))
          call auto_freeze(iscr,nfreeze,orb_info)
        end if

        call add_frozen_shell(iscr,len,orb_info)
        deallocate(iscr)

      end if

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
      call get_argument_value('calculate.routes','use_tr',
     &     lval=use_tr)

      return

      contains

*----------------------------------------------------------------------*
      subroutine set_keyword_status(keywd,status)
*----------------------------------------------------------------------*
*     set status of current keyword and all sub-levels
*----------------------------------------------------------------------*

      implicit none

      type(keyword), pointer ::
     &     keywd
      integer, intent(in) ::
     &     status

      type(keyword), pointer ::
     &     current
      integer ::
     &     level

      current => keywd
      level = 0

      key_loop: do
        current%status = status
        if (associated(current%down_h)) then
          ! go down
          current => current%down_h
          level = level+1
        else if (level.gt.0.and.associated(current%next)) then
          ! go next
          current => current%next
        else if (level.gt.0) then
          ! else find an upper level, where a next
          ! node exists:
          up_loop: do
            if (level.gt.0.and.associated(current%up)) then
              current => current%up
              level = level-1
              if (level.gt.0.and.associated(current%next)) then
                current => current%next
                exit up_loop
              end if
            else
              exit key_loop
            end if
          end do up_loop
        else
          exit key_loop
        end if

      end do key_loop

      return
      end subroutine set_keyword_status

*----------------------------------------------------------------------*
      subroutine unset_previous_keyword(keywd)
*----------------------------------------------------------------------*
*     set all previous keywords with same key as present keyword to
*     inactive (+ all sub-levels)
*----------------------------------------------------------------------*

      implicit none

      type(keyword), intent(in) ::
     &     keywd

      type(keyword), pointer ::
     &     current
      character ::
     &     key*(lenkey)

      ! no previous keyword -> return
      if (.not.associated(keywd%prev)) return

      current => keywd%prev
      key = trim(keywd%key)

      rev_loop: do
        if (trim(current%key).eq.key)
     &     call set_keyword_status(current,-1)
        if (associated(current%prev)) then
          current => current%prev
        else
          exit rev_loop
        end if
      end do rev_loop

      end subroutine unset_previous_keyword

      end subroutine process_input
