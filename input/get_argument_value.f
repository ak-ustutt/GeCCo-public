*----------------------------------------------------------------------*
      subroutine get_argument_value(
     &     context,argkey,keycount,argcount,
     &     ival,iarr,lval,larr,xval,xarr,str)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      logical, intent(out), optional ::
     &     lval, larr(*)
      integer, intent(out), optional ::
     &     ival, iarr(*)
      real(8), intent(out), optional ::
     &     xval, xarr(*)
      character, intent(out), optional ::
     &     str*(*)

      type(keyword), pointer ::
     &     curkey
      type(argument), pointer ::
     &     curarg
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target, type, dim, idx
      logical ::
     &     try_default, succ

      if (.not.associated(keyword_history%down_h))
     &     call quit(1,'get_argument_value','invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount

      try_default = .false.
      call find_active_node(keyword_history,curkey,
     &     context,icount_target)

      iargcount = 0
      iargcount_target = 1
      if (present(argcount)) iargcount_target = argcount

      succ = .false.

      ! try default instead
      if ((.not.associated(curkey).or..not.associated(curkey%arg_h))
     &     .and.iargcount_target.eq.1) then
        try_default = .true.
        call find_node(keyword_root,curkey,context)
      end if

      if (associated(curkey).and.associated(curkey%arg_h)) then
        curarg => curkey%arg_h

        arg_loop: do 

          if (trim(curarg%key).eq.trim(argkey)) iargcount = iargcount+1
          if (trim(curarg%key).eq.trim(argkey).and.
     &         iargcount.eq.iargcount_target) then
            dim = curarg%val%len
            type = curarg%val%type
            select case(type)
            case (vtyp_log)
              if (.not.(present(lval).or.present(larr)))
     &             call quit(1,'get_argument_value',
     &             'no l-value array present')
              if (allocated(curarg%val%lval)) then
                if (present(lval)) lval = curarg%val%lval(1)
                if (present(larr)) larr(1:dim) = curarg%val%lval(1:dim)
                succ = .true.
              end if
            case (vtyp_int)
              if (.not.(present(ival).or.present(iarr)))
     &             call quit(1,'get_argument_value',
     &             'no i-value array present')
              if (allocated(curarg%val%ival)) then
                if (present(ival)) ival = curarg%val%ival(1)
                if (present(iarr)) iarr(1:dim) = curarg%val%ival(1:dim)
                succ = .true.
              end if
            case (vtyp_rl8)
              if (.not.(present(xval).or.present(xarr)))
     &             call quit(1,'get_argument_value',
     &             'no r-value array present')
              if (allocated(curarg%val%xval)) then
                if (present(xval)) xval = curarg%val%xval(1)
                if (present(xarr)) xarr(1:dim) = curarg%val%xval(1:dim)
                succ = .true.
              end if
            case (vtyp_str)
              if (.not.(present(str)))
     &             call quit(1,'get_argument_value',
     &             'no r-value array present')
              if (allocated(curarg%val%cval)) then
                do idx = 1, dim
                  str(idx:idx) = curarg%val%cval(idx)
                end do
                succ = .true.
              end if
            end select
          end if

          if (succ) exit arg_loop
  
          if (associated(curarg%next)) then
            ! go to next argument
            curarg => curarg%next
          else
            if (try_default) exit arg_loop
            ! else try default (if applicable)
            try_default = .true.
            call find_node(keyword_root,curkey,context)

            if (.not.associated(curkey).or.
     &           .not.associated(curkey%arg_h)) exit arg_loop

            curarg => curkey%arg_h
          end if        
          
        end do arg_loop

      end if

      if (.not.succ)
     &     call quit(1,'get_argument_value',
     &     'Could not provide any value for '//trim(context)//
     &     '.'//trim(argkey))

      return
      end
