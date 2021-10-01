*----------------------------------------------------------------------*
      integer function is_keyword_set(context)
*----------------------------------------------------------------------*
*     return number of appearences of keyword in currently
*     active keyword_history
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      character, intent(in) ::
     &     context*(*)

      type(keyword), pointer ::
     &     current
      character ::
     &     curcontext*1024
      integer ::
     &     icount

      if (.not.associated(keyword_history%down_h))
     &     call quit(1,'is_keyword_set','invalid keyword history')
      current => keyword_history%down_h

      icount = 0
      key_loop: do 
        if (current%status.gt.0) then
          call keyword_get_context(curcontext,current)
          if (trim(context).eq.trim(curcontext)) icount = icount+1 
        end if

        if (current%status.gt.0.and.associated(current%down_h)) then
          ! go down
          current => current%down_h
        else if (associated(current%next)) then
          ! else stay within level
          current => current%next
        else
          ! else find an upper level, where a next
          ! node exists:
          up_loop: do
            if (associated(current%up)) then
              current => current%up
              if (associated(current%next)) then
                current => current%next
                exit up_loop
              end if
            else
              exit key_loop
            end if
          end do up_loop
        end if        

      end do key_loop

      is_keyword_set = icount

      return
      end
