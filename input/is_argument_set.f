*----------------------------------------------------------------------*
      integer function is_argument_set(context,argkey,keycount)
*----------------------------------------------------------------------*
*     return number of appearences of argument for keyword given 
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
     &     keycount

      type(keyword), pointer ::
     &     curkey
      type(argument), pointer ::
     &     curarg
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target

      if (.not.associated(keyword_history%down_h))
     &     call quit(1,'is_keyword_set','invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount

      call find_active_node(keyword_history,curkey,
     &     context,icount_target)

      iargcount = 0

      if (associated(curkey).and.associated(curkey%arg_h)) then
        curarg => curkey%arg_h

        arg_loop: do 

          if (trim(curarg%key).eq.trim(argkey)) iargcount = iargcount+1

          if (associated(curarg%next)) then
            curarg => curarg%next
          else
            exit arg_loop
          end if        
          
        end do arg_loop

      end if

      is_argument_set = iargcount

      return
      end
