*----------------------------------------------------------------------*
      subroutine get_argument_dimension(dim,
     &     context,argkey,keycount,argcount,type)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input
      implicit none

      integer, intent(out) ::
     &     dim
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      integer, intent(out), optional ::
     &     type

      type(keyword), pointer ::
     &     curkey
      type(argument), pointer ::
     &     curarg
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target

      if (.not.associated(keyword_history%down_h))
     &     call quit(1,'argument_dimension','invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount

      call find_active_node(keyword_history,curkey,
     &     context,icount_target)

      dim = -1
      if (present(type)) type = -1
      iargcount = 0
      iargcount_target = 1
      if (present(argcount)) iargcount_target = argcount

      if (associated(curkey).and.associated(curkey%arg_h)) then
        curarg => curkey%arg_h

        arg_loop: do 

          if (trim(curarg%key).eq.trim(argkey)) iargcount = iargcount+1
          if (iargcount.eq.iargcount_target) then
            dim = curarg%val%len
            if (present(type)) type = curarg%val%type
            exit arg_loop
          end if

          if (associated(curarg%next)) then
            curarg => curarg%next
          else
            exit arg_loop
          end if        
          
        end do arg_loop

      end if

      return
      end
