*----------------------------------------------------------------------*
      integer function is_argument_set(context,argkey,keycount)
*----------------------------------------------------------------------*
*     return number of appearences of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use keyword_trees, only:Node,inp_arg_from_context
      implicit none
      include "stdunit.h"

      integer,parameter::
     &     ntest= 00
      character(len=15),parameter ::
     &     i_am="is_argument_set"

      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount

      type(Node), pointer ::
     &     curarg

      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target,ii

      if (ntest.ge.100) then
          call write_title(lulog,wst_dbg_subr,i_am)
          write(lulog,*) 'looking for argument: ',argkey
          write (lulog,*) 'in context: ',context
      end if

      icount_target = 1
      if (present(keycount)) icount_target = keycount
      


      if (ntest.ge.100) then
         write(lulog,*)"looking for the ",icount_target,"th attribute"
      end if 
      
      iargcount=-1
      curarg=> inp_arg_from_context(context,argkey,latest=.false.,
     &     keycount=icount_target,argcount=iargcount)

      
      is_argument_set = -iargcount-1

      if (ntest.ge.100)then 
         write (lulog,'("found ",i3," occurences")') is_argument_set
         write (lulog,*) "of ",argkey," under ",context
      end if
      
      return
      end
