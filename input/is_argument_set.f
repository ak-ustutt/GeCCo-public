*----------------------------------------------------------------------*
      integer function is_argument_set(context,argkey,keycount)
*----------------------------------------------------------------------*
*     return number of appearences of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input2, only:find_active_node,input_doc,atr_name,
     &     NodeList,
     &     Node,hasChildNodes,getFirstChild,arg_tag,getChildNodes,
     &     getLength,getNodeName,item,getAttribute
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
     &     curkey,input_root
      type(Node), pointer ::
     &     curarg
      type(NodeList),pointer::
     &     nodes_list

      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target,ii

      if (ntest.ge.100) then
          call write_title(lulog,wst_dbg_subr,i_am)
          write(lulog,*) 'looking for argument: ',argkey
          write (lulog,*) 'in context: ',context
      end if
      input_root=>getFirstChild(input_doc)
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount
      
      
      if (ntest.ge.100) then
         write(lulog,*)"looking for the ",icount_target,"th attribute"
      end if 
      call find_active_node(input_root,curkey,
     &     context,icount_target)

      iargcount = 0
      if ( associated(curkey))then
         if (ntest.ge.100) then
            write(lulog,*) getAttribute(curkey,atr_name),
     &           "is associated"
         end if 
         if ( hasChildNodes(curkey)) then
            nodes_list => getChildNodes(curkey) 
            if (ntest.ge.100) then
               write(lulog,*) "nodeList of length:"
     &              ,getLength(nodes_list)
            end if 
            arg_loop: do ii=0,getLength(nodes_list)-1
               curarg=> item(nodes_list,ii)
               if (ntest.ge.100) then
                  write(lulog,*) "current_argument:",
     &                 getAttribute(curarg,atr_name)  
               end if 
               if (getNodeName(curarg) 
     &           .eq. arg_tag .and.
     &           getAttribute(curarg,atr_name)
     &           .eq. trim(argkey) ) iargcount = iargcount+1

            end do arg_loop
         end if 
      end if
      if (ntest.ge.100)then 
         write (lulog,'("found ",i3," occurences")') iargcount
         write (lulog,*) "of ",argkey," under ",context
      end if
      is_argument_set = iargcount

      return
      end
