*----------------------------------------------------------------------*
      integer function is_argument_set(context,argkey,keycount)
*----------------------------------------------------------------------*
*     return number of appearences of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input2, only:find_node,input_doc,atr_name,NodeList,
     &     Node,hasChildNodes,getFirstChild,arg_tag,getChildNodes,
     &     getLength,getNodeName,item,getAttribute
      implicit none

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
      
      input_root=getFirstChild(input_doc)
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount

      call find_node(input_root,curkey,
     &     context)

      iargcount = 0
      if ( associated(curkey))then
         if ( hasChildNodes(curkey)) then
            nodes_list => getChildNodes(curkey) 

            arg_loop: do ii=1,getLength(nodes_list)  

               if (getNodeName(item(nodes_list,ii)) 
     &           .eq. arg_tag .and.
     &           getAttribute(item(nodes_list,ii),atr_name)
     &           .eq. trim(argkey) ) iargcount = iargcount+1

            end do arg_loop
         end if 
      end if

      is_argument_set = iargcount

      return
      end
