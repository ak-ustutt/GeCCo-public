*----------------------------------------------------------------------*
      integer function is_keyword_set(context)
*----------------------------------------------------------------------*
*     return number of appearences of keyword in currently
*     active keyword_history
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*----------------------------------------------------------------------*

      use parse_input2, only : input_doc,getFirstChild,hasChildNodes,
     &     Node,getNodeName,key_tag, getParentNode,key_root_tag,
     &     keyword_get_context
      use FoX_dom, only : getNextSibling
      implicit none

      integer,parameter::
     &     ntest= 00
      character(len=14),parameter ::
     &     i_am="is_keyword_set"

      character, intent(in) ::
     &     context*(*)


      type(Node), pointer ::
     &     current,input_root

      character ::
     &     curcontext*1024

      integer ::
     &     icount


      input_root=getFirstChild(input_doc)
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')

      current => getFirstChild(input_root)

      icount = 0
      key_loop: do 
         if (getNodeName(current).eq. key_tag)then 
            call keyword_get_context(curcontext,current)
            if (trim(context).eq.trim(curcontext)) icount = icount+1 
         end if 
         
         ! 
         if (hasChildNodes(current)) then
            current => getFirstChild(current)
         else if (associated(getNextSibling(current))) then
            current => getNextSibling(current)
         else
            do while(getNodeName(getParentNode(current))
     &           .ne. key_root_tag .and. 
     &           .not. associated(getNextSibling(current)))
               current => getParentNode(current)
            end do 
            if (associated(getNextSibling(current))) then
               current => getNextSibling(current)
            else
               exit key_loop
            end if
         end if 
      end do key_loop

      is_keyword_set = icount

      return
      end
