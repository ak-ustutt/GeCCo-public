*----------------------------------------------------------------------*
      integer function is_keyword_set(context)
*----------------------------------------------------------------------*
*     return number of appearences of keyword in currently
*     active keyword_history
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*----------------------------------------------------------------------*

      use parse_input2, only : input_doc,getFirstChild,hasChildNodes,
     &     Node,getNodeName,key_tag,getParentNode,key_root_tag,atr_name,
     &     keyword_get_context,getAttribute,key_dsearch
      use FoX_dom, only : getNextSibling
      implicit none
      include 'stdunit.h'

      integer,parameter::
     &     ntest=00
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

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      if (ntest.ge.100)
     &     write (lulog,*) " association, input_doc",
     &     associated(input_doc)
      
      input_root=>getFirstChild(input_doc)
      if (ntest.ge.100)
     &     write (lulog,*) " association, input_root",
     &     associated(input_root)

      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')

      current => getFirstChild(input_root)

      icount = 0
      key_loop: do 
         curcontext=" "
         call keyword_get_context(curcontext,current)
! keyword_get_context only gets the context of the current keyword
         if (len_trim(curcontext).eq.0)then
            curcontext=getAttribute(current,atr_name)
         else
            curcontext=trim(curcontext)//"."//
     &           getAttribute(current,atr_name)
         end if 
         
         if (trim(context).eq.trim(curcontext)) icount = icount+1 
         
                  
         ! 
         current=>key_dsearch(current)
         if (.not.associated(current)) exit key_loop
      end do key_loop

      is_keyword_set = icount
      if (ntest.ge.100)
     &     write (lulog,*) i_am," has found",icount,
     &     "occurences of ",trim(context)
      return
      end
