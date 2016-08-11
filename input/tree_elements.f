!> defines the element classes
!! all are type(node) from FoX_dom
!! conceptual subclasses of element are keywords, arguments and root-elements.
!! Note: the document node is not a root elemen


      module tree_elements
      use FoX_dom
      implicit none

      private
      public :: Node
      public :: elem_getParentNode,elem_getComment,
     &     elem_getFirstChild,elem_getLastChild,
     &     elem_getNextSibling,elem_getPreviousSibling

      public :: getParentNode,getFirstChild,getLastChild,getNextSibling,
     &     getPreviousSibling
      public :: arg_tag,key_tag,key_root_tag
!      integer,parameter::
!     &     ELEMENT_NODE=1

      character(len=*),parameter::
     &     key_tag="keyword",
     &     arg_tag="argument",
     &     key_root_tag="key_root",
     &     comment_tag="#comment"
      contains
*----------------------------------------------------------------------*
!>     tests if a given node is an element
!!     a element is either an keyword, an argument or a root element.
!!   can handle unassociated nodes 
*----------------------------------------------------------------------*
      function is_element(curnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="is_element"

      logical::
     &     is_element
      type(node),pointer,intent(in)::
     &     curnode

      if (associated(curnode))then
         if(getNodeType(curnode).eq.ELEMENT_NODE)then
            select case(getNodeName(curnode))
            case(key_tag)
               is_element=.True.
            case(arg_tag)
               is_element=.True.
            case(key_root_tag)
               is_element=.True.
            case default
               is_element=.False.
            end select
         else
            is_element=.False.
         end if
      else
         is_element=.False.
      end if
      return
      end function

*----------------------------------------------------------------------*
!>    retrieves the first element that is anchestor to the current 
!!     element with tag tag.
!!    returns null if no such element exists 
*----------------------------------------------------------------------*
      function elem_getParentNode(curnode, tag)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_getParent"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode
      character(len=*),intent(in),optional::
     &     tag
      

      nxtnode=>getParentNode(curnode)

      do while (associated(nxtnode))
         if (present(tag) ) then
            !     node name is tag for element, #document for documents, 
            !     #comment for comments #text for text nodes
            if (getNodeName(nxtnode).eq. tag ) return 
         else if (is_element(nxtnode))then
            return
         else 
            continue
         end if

         nxtnode=> getParentNode(nxtnode)
      end do
      end function

*----------------------------------------------------------------------*
!>    retrieves first subnode with specified tag
!!    returns null() if no such element exists
*----------------------------------------------------------------------*
      function elem_getFirstChild(curnode,tag)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest=00
      character(len=*),parameter ::
     &     i_am="elem_getFirstChild"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode
      character(len=*),intent(in),optional::
     &     tag

      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (present(tag)) write(lulog,*) " of type:",trim(tag)
      end if 

      nxtnode=> null()
      if(.not. hasChildNodes(curnode) ) return 
      nxtnode=> getFirstChild(curnode)
      if (present(tag))then
         if (getNodeName(nxtnode).eq. trim(tag)) return
         nxtnode=> elem_getNextSibling(nxtnode,tag)
      else
         if (is_element(nxtnode))return
         nxtnode=> elem_getNextSibling(nxtnode)
      end if
      return
      end function

*----------------------------------------------------------------------*
!>    retrieves next sibling with a specified tag
!!
*----------------------------------------------------------------------*
      function elem_getNextSibling(curnode,tag) result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_getNextSibling"

      type(Node), pointer::
     &     nxtnode
      character(len=*),intent(in),optional::
     &     tag
      type(Node), pointer,intent(in)::
     &     curnode
      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (present(tag)) write(lulog,*) " of type:",trim(tag)
      end if 
      nxtnode=> getNextSibling(curnode)
      do while (associated(nxtnode))
         if (present(tag))then
            if (getNodeName(nxtnode) .eq. trim(tag)) return
            continue
         else
            if (is_element(nxtnode))return
            continue
         end if
         nxtnode=> getNextSibling(nxtnode)
      end do 
      return
      end function

*----------------------------------------------------------------------*
!>    retrieves last direct subnode with specified tag
!!    returns null() if no such element exists
*----------------------------------------------------------------------*
      function elem_getLastChild(curnode,tag)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_getLastChild"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode
      character(len=*),intent(in),optional::
     &     tag

      nxtnode=> null()
      if(.not. hasChildNodes(curnode) ) return 

      nxtnode=> getLastChild(curnode)
      if (present(tag))then
         if (getNodeName(nxtnode).eq. trim(tag)) return
         nxtnode=> elem_getPreviousSibling(nxtnode,tag)
      else
         if (is_element(nxtnode))return
         nxtnode=> elem_getPreviousSibling(nxtnode)
      end if
      return
      end function


*----------------------------------------------------------------------*
!>    retrieves previous sibling with a specified tag
!!
*----------------------------------------------------------------------*
      function elem_getPreviousSibling(curnode,tag) result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_getPreviousSibling"

      type(Node), pointer::
     &     nxtnode
      character(len=*),intent(in),optional::
     &     tag
      type(Node), pointer,intent(in)::
     &     curnode

      nxtnode=> getNextSibling(curnode)

      do while (associated(nxtnode))
         if (present(tag))then
            if (getNodeName(nxtnode) .eq. trim(tag)) return
            continue
         else
            if (is_element(nxtnode))return
            continue
         end if
         nxtnode=> getPreviousSibling(nxtnode)
      end do 
      return
      end function

*----------------------------------------------------------------------*
!>    retrieves the associated comment
!!
!!    returns null if no comment is available
*----------------------------------------------------------------------*
      function elem_getComment(curnode) result(commentnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_getComment"

      type(Node), pointer::
     &     commentnode
      type(Node), pointer,intent(in)::
     &     curnode
      commentnode=>null()
      if(.not. hasChildNodes(curnode)) return
      commentnode=> getFirstChild(curnode)
      do while (associated(commentnode))
         if( getNodeName(commentnode) .eq. comment_tag) return
         commentnode => getNextSibling(commentnode)
      end do
      return
      end function
      end module 
