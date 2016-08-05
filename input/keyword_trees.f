      module keyword_trees
      use FoX_dom
!,only:DOMException,NodeList,
!     &     getAttribute,hasAttribute,setAttribute,
!     &     getNodeName,
!     &     getElementsByTagName,getLength,item,
!     &     parseString,parseFile,
!     &     getExceptionCode,getOwnerDocument,createElement,appendChild
      use tree_elements
      implicit none


      private

      public :: Node
      
      public :: arg_tag,key_tag
      public :: atr_name,atr_kind,atr_len,atr_val

! printing function needs knowledge of status
      public :: atr_stat, status_active, status_inactive

!     attribute accessors
      public :: getAttribute,hasAttribute,setAttribute

      public :: elem_getParentNode,elem_getComment,
     &     elem_getFirstChild,elem_getLastChild,
     &     elem_getNextSibling,elem_getPreviousSibling

      public :: getNodeName
      public :: tree_t
      public :: getRoot,getLevel
      
      public :: create_keyword_trees

      public :: fetch_input_keyword_tree,fetch_registry_keyword_tree

      public :: tree_iterate,tree_goback_to_element
      
      public :: tree_get_arg_from_context,tree_get_key_from_context

      public :: getSubNode

      public :: tree_create_new_element
      public :: max_name_len
      type tree_t
        type(node),pointer :: root=> null()
        integer  :: curlevel= 0
        type(node),pointer :: curnode => null()
      end type 
      
      character,parameter ::
     &     context_sep="."
      
      character,parameter::     !attribute names
     &     atr_name*4="name",
     &     atr_kind*4="type",
     &     atr_len*6="length",
     &     atr_val*5="value",
     &     atr_stat*6="status"


      character,parameter::
     &     status_active="A",
     &     status_inactive="I"

      integer,parameter::
     &     max_context_lvl=16,    ! 16 level deep keyword context should be enough
     &     max_context_len=256

      integer,parameter::
     &     max_name_len=16





      type(node),pointer :: 
     &     registry_doc,        !root of the tree that represents the keyword registry
     &     input_doc            !document element of input

      contains 

!=====================================================================
! accessors for the module variables.
!=====================================================================
*----------------------------------------------------------------------*
!>    returns the keyword_root element of the registry
!!
*----------------------------------------------------------------------*
      function fetch_registry_keyword_tree() result(keytree)
*----------------------------------------------------------------------*
      implicit none 
      include "stdunit.h"
      character(len=*),parameter ::
     &     i_am="fetch_registry_keyword_tree"
      integer, parameter ::
     &     ntest=00
      type(tree_t)::
     &     keytree
      type(NodeList),pointer::
     &     root_list

      if (.not. associated(registry_doc)) call quit(1,i_am,
     &     "no registry tree found")
      root_list => getElementsByTagName(registry_doc, key_root_tag)

      if(getLength(root_list).ne.1)call quit(1,i_am,
     &     "found more/less than one root element")

      keytree%root=> item(root_list, 0) ! 0-indexed ... that's no true Fortran
      keytree%curlevel=0
      keytree%curnode=> keytree%root
      end function



*----------------------------------------------------------------------*
!>    returns the keyword_root element of the input
!!
*----------------------------------------------------------------------*
      function fetch_input_keyword_tree() result(keytree)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="fetch_input_keyword_tree"
      integer, parameter ::
     &     ntest=00
      type(tree_t)::
     &     keytree

      if (.not. associated(input_doc)) call quit(1,i_am,
     &     "no input tree found")

      keytree%root=> getFirstChild(input_doc)
      
      if (.not. associated(keytree%root)) call quit(1,i_am,
     &     "no root element found for input")

      keytree%curlevel=0
      keytree%curnode=> keytree%root

      end function

!=======================================================================
! initialization subroutines
!=======================================================================

*----------------------------------------------------------------------*
!!    parses the keyword_file and initializes the input
!!  
!!    @param file name of the input file
*----------------------------------------------------------------------*
      subroutine create_keyword_trees(file)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=*),intent(in)::
     &     file
      character(len=13),parameter ::
     &     i_am="create_keyword_trees"
      integer, parameter ::
     &     ntest=00
      type(DOMException)::
     &     ex
      type(tree_t) ::
     &     registry
      type(node),pointer ::
     &     reg_root,dummy
      if (ntest.gt.100)then 
         call write_title(lulog,wst_dbg_subr,i_am)
      end if 

      call reg_init_()
      call inp_start_()

      write(lulog,*) "reading keyword_file",trim(file) 
      registry_doc => parseFile(trim(file), ex=ex)
      !! @TODO more special error handling 
      if (getExceptionCode(ex) .ne. 0)
     &     call quit(1,i_am,"could not read keyword registry")


      registry = fetch_registry_keyword_tree() 
      reg_root=> getRoot(registry)
      call setAttribute(reg_root,atr_name,"registry")
      end subroutine
*----------------------------------------------------------------------*
!!    initializes the module variables relating to the registry
*----------------------------------------------------------------------*
      subroutine reg_init_()
*----------------------------------------------------------------------*
      implicit none
      registry_doc =>null()
      end subroutine

*----------------------------------------------------------------------*
!!    initializes the module variables relating to the input
*----------------------------------------------------------------------*
      subroutine inp_init_()
*----------------------------------------------------------------------*
      implicit none
      input_doc=> null()
      end subroutine

*----------------------------------------------------------------------*
!!    creates the document to hold the input variables
*----------------------------------------------------------------------*
      subroutine inp_start_()
*----------------------------------------------------------------------*
      implicit none
      type(Node),pointer::
     &     input_root
      call inp_init_()
      
      !easiest way to create a document
      input_doc=> parseString(
     &     "<"//key_root_tag//" "//atr_name//'="input">'
     &     //"</"//key_root_tag//">")
      end subroutine inp_start_


!======================================================================!
!     routinen for keyword transplantation
!======================================================================!
*----------------------------------------------------------------------*
!>     creates an element of the same type and at the same position 
!!     as template on tree 
!!     
!!     a context corresponding to template's context has to exist on tree.
!!     @param 
!!     @param template registry element that provides the template
*----------------------------------------------------------------------*
      function tree_create_new_element( tree, template) 
     &     result(new_elem)
*----------------------------------------------------------------------*
      implicit none

      include "stdunit.h"
      integer,parameter::
     &     ntest=00
      character(len=*),parameter::
     &     i_am = "tree_create_new_element"
      type(tree_t),intent(inout)::
     &     tree
      type(Node),pointer::
     &     new_elem
      type(Node),pointer,intent(in)::
     &     template
      character(len=max_context_len)::
     &     context
      
      type(Node),pointer::
     &     doc
      type(Node),pointer::
     &     new_parent ! not a new parent, but newly a parent

      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*)" template Name:",getAttribute(template,atr_name)
      end if 
      call  key_get_context(template,context) !using context to synchronize inp and registry

      if (ntest.ge.100) then
         write (lulog,*) "creating a new keyword/argument on context:",
     &   context
      end if

      doc=> getOwnerDocument(tree%root)
      new_elem=>createElement(doc,getNodeName(template))
      new_parent=> tree_get_key_from_context(tree,trim(context),
     &     latest=.True.)

      if (ntest.ge.100) then
         write (lulog,*) "tag, name of parent:",
     &        getNodeName(new_parent),getAttribute(new_parent,atr_name)
      end if
      
      ! Oh, an adoption so cuuute
      new_elem=>appendChild(new_parent, new_elem)

      return
      end function




*----------------------------------------------------------------------*
!>    returns the context of a given keyword
!!    @param[out] curcontext the returned context
!!    @param[in] keywd keywd which context is to be described
!!    @param[in] full optional logical if true, the context includes the keywords name
*----------------------------------------------------------------------*
      subroutine key_get_context(keywd,curcontext,full)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=19),parameter::
     &     i_am="keyword_get_context"
      integer,parameter::
     &     ntest=00
      integer, parameter ::
     &     maxlen  = 256
      character(len=*),intent(inout) ::
     &     curcontext
      type(node), pointer,intent(in)::
     &     keywd 
      logical,intent(in),optional::
     &     full

      character(len=max_name_len)::
     &     curname
      
      type(node), pointer ::
     &     internal

      curcontext=" "
      if (present(full))then
         if (full) curcontext=getAttribute(keywd,atr_name)//context_sep
      end if 
      internal=> getParentNode(keywd)
      do while (getNodeName(internal) .ne. key_root_tag)
         curcontext=getAttribute(internal,atr_name)//
     &   context_sep//trim(curcontext)
         internal=> getParentNode(internal)
      end do
      curcontext(len_trim(curcontext):len_trim(curcontext))=" "
      end subroutine





!=======================================================================
! basic accessors
!=======================================================================
*----------------------------------------------------------------------*
!>     get the root element of a tree
*----------------------------------------------------------------------*
      function getRoot(tree)
*----------------------------------------------------------------------*
      type(tree_t),intent(in)::
     &     tree
      type(node),pointer::
     &     getRoot

      getRoot=> tree%root
      end function
*----------------------------------------------------------------------*
!>    how far below the root element is the current element?
*----------------------------------------------------------------------*
      function getLevel(tree)
      type(tree_t),intent(in)::
     &     tree
      integer::
     &     getLevel

      getLevel= tree%curlevel
      end function



!=======================================================================
!iterate over the tree
!=======================================================================

*----------------------------------------------------------------------*
!>    general iterative depth first walk
!!
!!    @param[in] curkey pointer to parent keyword
!!    @param[inout] level 
!!    @return next element
*----------------------------------------------------------------------*
      function tree_iterate(tree) result(nxtnode) 
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="tree_iterate"

      type(tree_t),intent(inout)::
     &     tree
      type(node),pointer::
     &     nxtnode
      type(node),pointer::
     &     sibling, child, parent
      integer::
     &     level

      nxtnode=>tree %curnode
      level =tree%curlevel

      sibling=> elem_getNextSibling(nxtnode)
      child=> elem_getFirstChild(nxtnode)


      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if 

      if (associated(child))then
         nxtnode=> child
         level=level+1
      else if (associated(sibling))then
         nxtnode => sibling
      else
         up_loop: do 
            nxtnode => elem_getParentNode(nxtnode)  ! makes problems if there are non element nodes in the tree.
            level=level-1 ! 
            if (level .eq. 0) then
               
            end if
            if (.not.associated(nxtnode)) then 
               return
            end if 
                  
            sibling=>elem_getNextSibling(nxtnode)
            if (associated(sibling))then 
               nxtnode=> sibling
               exit up_loop
            end if 
         end do up_loop
      end if

      

      if (ntest .ge. 100 .and. associated(nxtnode)) then
         write (lulog,*) " found node",getNodeName(nxtnode)
      else if (ntest.ge.100)then
         write (lulog,*) " no node found"
      end if
      tree%curnode=> nxtnode
      tree%curlevel=level
      end function

!======================================================================
!backtracking search
!======================================================================
*----------------------------------------------------------------------*
!!     finds an element that is subelement to any keyword in the current context
!!     sets its context to current keyword
!!     also returns the element
*----------------------------------------------------------------------*
      function tree_goback_to_element(tree,name,tag) 
     &     result(nxtnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="tree_goback_to_element"
      integer,parameter::
     &     ntest=00

      type(node),pointer ::
     &     parnode,nxtnode
      type(tree_t),intent(inout)::
     &     tree
      character(len=*),intent(in) ::
     &     name
      character(len=*),intent(in) ::
     &     tag
      integer ::
     &     icount
      nxtnode=> null()
      parnode=> tree%curnode

      tree%curlevel=tree%curlevel+1
      do while (tree%curlevel .ge.0)
         icount=1
         nxtnode=> getSubnode(parnode,name,tag=tag,latest=.False. , 
     &        icount=icount)
         if (.not.associated(nxtnode))then 
            parnode=> elem_getParentNode(parnode)
            tree%curlevel=tree%curlevel-1
         else
            if (tag .eq. arg_tag)then
               tree%curnode=> parnode
            else
               tree%curnode=> nxtnode
            end if
            exit
         end if
      end do 
      

      end function

!======================================================================
!guided search
!======================================================================


*----------------------------------------------------------------------*
!>    resolves a context and name to an argument in a tree
!!
!!    does not change current element pointer of that tree
!!    @param[in] context context of the Argument
!!    @param[in] name name of the Argument
!!    @param[in] latest should the Arguments be searched in reverse order
!!    @param[inout] keycount, index of the keyword in identical context that will be accessed
!!    @param[inout] argcount  index of the argument under the keyword 
*----------------------------------------------------------------------*
      function tree_get_arg_from_context(tree,context,name,
     &     latest,keycount,argcount)  result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="tree_get_arg_from_context"

      type(node), pointer ::
     &     finnode
      type(tree_t),intent(in) ::
     &     tree
      character(len=*),intent(in) ::
     &     context,name
      logical,intent(in)::
     &     latest
      integer,intent(inout),optional::
     &     keycount,argcount
      integer::
     &     ikeycount,iargcount
      if (.not.present(keycount)) then
         ikeycount=1
      else
         ikeycount=keycount
      end if

      if (.not.present(argcount)) then
         iargcount=1
      else
         iargcount=argcount
      end if

      finnode=> arg_from_context(context,name,tree%root,latest,
     &     ikeycount,iargcount)

      if(present(argcount))then
         argcount=iargcount
      end if
      end function

*----------------------------------------------------------------------*
!>    resolves a context to a keyword in a tree
!!
!!    does not change current element pointer of that tree
!!    @param[in] context context of the Keyword
!!    @param[in] latest should the Keyword be searched in reverse order
!!    @param[inout] keycount, index of the keyword in identical context that will be accessed
*----------------------------------------------------------------------*
      function tree_get_key_from_context(tree,context,
     &     latest,keycount)  result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="arg_from_context"

      type(node), pointer ::
     &     finnode
      type(tree_t),intent(in) ::
     &     tree
      character(len=*),intent(in) ::
     &     context
      logical,intent(in),optional::
     &     latest
      integer,intent(inout),optional::
     &     keycount
      integer::
     &     ikeycount
      logical::
     &     ilatest
      if (.not.present(keycount)) then
         ikeycount=1
      else
         ikeycount=keycount
      end if

      if(.not. present(latest))then
         ilatest=.False.
      else
         ilatest=latest
      end if

      finnode=> key_from_context(context,tree%root,latest,
     &     ikeycount)

      if (present(keycount))keycount=ikeycount

      end function
*----------------------------------------------------------------------*
!>    resolves a context and name to an argument
!!    @param[in] context context of the Argument
!!    @param[in] name name of the Argument
!!    @param[in] latest should the Arguments be searched in reverse order
!!    @param[inout] keycount, index of the keyword in identical context that will be accessed
!!    @param[inout] argcount  index of the argument under the keyword 
!!      keycount and argcount are 0 if they were positive and an argument was found
*----------------------------------------------------------------------*
      function arg_from_context(context,name,tree_root,latest,keycount,
     &     argcount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      include "stdunit.h"
      character(len=*),parameter ::
     &     i_am="arg_from_context"
      integer, parameter ::
     &     ntest = 00

      type(node), pointer ::
     &     finnode
      character(len=*),intent(in) ::
     &     context,name
      type(node), pointer,intent(in) ::
     &     tree_root
      logical,intent(in)::
     &     latest
      integer,intent(inout)::
     &     keycount,argcount

      type(node), pointer::
     &     curkey

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' argument = "',name,'"'
         write(lulog,*) ' context = "',trim(context),'"'
        if (.not.latest) write(lulog,*) ' forward search'
        if (latest) write(lulog,*) ' backward search'
      end if 
      finnode=>null()
      curkey=>key_from_context(context,tree_root,latest,keycount)
      
      if (associated(curkey))then 
         finnode=>getSubNode(curkey,name,arg_tag,latest,argcount)
      end if

      if(ntest.ge.100.and. associated(finnode))
     &     write (lulog,*) "success"
      if(ntest.ge.100.and. .not. associated(finnode))
     &     write (lulog,*) "No success"

      end function














*----------------------------------------------------------------------*
!>    resolves a context to a keyword
!!
!!    
!!    @param context of the keyword including that keywords name
!!    @param[in] tree_root starting point of search
!!    @param[in] latest for reversed search direction
!!    Note: only traverses not Inactive keywords
*----------------------------------------------------------------------*
      function key_from_context(context,tree_root,latest,keycount) 
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      include "stdunit.h"
!!   traverses the tree with a stack reliant depth first search
!!
!!   for the context only the beginning of every name is saved (ipst)
!!   for the elements: on every level, the index of last checked node is saved (sikeykcount). 
!!      and increase upon revisiting this level.
!!   as getSubnode reduces its icount parameter, sikeycount has to be given as ikeycount to this function
!!   if the last level is reached keycount is reduced by the amount of found occurences if it reaches 0 this is the returned key.
      character(len=*),parameter ::
     &     i_am="key_from_context"

      integer, parameter ::
     &     ntest = 100
      integer,parameter ::
     &     max_stack=max_context_lvl

      type(node), pointer ::
     &     finnode
      character(len=*),intent(in) ::
     &     context
      type(node), pointer,intent(in) ::
     &     tree_root
      logical,intent(in)::
     &     latest
      integer,intent(inout)::
     &     keycount
      type(node), pointer ::
     &     curnode

      integer::
     &     ipst, ipnd, len,
     &     ikeycount, sikeycount            !keycount in this sublevel and saved keycount
      integer ::
     &     st_stack(max_stack), ist_stack,   !stack for starting indices 
     &     keyc_stack(max_stack), ikeyc_stack,  ! stack for sublevel keycounts
     &     ierr
      logical ::
     &     ctxt_end
      
      ierr=0
      
      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' context = "',trim(context),'"'
        if (.not.latest) write(lulog,*) ' forward search'
        if (latest) write(lulog,*) ' backward search'
      end if 
      
      call stack_init(st_stack,ist_stack)
      call stack_init(keyc_stack,ikeyc_stack)
    
      curnode=>null()
      finnode => tree_root
      len = len_trim(context)
      if (len.eq.0) return !nothing to do

      ipst = 1
      ctxt_end=.False.       !flag if the current name is the last in the context
      sikeycount=1


      do
         ! finding next slice in context
         ipnd = index(context(ipst:),context_sep)+ipst-2
         if (ipnd.lt.ipst) ctxt_end=.True.
         if(ctxt_end) ipnd=len


         if (ntest.ge.100) then
            write(lulog,*) ' current subkeyword: "',
     &           context(ipst:ipnd),'" idx:',ipst,ipnd
            write(lulog,*) ' in final level',ctxt_end
            write(lulog,*) ' keys to be found:',keycount
         end if
 
         ! comparing to Nodes
         if (ctxt_end)then
            ! key_getSubkey looks for the keycount'th subnode
            !  AND reduces keycount by the amount of found nodes!!!
            ! key_tag ensures that we only look for keywords
            curnode=> getSubnode(finnode,context(ipst:ipnd),key_tag,
     &           latest,keycount)

            if (ntest.ge.100) then
               if (associated(curnode))write(lulog,*) 
     &              ' looking at "',getAttribute(curnode,atr_name),'"'
               write(lulog,*) ' in final level',ctxt_end
               write(lulog,*) ' keys still to be found:',keycount
            end if

            if (keycount.eq.0)then
               finnode=>curnode
               exit
            end if
         else            
            ikeycount=sikeycount
            curnode=> getSubnode(finnode, context(ipst:ipnd), key_tag,
     &           latest, ikeycount)
         end if


         
         ! traversing tree
         if (associated(curnode))then 
            ! go one level down
            if (ntest .gt. 100) then
               write(lulog,*) "found",getAttribute(curnode,atr_name)
               write(lulog,*) "going one level down"
               write(lulog,*) "saving:",ipst,sikeycount
            end if 

            call stack_push(st_stack,ist_stack,ipst,ierr)
            call stack_push(keyc_stack,ikeyc_stack,sikeycount,ierr)
            if (ierr.gt.0) call quit (1,i_am,"stack exceeded")
      
            finnode=>curnode
            ipst=ipnd+2
            sikeycount=1
         else if(ipst.eq.1)then
            ! go one level up: can't
            finnode=> null()
            exit
         else
            ! go one level up
            call stack_pop(st_stack,ist_stack,ipst,ierr)
            call stack_pop(keyc_stack,ikeyc_stack,sikeycount,ierr)
            if (ierr.gt.0) call quit(1,i_am,
     &           " trying to access negative stack. Why, just why?")
            ctxt_end=.False.
            finnode=> elem_getParentNode(finnode)
            sikeycount=sikeycount+1

            if (ntest .gt. 100) then
               write(lulog,*) "going one level up"
               write(lulog,*) "starting at:",ipst,sikeycount
            end if 
         end if 
      end do

      call stack_del(st_stack,ist_stack)
      call stack_del(keyc_stack,ikeyc_stack)
      
      if(ntest.ge.100.and. associated(finnode))
     &     write (lulog,*) "success"
      if(ntest.ge.100.and. .not. associated(finnode))
     &     write (lulog,*) "No success"
      return 

      contains
*----------------------------------------------------------------------*
!>    implementing a simple stack for our starting indices
!!
!!    need the external parameter max_stack
!!    @param curkey pointer to current keyword
!!    @param[inout] level on in old level on out new level
*----------------------------------------------------------------------*
      subroutine stack_init(stack,stack_pointer)
      implicit none
      integer,intent(inout)::
     &     stack(max_stack),stack_pointer
      stack_pointer=0
!     stack =0 not neccessary
      end subroutine 

*----------------------------------------------------------------------*
      subroutine stack_push(stack,stack_pointer,i,ierr)
      implicit none
      integer,intent(inout)::
     &     stack(max_stack),stack_pointer
      integer,intent(in)::
     &     i
      integer, intent(inout)::
     &     ierr
  
      stack_pointer=stack_pointer+1
      if (stack_pointer.gt.max_stack) then 
         ierr=ierr+1
         return
      end if 
      stack(stack_pointer)=i
      end subroutine

*----------------------------------------------------------------------*
      subroutine stack_pop(stack,stack_pointer,i,ierr)
      implicit none
      integer,intent(inout)::
     &     stack(max_stack),stack_pointer
      integer,intent(out)::
     &     i
      integer, intent(out)::
     &     ierr
  
      if (stack_pointer.le.0) then 
         ierr=ierr+1
         return
      end if 
      i=stack(stack_pointer)
      stack_pointer=stack_pointer-1
      end subroutine

*----------------------------------------------------------------------*
      subroutine stack_del(stack,stack_pointer)
      integer,intent(inout)::
     &     stack(max_stack),stack_pointer
      continue
      return
      end subroutine 
      end function





*----------------------------------------------------------------------*
!>    retrieves an element below curkey
!!
!!    @param curkey pointer to parent keyword
!!    @param name name of the argument
!!    @param tag tag of the element
!!    @param latest if ocurrences should be counted backwards
!!    @param icount the icountth occurrence is retrieved
*----------------------------------------------------------------------*
      function getSubNode(curkey,name,tag,latest,icount)
     &     result(nextnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=*),parameter::
     &     i_am="getSubNode"
      integer,parameter::
     &     ntest=00

      type(Node),pointer::
     &     nextnode

      type(Node),pointer,intent(in)::
     &     curkey
      character(len=*),intent(in)::
     &     name ,tag
      logical,intent(in)::
     &     latest
      integer,intent(inout)::
     &     icount

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) " looking for:",trim(name)
         write(lulog,*) " of type:",trim(tag)
         if (latest)write(lulog,*) " reversed order"
         write(lulog,*)  " number to be found:",icount
      end if 
   
      if (latest)then
         nextnode=>elem_getLastChild(curkey, tag=tag)
      else 
         nextnode=>elem_getFirstChild(curkey,tag=tag)
      end if 

      do 
         if (.not. associated(nextnode)) exit
      if (ntest .gt. 100) then
         write(lulog,*)  " next node:",getAttribute(nextnode,atr_name)
         write(lulog,*)  " type:",getNodeName(nextnode)
         write(lulog,*)  " status:",getAttribute(nextnode,atr_stat)
         write(lulog,*)  " number to be found:",icount
      end if 

         if ( ( getAttribute(nextnode,atr_name) .eq. trim(name) ).and.
     &        ( getAttribute(nextnode,atr_stat) .ne. status_inactive))
     &        then
               icount=icount-1
            if (icount.eq.0) then
!     nextnode points correctly
               if (ntest .gt. 100)  write(lulog,*)"success"
               return
            end if
         end if 
         if (latest)then
            nextnode=>elem_getPreviousSibling(nextnode,tag=tag)
         else 
            nextnode=>elem_getNextSibling(nextnode,tag=tag)
         end if
      end do 
      return 
      end function

      end module
