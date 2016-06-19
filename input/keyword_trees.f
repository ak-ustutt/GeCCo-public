!> This module works mostly on the following two objects:
!! The keyword DOM-tree  and the input DOM-tree
!! both have the following  general structure 
!!   * A document: (According to W3C specification every xml tree lives upon a document)
!!   * An element with the key_root_tag as tag;
!!               -# in the input tree this is the only child of the document node
!!               -# in the registry tree it can be a lower child. it has to be unique
!!   * Elements with the key or argument tag which are  Subnodes to the key_root element
!!   * Keywords have a name attribute and possibly keywords and /or arguments as children
!!   * Arguments have the following attributes name, len, kind
!!          -# atr_name contains the name
!!          -# atr_len contains the length
!!          -# atr_kind contains the number of the kind of argument(see 'par_vtypes.h')
!!          -# atr_active contains the active/inactive status
!!          -# note that all attributes are saved as strings. they can be converted via rts


!! General Idea is that keyword_trees holds the information storage and parse_input requests containers(Nodes) to stor e informations in it. 
!! this separation is not ~completely~ implemented


!!    
!!
      module keyword_trees
      use FoX_dom
      implicit none

      private
      !Things to export :
      ! information container
      public :: Node
      ! some constants (optimally these would not be exported, but this way it is easier)
      public :: arg_tag,key_tag
      ! constants for atr_names (could be used only in parse_input and the get_arg routines)
      public :: atr_name,atr_kind,atr_len,atr_val
      
      ! printing function needs knowledge of status
      public :: atr_stat, status_active, status_inactive

      ! attribute accessors
      public :: getAttribute,hasAttribute,setAttribute


      !some low level functions for the parse_input search function
      public :: getNodeName,key_root_tag , getParentNode,getFirstChild

      public :: filtered_dsearch


      public :: key_getFirstSubkey,key_getFirstArgument
     &     ,iterate_siblingargs
      public :: key_getArgument,key_getSubkey
      ! gives parse_input a start node so every keyword has only to be searched relative to the last.
      public :: inp_fetch_root,reg_fetch_root

      ! navigation routines for the is_keyword and is_argument set versions
      public :: inp_key_from_context,reg_key_from_context
      public :: inp_arg_from_context,reg_arg_from_context

      ! entrance routines for some actions 
      public :: inp_create_new_element
      public :: inp_postprocess,reg_import
   


      include 'stdunit.h'
      include 'par_vtypes.h'

      !! @TODO make this usable on systems where / is not directory separator
      !! @TODO maybe make this configurable in the installation process
      character,parameter::     !attribute names
     &     atr_name*4="name",
     &     atr_kind*4="type",
     &     atr_len*6="length",
     &     atr_val*5="value",
     &     atr_stat*6="status"

      integer,parameter::
     &     name_len=16

      character,parameter::
     &     status_active="A",
     &     status_inactive="I"

      integer,parameter::
     &     max_context_lvl=16,    ! 16 level deep keyword history should be enough
     &     max_context_len=256
      character,parameter::
     &     context_sep="."

      character,parameter::      !tags 
     &     key_root_tag*8="key_root",
     &     key_tag*7="keyword",
     &     arg_tag*8="argument"


      type(Node), pointer :: 
     &     registry_doc,         !root of the tree that represents the keyword registry
     &     input_doc,           !document element of input
     &     history_pointer      ! to  this "calculate" block the calculation has advanced

      

      contains 

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
! Subroutines coupled to the module variables
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
!>    returns the keyword_root element of the registry
!!
*----------------------------------------------------------------------*
      function reg_fetch_root() result(root)
*----------------------------------------------------------------------*
      implicit none 
      character(len=14),parameter ::
     &     i_am="reg_fetch_root"
      integer, parameter ::
     &     ntest=00
      type(Node),pointer::
     &     root
      type(NodeList),pointer::
     &     root_list
      
      if (.not. associated(registry_doc)) call quit(1,i_am,
     &     "no registry tree found")
      root_list => getElementsByTagName(registry_doc, key_root_tag)

      if(getLength(root_list).ne.1)call quit(1,i_am,
     &     "found more/less than one root element")

      root=> item(root_list, 0) ! 0-indexed ... that's no true Fortran

      end function


*----------------------------------------------------------------------*
!>    returns the keyword_root element of the input
!!
*----------------------------------------------------------------------*
      function inp_fetch_root() result(root)
*----------------------------------------------------------------------*
      implicit none 
      character(len=14),parameter ::
     &     i_am="reg_fetch_root"
      integer, parameter ::
     &     ntest=00
      type(Node),pointer::
     &     root

      if (.not. associated(input_doc)) call quit(1,i_am,
     &     "no input tree found")
      root=> getFirstChild(input_doc)
      
      if (.not. associated(root)) call quit(1,i_am,
     &     "no root element found for input")

      end function


*----------------------------------------------------------------------*
!>    returns the keyword_root element of the input
!!
*----------------------------------------------------------------------*
      function inp_fetch_doc() result(ret)
*----------------------------------------------------------------------*
      implicit none 
      character(len=13),parameter ::
     &     i_am="reg_fetch_doc"
      integer, parameter ::
     &     ntest=00
      type(Node),pointer::
     &     ret
      if (.not. associated(input_doc)) call quit(1,i_am,
     &     "no input tree found")

      ret=> input_doc
 
      if (.not. associated(ret)) call quit(1,i_am,
     &     "no root element found for input")

      end function





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
      history_pointer=> null()
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
     &     "<"//key_root_tag//"></"//key_root_tag//">")
      input_root=> inp_fetch_root()
      call setAttribute(input_root,atr_name,"input")
      end subroutine inp_start_

*----------------------------------------------------------------------*
!!    parses the keyword_file and initializes the input
!!  
!!    @param file name of the input file
*----------------------------------------------------------------------*
      subroutine reg_import(file)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=*),intent(in)::
     &     file
      character(len=13),parameter ::
     &     i_am="reg_import"
      integer, parameter ::
     &     ntest=00
      type(DOMException)::
     &     ex
      type(Node),pointer ::
     &     reg_root
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


      reg_root=> reg_fetch_root() 
      
      call setAttribute(reg_root,atr_name,"registry")
      end subroutine

*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!  output subroutines
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*




*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!  entrance for input postprocessing
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*

*----------------------------------------------------------------------*
!>   wrapper to uncouple set_input_status_ from module variables
!!
!!  @param  one_more logical to show if there are unprocessed blocks remaining.
!!          input status not used
*----------------------------------------------------------------------*
      subroutine inp_postprocess(one_more)
*----------------------------------------------------------------------*
      implicit none
      logical, intent(inout)::
     &     one_more
      type(Node),pointer::
     &     input_root, history_pointer

      input_root=>  inp_fetch_root()
      call tree_set_input_status_(input_root, history_pointer,one_more)
      end subroutine 
*----------------------------------------------------------------------*
!>     creates an element of the same type and at the same position 
!!     as template on input
!!     
!!     @param template registry element that provides the template
*----------------------------------------------------------------------*
      function inp_create_new_element(template) result(new_elem)
*----------------------------------------------------------------------*
      implicit none

      integer,parameter::
     &     ntest=00
      character(len=*),parameter::
     &     i_am = "create_empty_element"

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

      doc=> inp_fetch_doc()
      new_elem=>createElement(doc,getNodeName(template))

      new_parent=> inp_key_from_context(trim(context),latest=.True.)
      if (ntest.ge.100) then
         write (lulog,*) "tag, name of parent:",
     &        getNodeName(new_parent),getAttribute(new_parent,atr_name)
      end if
      
      ! Oh, an adoption so cuuute
      new_elem=>appendChild(new_parent, new_elem)
      call setAttribute(new_elem,atr_stat,status_active)

      return
      end function









*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
! self contained subroutines
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!   postprocessing subroutinen
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*

*----------------------------------------------------------------------*
!>     set all previous keywords with same key as present keyword to
!!     inactive (+ all sub-levels)
*----------------------------------------------------------------------*
      subroutine key_unset_previous_keywords(keywd)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=23),parameter ::
     &     i_am="unset_previous_keywords"
      integer, parameter ::
     &     ntest=00
      type(Node), pointer,intent(in) ::
     &     keywd     
      type(Node), pointer ::
     &     current
      character(len=name_len) ::
     &     name
      integer :: i
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of keyword:"
     &        ,associated(keywd)
      end if 

      current=>getPreviousSibling(keywd)
      if(.not.associated(current)) return 
      name= trim(getAttribute(keywd,atr_name))

      if (ntest.ge.100) then
         write (lulog,*) "saved name:",trim(name).eq.  
     &        trim(getAttribute(keywd,atr_name))
      end if 

      do 
         if (ntest.ge.100) then
            write (lulog,*) "current name:"
     &           ,getAttribute(current,atr_name)
         end if 
         
        if (trim(getAttribute(current,atr_name)).eq.trim(name))
     &        call tree_set_status(current,status_inactive)
        current =>getPreviousSibling(current)
        if (associated(current) )then
           continue
        else
          exit 
        end if
      end do 
      end subroutine

*----------------------------------------------------------------------*
!!    toggles the status of all keywords set in the input:
!!
!!    advances to next 'calculate' Block 
!!   for any given toplevel context only the last block is set active
!!   all other keywords are set inactive.
!!   @param in_root root element of the input DOM-tree
!!   @param history pointer to the last active history file (may be null())
!!   @param one_more true if active blocks were found
*----------------------------------------------------------------------*
      subroutine tree_set_input_status_(in_root, history,one_more)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="set_input_status"
      type(Node),pointer,intent(in)::
     &     in_root
      type(Node),pointer,intent(inout)::
     &     history
      logical, intent(inout)::
     &     one_more
      type(Node),pointer::
     &     current,nxtkey
      integer::
     &     i
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (associated(history))
     &        write (lulog,*) " starting at block: ",
     &        getAttribute(history,atr_name)
              write (lulog,*) " root: ",getAttribute(in_root,atr_name)
      end if 
      one_more=.true.
      if (.not.associated(history)) then
! as the relevant subroutines only work on active elements, all elements were by default active
! in the input creation
         call tree_set_status(in_root,status_inactive)
!      if (ntest.ge.100) then
!         call show_input(lulog)
!      end if 

         history =>getFirstChild( in_root)
         if (.not.associated(  history)) then
            call quit(0,i_am,'not a single keyword given?')
         end if
         if(ntest.ge.100)
     &        write (lulog,*) "history_pointer initiated:",
     &           getAttribute(history,atr_name)
      else 
         current=> getNextSibling(history)
         if (associated( current ) )then
            history =>current
            if (associated(history).and. ntest.ge.100)
     &           write (lulog,*) " advancing history to block: ",
     &           getAttribute(history,atr_name)
            
         else
            if(ntest.ge.100)
     &           write (lulog,*) " no more history: "
            one_more = .false.
            return
         end if
      end if 

      current => history

      do 
         call tree_set_status(current,status_active)
         call key_unset_previous_keywords(current)

         if (trim(getAttribute(current,atr_name)).eq.'calculate') exit

         current=> getNextSibling( current)
         if (associated(  current ) ) then
            history => current
            if (ntest.ge.100)
     &           write (lulog,*) " advancing to block: ",
     &           getAttribute(history,atr_name)
         else
            exit
         end if
      end do
      if (ntest.ge.100) then
         call show_input(lulog)
      end if 
      end subroutine 

*----------------------------------------------------------------------*
!>    sets the given keyword and all subkeywords to a given status
!!
!!   @param keywd pointer to a keyword
!!   @param status integer with the status we set it to 
*----------------------------------------------------------------------*
      subroutine tree_set_status(keywd,status)
*----------------------------------------------------------------------*
      implicit none
      type(Node), pointer, intent(in)::
     &     keywd
      character,intent(in)::
     &     status
      type(Node),pointer::
     &     current 
      integer::
     &     level
      level=0
      current => keywd

      call setAttribute(current,atr_stat,status)
      ! goto sublevel if not going to sublevel, do nothing
      current=> filtered_dsearch(current,level,only_keys=.True.
     &        , key_root_stop=.True.)

      do while(associated(current).and.level.gt.0)
         call setAttribute(current,atr_stat,status)
         current=> filtered_dsearch(current,level,only_keys=.True.
     &        , key_root_stop=.True.)
      end do
      end subroutine











*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
!     Some general use subroutines
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
!>    returns the context of a given keyword
!!    @param[out] curcontext the returned context
!!    @param[in] keywd keywd which context is to be described
!!    @param[in] full optional logical if true, the context includes the keywords name
*----------------------------------------------------------------------*
      subroutine key_get_context(keywd,curcontext,full)
*----------------------------------------------------------------------*
      implicit none
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

      character(len=name_len)::
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














*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!     Iterating search functions
!     if called repeatedly these functions will iterate over the remaining tree in a depths first search
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
      function getParentKey(curnode)  result(nxtkey)
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="getParentKey"

      type(Node), pointer::
     &     nxtkey
      type(Node), pointer,intent(in)::
     &     curnode
      
      nxtkey=> getParentNode(curnode)
      if (getNodeName(nxtkey).eq. key_tag) return
      nxtkey => null()

      end function
      
*----------------------------------------------------------------------*
!>  retrieves first subnode that is a keyword 
!!  returns null() if no such element exists
*----------------------------------------------------------------------*
      function key_getFirstSubkey(curkey)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblingkeys"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curkey

      nxtnode => elem_getFirstSubnode(curkey, key_tag)

      end function
*----------------------------------------------------------------------*
!>  retrieves first subnode that is a keyword 
!!  returns null() if no such element exists
*----------------------------------------------------------------------*
      function key_getFirstArgument(curkey)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblingkeys"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curkey

      nxtnode => elem_getFirstSubnode(curkey, arg_tag)

      end function

*----------------------------------------------------------------------*
!>    retrieves first subnode with specified tag
!!    returns null() if no such element exists
*----------------------------------------------------------------------*
      function elem_getFirstSubnode(curnode,tag)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblingkeys"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode
      character(len=*),intent(in)::
     &     tag

      nxtnode=> null()
      if(.not. hasChildNodes(curnode) ) return 

      nxtnode=> getFirstChild(curnode)
      if (getNodeName(nxtnode).eq. trim(tag)) return

      nxtnode=> iterate_siblings(nxtnode,tag)

      end function


*----------------------------------------------------------------------*
!>    iterates over all siblings that are keys
!!    returns null() if no such key exists 
*----------------------------------------------------------------------*
      function iterate_siblingkeys(curnode)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblingkeys"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode

      nxtnode=> iterate_siblings(curnode,key_tag)
      end function 

*----------------------------------------------------------------------*
!>    iterates over all siblings that are arguments
!!    returns null() if no such argument exists 
*----------------------------------------------------------------------*
      function iterate_siblingargs(curnode)  result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblingargs"

      type(Node), pointer::
     &     nxtnode
      type(Node), pointer,intent(in)::
     &     curnode

      nxtnode=> iterate_siblings(curnode,arg_tag)
      end function 


*----------------------------------------------------------------------*
!>    iterates over all siblings, returns only siblings with specified tag
*----------------------------------------------------------------------*
      function iterate_siblings(curnode,tag) result(nxtnode)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="iterate_siblings"

      type(Node), pointer::
     &     nxtnode
      character(len=*),intent(in)::
     &     tag
      type(Node), pointer,intent(in)::
     &     curnode

      nxtnode=> getNextSibling(curnode)
      do while (associated(nxtnode))
         if (getNodeName(nxtnode) .eq. trim(tag)) exit 
         nxtnode=> getNextSibling(nxtnode)
      end do 
      end function
*----------------------------------------------------------------------*
!>    keyword_specific version of the dsearch function
!!
!!
!!    @param curkey pointer to current keyword
!!    @param[inout] level on in old level on out new level
*----------------------------------------------------------------------*
      function filtered_dsearch(curkey, level, 
     &     only_keys, key_root_stop) 
     &     result(nxtkey) 
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="filtered_dsearch"

      type(Node), pointer::
     &     nxtkey

      type(Node), pointer,intent(in)::
     &     curkey
      integer,intent(inout),optional ::
     &     level
      logical,intent(in),optional ::
     &     only_keys, key_root_stop

      integer::
     &     dlevel !levelchange
    
      logical ::
     &     ionly_keys, ikey_root_stop, show_keys, show_args
      
      ionly_keys=.True.
      if(present(only_keys))ionly_keys=only_keys
      ikey_root_stop=.True.
      if(present(key_root_stop))ikey_root_stop=key_root_stop

      show_keys=.True.
      show_args=.not. ionly_keys

      
      dlevel=0
      nxtkey=>curkey

      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) "starting at:",
     &        getAttribute(nxtkey,atr_name)
         if (present(level))write(lulog,*) "level:", level
      end if 


      main_loop: do
         nxtkey=>elem_dsearch(nxtkey,dlevel)
         if(associated(nxtkey).and. ntest .ge.100) then 
            if (getNodeName(nxtkey) .eq. arg_tag .or. 
     &           getNodeName(nxtkey) .eq. key_tag) then
               write(lulog,*)" returning:",getAttribute(nxtkey,atr_name)
            else 
               write(lulog,*)" returning:", getNodeName(nxtkey)
            end if 
         else if (ntest.ge.100) then
            write(lulog,*) "not associated"
         end if


         if (.not.associated(nxtkey))exit main_loop

         if ( ikey_root_stop .and.
     &        getNodeName(nxtkey).eq.trim(key_root_tag) ) then
            nxtkey=>null()
            exit main_loop
         else if  (show_keys .and.
     &           getNodeName(nxtkey).eq. trim(key_tag))then
            exit main_loop 
         else  if (show_args .and. 
     &           getNodeName(nxtkey).eq. trim(key_tag))then
            exit main_loop
         end if
      end do main_loop
      if (present(level))level=level+dlevel



      if(associated(nxtkey).and. ntest .ge.100) then
         if (getNodeName(nxtkey) .eq. arg_tag .or. 
     &        getNodeName(nxtkey) .eq. key_tag) then
         write(lulog,*) " returning:", getAttribute(nxtkey,atr_name)
         else 
            write(lulog,*) " returning:", getNodeName(nxtkey)
         end if
         if (present(level))write(lulog,*) "exitlevel:", level
      end if 
      end function







*----------------------------------------------------------------------*
!>    general iterative depth first walk
!!
!!    @param[in] curkey pointer to parent keyword
!!    @param[inout] level 
!!    @return next element
*----------------------------------------------------------------------*
      function elem_dsearch(curkey,level) result(nxtkey)
*----------------------------------------------------------------------*
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="elem_dsearch"
      type(Node),pointer::
     &     nxtkey

      type(Node),pointer,intent(in)::
     &     curkey
      integer, intent(inout)::
     &     level
      
      type(Node),pointer::
     &     sibling


      sibling=>getNextSibling(curkey)
      nxtkey=>curkey

      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if 

      if (hasChildNodes(nxtkey))then
         nxtkey=> getFirstChild(nxtkey)
         level=level+1
      else if (associated(sibling))then
         nxtkey => sibling
      else
         up_loop: do 
            nxtkey => getParentNode(nxtkey)
            level=level-1

            if (.not.associated(nxtkey)) then !  
               return
            end if 
                  
            sibling=>getNextSibling(nxtkey)
            if (associated(sibling))then 
               nxtkey=> sibling
               exit up_loop
            end if 
         end do up_loop
      end if
      if (ntest .ge. 100 .and. associated(nxtkey)) then
         write (lulog,*) " found node",getNodeName(nxtkey)
      else if (ntest.ge.100)then
         write (lulog,*) " no node found"
      end if 
      end function


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!     guided search functions 
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
! external subroutines should only be concerned with  (sadly they aren't)
!     reg_key_from_context(context,latest,keycount)
!     inp_key_from_context(context,latest,keycount)
!     reg_arg_from_context(context,name,latest,keycount,argcount)
!     inp_arg_from_context(context,name,latest,keycount,argcount)
! These functions only see active keywords


*----------------------------------------------------------------------*
!>    resolves an argument(context and name) in the registry
!!    
!!    @param context string of the form <key>.<key>.<key>
!!    @param name name of the argument
!!    @param latest optional if true the tree is searched in reverse order
!!    @param keycount looking under the keycount'th appearance of the keyword !active only 
!!    @param argcount looking for the argcount'th appearance under the specified keyword
!!    @return pointer to the argument node or null() if no node was found
!!    NOTE, if keycount or argcount are set to an negative value, 
!!    the return will be null() and the attribute will be reduced by the number of found nodes
!!    its not a bug, it's a feature
*----------------------------------------------------------------------*
      function reg_arg_from_context(context,name,latest,keycount,
     &     argcount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="reg_key_from_context"
      integer, parameter ::
     &     ntest = 00

      
      type(node), pointer ::
     &     finnode
      character(len=*),intent(in) ::
     &     context,name
      logical,intent(in),optional::
     &     latest
      integer,intent(in),optional::
     &     argcount,keycount

      logical::
     &     reversed
      integer::
     &     iargcount,ikeycount

      iargcount=1
      if(present(argcount)) iargcount=argcount

      ikeycount=1
      if(present(keycount)) ikeycount=keycount

      reversed=.false.
      if (present(latest)) reversed = latest
      
      finnode=>reg_fetch_root()
      finnode=> arg_from_context(context,name,finnode,reversed,
     &     ikeycount,iargcount)

      end function

*----------------------------------------------------------------------*
!>    resolves an argument (context and name) in the input
!! 
!!    see reg_arg_from_context   
*----------------------------------------------------------------------*
      function inp_arg_from_context(context,name,latest,keycount,
     &     argcount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="inp_arg_from_context"
      integer, parameter ::
     &     ntest = 00

      
      type(node), pointer ::
     &     finnode
      character(len=*),intent(in) ::
     &     context,name
      logical,intent(in),optional::
     &     latest
      integer,intent(inout),optional::
     &     argcount,keycount
      logical::
     &     reversed
      integer::
     &     iargcount,ikeycount

      reversed=.false.
      if (present(latest)) reversed = latest

      iargcount=1
      if (present(argcount)) iargcount = argcount

      ikeycount=1
      if (present(keycount)) ikeycount = keycount

      finnode=> inp_fetch_root()
      finnode=> arg_from_context(context,name,finnode,reversed,
     &     ikeycount,iargcount)

      if (present(argcount)) argcount = iargcount
      if (present(keycount)) keycount = ikeycount
      end function





*----------------------------------------------------------------------*
!>    resolves a context string in the registry
!! 
!!    see reg_arg_from_context
*----------------------------------------------------------------------*
      function reg_key_from_context(context,latest,keycount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="reg_key_from_context"
      integer, parameter ::
     &     ntest = 00

      
      type(node), pointer ::
     &     finnode
      character(len=*),intent(in) ::
     &     context
      logical,intent(in),optional::
     &     latest
      integer,intent(in),optional::
     &     keycount
      logical::
     &     reversed
      integer::
     &     ikeycount

      reversed=.false.
      if (present(latest)) reversed = latest

      ikeycount=1
      if (present(latest))ikeycount=keycount


      finnode=>reg_fetch_root()
      finnode=> key_from_context(context,finnode,reversed,ikeycount)

      end function

*----------------------------------------------------------------------*
!>    resolves a context string in the input
!! 
!!    see reg_arg_from_context
*----------------------------------------------------------------------*
      function inp_key_from_context(context,latest,keycount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="inp_key_from_context"
      integer, parameter ::
     &     ntest = 00

      
      type(node), pointer ::
     &     finnode,input_root
      character(len=*),intent(in) ::
     &     context
      logical,intent(in),optional::
     &     latest
      integer,intent(inout),optional::
     &     keycount
      logical::
     &     reversed
      integer::
     &     ikeycount

      reversed=.false.
      if (present(latest)) reversed = latest
      ikeycount=1
      if (present(keycount))ikeycount=keycount

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (latest) write (lulog,*) "reverse search"
         if (.not.latest) write (lulog,*) "normal search"
         write (lulog,'(x,i3,"th occurrence requested")') 
      end if

      input_root=> inp_fetch_root()
      finnode=> key_from_context(context,input_root,reversed,ikeycount)

      if (ntest.gt.100.and. associated(finnode))
     &     write (lulog,*)"success"
      if (ntest.gt.100.and. .not. associated(finnode))
     &     write (lulog,*)"no success"
      if (present(keycount))keycount=ikeycount
      end function

*----------------------------------------------------------------------*
!>    resolves a context and name to a keyword
!!    @param curkey pointer to current keyword
!!    @param[inout] level on in old level on out new level
*----------------------------------------------------------------------*
      function arg_from_context(context,name,tree_root,latest,keycount,
     &     argcount)
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
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
         finnode=>key_getArgument(curkey,name,latest,argcount)
      end if

      if(ntest.ge.100.and. associated(finnode))
     &     write (lulog,*) "success"
      if(ntest.ge.100.and. .not. associated(finnode))
     &     write (lulog,*) "No success"

      end function

*----------------------------------------------------------------------*
!>    resolves a context to a keyword
!!
!!    traverses the tree with a stack reliant depth first search
!!    @param context of the keyword including that keywords name
!!    @param[in] tree_root starting point of search
!!    @param[in] latest for reversed search direction

*----------------------------------------------------------------------*
      function key_from_context(context,tree_root,latest,keycount) 
     &     result(finnode)
*----------------------------------------------------------------------*
      implicit none 
      character(len=*),parameter ::
     &     i_am="key_from_context"

      integer, parameter ::
     &     ntest = 00
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
      ctxt_end=.False.
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
            curnode=> key_getSubkey(finnode,context(ipst:ipnd),latest,
     &           keycount)

            if (ntest.ge.100) then
               if (associated(curnode))write(lulog,*) 
     &              ' looking at "',getAttribute(curnode,atr_name),'"'
               write(lulog,*) ' in final level',ctxt_end
               write(lulog,*) ' keys to be found:',keycount
            end if

            if (keycount.eq.0)then
               finnode=>curnode
               exit
            end if
         else            
            ikeycount=sikeycount
            curnode=>key_getSubkey(finnode,context(ipst:ipnd),latest,
     &           ikeycount)
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
            if (ierr.gt.0) call quit (1,i_am,"stack exceeded:"//
     &           str(ist_stack))
      
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
            finnode=> getParentNode(finnode)
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
!>    retrieves the subkeyword below curkey
!!
!!    @param curkey pointer to parent keyword
!!    @param name name of the subkeyword
!!    @param latest if ocurrences should be counted backwards
!!    @param icount the icountth occurrence is retrieved
*----------------------------------------------------------------------*
      function key_getSubkey(curkey,name,latest,icount) result(nextkey)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="key_getSubkey"
      integer,parameter::
     &     ntest=00

      type(Node),pointer::
     &     nextkey

      type(Node),pointer,intent(in)::
     &     curkey
      character(len=*),intent(in)::
     &     name 
      logical,intent(in)::
     &     latest
      integer,intent(inout)::
     &     icount 

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) " looking for:",trim(name)
         write(lulog,*)  " number to be found:",icount
         if (latest)write(lulog,*) " reversed order"
      end if
      nextkey=> key_getSubnode(curkey,name,key_tag,latest,icount)
      end function








*----------------------------------------------------------------------*
!>    retrieves an argument below curkey
!!
!!    @param curkey pointer to parent keyword
!!    @param name name of the argument
!!    @param latest if ocurrences should be counted backwards
!!    @param icount the icountth occurrence is retrieved
*----------------------------------------------------------------------*
      function key_getArgument(curkey,name,latest,icount)result(nextarg)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="key_getArgument"
      integer,parameter::
     &     ntest=00

      type(Node),pointer::
     &     nextarg

      type(Node),pointer,intent(in)::
     &     curkey
      character(len=*),intent(in)::
     &     name 
      logical,intent(in)::
     &     latest
      integer,intent(inout)::
     &     icount 

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) " looking for:",trim(name)
         write(lulog,*)  " number to be found:",icount
         if (latest)write(lulog,*) " reversed order"
      end if
      nextarg=> key_getSubnode(curkey,name,arg_tag,latest,icount)
      
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
      function key_getSubNode(curkey,name,tag,latest,icount)
     &     result(nextnode)
      implicit none
      character(len=*),parameter::
     &     i_am="key_getSubNode"
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
   

      if (latest)then
         nextnode=>getLastChild(curkey)
      else 
         nextnode=>getFirstChild(curkey)
      end if 

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) " looking for:",trim(name)
         write(lulog,*) " of type:",trim(tag)
         if (latest)write(lulog,*) " reversed order"
         write(lulog,*)  " number to be found:",icount
      end if 

      do 
         if (.not. associated(nextnode)) exit
      if (ntest .gt. 100) then
         write(lulog,*)  " next node:",getAttribute(nextnode,atr_name)
         write(lulog,*)  " type:",getNodeName(nextnode)
         write(lulog,*)  " status:",getAttribute(nextnode,atr_stat)
         write(lulog,*)  " number to be found:",icount
      end if 

         if ( (getNodeName(nextnode).eq.trim(tag) ) .and.
     &        ( getAttribute(nextnode,atr_name) .eq. trim(name) ).and.
     &        (getAttribute(nextnode,atr_stat).eq. status_active))
     &        then
               icount=icount-1
            if (icount.eq.0) then
!     nextnode points correctly
               if (ntest .gt. 100)  write(lulog,*)"success"
               return
            end if
         end if 
         if (latest)then
            nextnode=>getPreviousSibling(nextnode)
         else 
            nextnode=>getNextSibling(nextnode)
         end if
      end do 
      return 
      end function

      end module
