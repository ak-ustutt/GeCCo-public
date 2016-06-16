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
!!          -# name contains the name
!!          -# len contains the length
!!          -# kind contains the number of the kind of argument(see 'par_vtypes.h')
!!          -# note that all attributes are saved as strings. they can be converted via rts

      module parse_input2
      use FoX_dom
      use FoX_common, only:str,rts
      implicit none
      include 'stdunit.h'
      include 'par_vtypes.h'
      integer,parameter::
     &     file_loc_len=28
      character(len=file_loc_len),parameter ::
     &     rel_file_loc="/data/keyword_registry3.xml"
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
      integer,parameter::
     &     DOCUMENT_TYPE=9
      character,parameter::
     &     status_active="A",
     &     status_inactive="I"

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
     &     "input tree found")
      root=> getFirstChild(input_doc)
      
      if (.not. associated(root)) call quit(1,i_am,
     &     "no root element found for input")

      end function


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
     &     input_root
      input_root=>  inp_fetch_root()
      call tree_set_input_status_(input_root, history_pointer,one_more)
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
         write (lulog,*) "reading file:",trim(file)
      end if 

      call reg_init_()
      call inp_start_()

      
      registry_doc => parseFile(trim(file), ex=ex)

      !! @TODO more special error handling 
      if (getExceptionCode(ex) .ne. 0)
     &     call quit(1,i_am,"could not read keyword registry")


      reg_root=> reg_fetch_root() 
      
      call setAttribute(reg_root,atr_name,"registry")
      end subroutine



*----------------------------------------------------------------------*
!!    returns name of the keyword_file
*----------------------------------------------------------------------*
      function get_keyword_file() result(file_name)
*----------------------------------------------------------------------*
      implicit none
      integer,parameter::
     &     ntest= 00
      character(len=16),parameter ::
     &     i_am="get_keyword_file"
      character(len=256)::
     &     path_name,file_name
      integer ::
     &     len


      call get_environment_variable( "GECCO_DIR", value=path_name,
     &     length = len)

      if (len.EQ.0)
     &     call quit(0,i_am,
     &     "Please, set the GECCO_DIR environment variable.")
      if (len .gt.(256-file_loc_len)  )
     &     call quit(0,i_am,
     &     "GECCO_DIR to long, cannot set keyword_registry")
      file_name=trim(path_name)//rel_file_loc
      end function



*----------------------------------------------------------------------*
!>    wrapper for keyword_list for the input
*----------------------------------------------------------------------*
      subroutine inp_show(unit)
*----------------------------------------------------------------------*
      integer,intent(in)::
     &     unit

      type(Node),pointer::
     &     input_root
      input_root => inp_fetch_root()
      call keyword_list(unit,input_root," ",show_args=.True.)
      end subroutine

*----------------------------------------------------------------------*
!>    wrapper for keyword_list for the input
*----------------------------------------------------------------------*
      subroutine reg_show(unit)
*----------------------------------------------------------------------*
      integer,intent(in)::
     &     unit
      type(Node),pointer::
     &     key_root
      key_root => reg_fetch_root()
      call keyword_list(unit,key_root," ",show_args=.True., 
     &     show_status=.False.)
      end subroutine












*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
! self contained subroutines
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

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



      do i=1,2
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
      character(len=16),parameter ::
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
         call inp_show(lulog)
      end if 
      call tree_set_status(in_root,status_inactive)
      one_more=.true.
      if (.not.associated(history)) then
         history =>getFirstChild( in_root)
         
         if (.not.associated(  history)) then
            call quit(0,i_am,'not a single keyword given?')
         end if
         if(ntest.ge.100)
     &        write (lulog,*) "history_pointer initiated"
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
         call inp_show(lulog)
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
      current => keywd
      
      do while(associated(current))
         call setAttribute(current,atr_stat,status)
         current=> key_dsearch(current)
      end do
      end subroutine





*----------------------------------------------------------------------*
!>   finds the icount node of a given context that is active.
!!    for keywords only
!!    
!!   @param[in] tree_root root element of the input_tree 
!!   @param[out] finnode found node (null() if no such node exists)
!!   @param[in] context context of the searched node
!!   @param[in] icount
*----------------------------------------------------------------------*
      subroutine find_active_node(tree_root,finnode,context,icount)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter ::
     &     ntest=00
      character(len=16)::
     &     i_am="find_active_node"
      type(Node), pointer ,intent(in)::
     &     tree_root
      type(Node), pointer ,intent(out)::
     &     finnode
      character ,intent(in)::
     &     context*(*)
      integer,intent(in)::
     &     icount

      character ::
     &     curcontext*1024     
      type(Node), pointer ::
     &     current
      integer ::
     &     jcount
      character::
     &     status
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "looking for",icount,"th",trim(context) 
      end if

      if (.not.hasChildNodes(tree_root))
     &     call quit(1,i_am,'invalid keyword tree')

      current => getFirstChild(tree_root)

      finnode => null()

      jcount = 0
      key_loop: do 
        call rts(getAttribute(current,atr_stat),status)
        if (ntest.ge.100)then 
           write (lulog,*) "keyword:",getAttribute(current,atr_name)
           write (lulog,*) "status:",status
        end if

        if (status.eq.status_active) then
          call keyword_get_context(curcontext,current)
! keyword_get_context only gets the context of the current keyword
          if (len_trim(curcontext).eq.0)then
             curcontext=getAttribute(current,atr_name)
          else
             curcontext=trim(curcontext)//"."//
     &            getAttribute(current,atr_name)
          end if 
          
          if (trim(context).eq.trim(curcontext)) jcount = jcount+1 
          if (icount.eq.jcount) then
            finnode => current
            exit key_loop
         end if
        end if
        current => key_dsearch(current)
        if (.not.associated(current)) exit key_loop
      end do key_loop
            if (ntest.ge.100)
     &     write (lulog,*) i_am," has found",jcount,
     &     "occurences of ",trim(context),"and returning:",
     &     associated(finnode)
      return
      end subroutine



*----------------------------------------------------------------------*
!>    find keyword in first sublevel or up
!!
!!
!!    @param cur_node pointer to current node
!!    @param[out] nxt_node found node or unassociated
*----------------------------------------------------------------------*
      subroutine next_node(cur_node,nxt_node,key)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      character(len=9),parameter::
     &     i_am="next_node"
      integer, parameter ::
     &     ntest =00

      type(node), pointer ::
     &     cur_node
      type(node), pointer ::
     &     nxt_node
      character, intent(in) ::
     &     key*(*)

      type(node), pointer ::
     &     current,listnode
      type(NodeList), pointer ::
     &     keylist
      integer ::
     &     ii
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' start = "',getAttribute(cur_node,atr_name),'"'
         write(lulog,*) ' search for "',trim(key),'"'
      end if

      nxt_node=>null()
      current => cur_node



      if (hasChildNodes(current))then
         keylist=>getChildNodes(current)
      else
         current=>getParentNode(current)
         keylist=>getChildNodes(current)
      end if 
   
      node_loop: do
         do ii=0,getLength(keylist)-1      ! No,No,No
            listnode=>item(keylist,ii)
            if (getNodeName(listnode) 
     &           .eq. key_tag .and. 
     &           getAttribute(listnode, atr_name)
     &           .eq. key ) then 
               nxt_node=> listnode
               exit node_loop
            end if

         end do

         if (getNodeName(current).ne. key_root_tag)then
            current=>getParentNode(current)
            keylist=>getChildNodes(current)
         else
            exit node_loop
         end if
      end do node_loop

      if (ntest.ge.100) then
         if (associated(nxt_node)) write(lulog,*) 'success'
         if (.not.associated(nxt_node)) write(lulog,*) 'no success'
      end if

      end subroutine next_node



*----------------------------------------------------------------------*
!Output subroutine; to be implemented later
*----------------------------------------------------------------------*
      subroutine keyword_list(luwrt,tree_root,
     &     context,n_descent,show_args,show_status)
*----------------------------------------------------------------------*
      include 'stdunit.h'
      character(len=7),dimension(8),parameter::
     &     type_array=(/"logical","integer","unknown","real",
     &     "unknown","unknown","unknown","string"/)
      integer, parameter ::
     &     ntest=00
      character(len=17),parameter ::
     &     i_am="keyword_list"

      integer, intent(in)::
     &     luwrt
      type(Node),pointer, intent(in)::
     &     tree_root
      character(len=*),intent(in)::
     &     context
      integer,optional,intent(in)::
     &     n_descent
      logical, optional, intent(in)::
     &     show_args,show_status

      character(len=64)::
     &     fmtstr
      logical:: 
     &     args_vis,status_vis
      integer::
     &     level, 
     &      type, dim,
     &     ii
      character::
     &     status

      type(Node),pointer::
     &     curkey,curarg
      type(Nodelist),pointer::
     &     child_list
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of key_root:"
     &        ,associated(tree_root)

      end if

      args_vis=.false.
      if (present(show_args)) args_vis=show_args
      status_vis=.True.
      if (present(show_status)) status_vis=show_status

      if (ntest.ge.100) then
         write (lulog,*) "show_status?",status_vis
         write (lulog,*) "show_arguments?",args_vis
      end if
!      if (present(n_descent)) levels=n_descent

      curkey=>null()
      curarg=>null()
      if (len_trim(context).eq.0)then 
         curkey=> getFirstChild(tree_root)
         if (.not. associated(curkey))then
            write(luwrt,*) "no keywords_set"
            return 
         end if
         level=level+1
      else
         call find_node(tree_root,curkey,context)
         if (.not. associated(curkey)) call quit(1,i_am,
     &        "starting node not found")
      end if 

      level=0
      key_loop: do 
         if (status_vis)then 
            status=getAttribute(curkey,atr_stat)
            if (ntest.ge.100)then 
               write (lulog,*) "status:",getAttribute(curkey,atr_stat)
            end if 
            if (status.eq.status_active) then
               write(fmtstr,'("(""A"",",i3,"x,a)")') 2*level+1
            else
               write(fmtstr,'("(""I"",",i3,"x,a)")') 2*level+1
            end if
         else 
             write(fmtstr,'("("">"",",i3,"x,a)")') 2*level+1
         end if
!!       @TODO print comments out
         write(luwrt,fmtstr) getAttribute(curkey,atr_name)
         
         if (hasChildNodes(curkey) .and. args_vis)then 
            child_list=>getChildNodes(curkey)
            arg_loop :do ii=0,getLength(child_list)-1
               curarg=> item(child_list,ii)
               if (.not. getNodeName(curarg).eq. arg_tag )cycle arg_loop
               call rts(getAttribute(curarg,atr_kind),type)
               call rts(getAttribute(curarg,atr_len),dim)
               write(fmtstr,'("(x,",i3,"x,a,x,i2,x,a)")') 2*level+4
               write(luwrt,fmtstr) getAttribute(curarg,atr_name)//" "
     &              //trim(type_array(type))//" of len",dim,": "// 
     &              getAttribute(curarg,atr_val)
            end do arg_loop 
         end if
         curkey=> key_dsearch(curkey,level)
         if (.not.associated(curkey)) exit key_loop
         if (present(n_descent))then
            if (level.gt.n_descent) exit key_loop
         end if
         
      end do key_loop
      end subroutine











*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
!     Some general use subroutines
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
!>    returns the context of a given keyword
*----------------------------------------------------------------------*
      subroutine keyword_get_context(curcontext,current,full)
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
     &     current 
      logical,intent(in),optional::
     &     full

      character(len=name_len)::
     &     curname
      
      type(node), pointer ::
     &     internal

      curcontext=" "
      if (present(full))then
         if (full) curcontext=getAttribute(current,atr_name)//"."
      end if 
      internal=> getParentNode(current)
      do while (getNodeName(internal) .ne. key_root_tag)
         curcontext=getAttribute(internal,atr_name)//
     &   "."//trim(curcontext)
         internal=> getParentNode(internal)
      end do
      curcontext(len_trim(curcontext):len_trim(curcontext))=" "
      end subroutine

*----------------------------------------------------------------------*
!>     look whether keyword cur_key hosts an argument with key "key"
!!     and return the corresponding node
!!    @param[out] arg found argument node (null if no node was found)
!!    @param[in] cur_key current keyword, where children are searched
!!    @param[in] key name of the argument
*----------------------------------------------------------------------*
      subroutine arg_node(arg,cur_key,key)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      character(len=9),parameter::
     &     i_am="arg_node"
      integer, parameter ::
     &     ntest =00

      type(Node), pointer ::
     &     cur_key
      character, intent(in) ::
     &     key*(*)
      type(Node), pointer, intent(out) ::
     &     arg
      type(NodeList),pointer ::
     &     nodes_list
      type(Node), pointer ::
     &     curnode
      integer ::
     &     ii
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' looking for key "',key,'"'
         write(lulog,*) ' association status:',associated(cur_key)
         write(lulog,*) ' current_keyword "',
     &        getAttribute(cur_key,atr_name),'"'
      end if

      arg => null()
      if (hasChildNodes(cur_key))then
         nodes_list=> getChildNodes(cur_key)
         do ii=0,getLength(nodes_list)-1 ! this ain't no true Fortran
            curnode=>item(nodes_list,ii)
            if (ntest.ge.100) then
               write(lulog,*) ' looking at "',
     &              getAttribute(curnode,atr_name),'"'
               write(lulog,*) ' of type:',getNodeName(cur_key)
            end if 
            if (getNodeName(curnode)
     &           .eq. arg_tag .and. 
     &           getAttribute(curnode, atr_name)
     &           .eq. key ) then
               arg => curnode
               exit
            end if 
         end do 
      end if 


      if (ntest.ge.100) then
        if (associated(arg)) write(lulog,*) 'success'
        if (.not.associated(arg)) write(lulog,*) 'no success'
      end if
      end subroutine



*----------------------------------------------------------------------*
!>     navigates to a specific keyword
!!
!!     @param[in] tree_root root element from which the search starts, must be associated.
!!     @param[out] finnode either pointing to the node given by context or not associated if no keyword given by context exists.
!!     @param[in] context is a string as e.g. "key.subkey.subsubkey"
!!     @param[in] latest optional parameter if .True. subkeys on the same level are searched in reverse order
*----------------------------------------------------------------------*
      subroutine find_node(tree_root,finnode,context,latest)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      character(len=9),parameter ::
     &     i_am="find_node"
      integer, parameter ::
     &     ntest = 00

      type(node), pointer ::
     &     tree_root
      type(node), pointer ::
     &     finnode
      character, intent(in) ::
     &     context*(*)
      logical, intent(in), optional ::
     &     latest

      logical :: 
     &     forward
      integer ::
     &     ipst, ipnd, len, ii

      type(NodeList), pointer ::
     &     nodes_list
      type(node), pointer ::
     &     current,tmpnode
      logical ::
     &     found
      
      forward = .true.
      if (present(latest)) forward = .not.latest
      print *, "find_node called on", getAttribute(tree_root,atr_name)
      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' context = "',trim(context),'"'
        if (forward) write(lulog,*) ' forward search'
        if (.not.forward) write(lulog,*) ' backward search'
      end if 
      ii=1
      finnode => null()
      current => tree_root
      finnode =>  key_from_context(context,tree_root,.not.forward,
     &     ii)

      if (ntest.ge.100) then
        if (associated(finnode)) write(lulog,*) 'success'
        if (.not.associated(finnode)) write(lulog,*) 'no success'
      end if
      end subroutine









*----------------------------------------------------------------------*
      subroutine get_argument_dimension_core(curarg,num,type,succ)
*----------------------------------------------------------------------*
      use Fox_common, only: rts
      include 'par_vtypes.h'
      include 'stdunit.h'
      integer,parameter::
     &     ntest=00
      character(len=27),parameter ::
     &     i_am="get_argument_dimension_core"

      type(Node), pointer,intent(inout) ::
     &     curarg
      integer , intent(out)::
     &     num,type
      logical , intent(out)::
     &     succ
     
      logical ::
     &     lval 
      logical,allocatable::
     &     larr(:)
      integer::
     &     ival
      integer,allocatable::
     &     iarr(:)
      real(8)::
     &     xval
      real(8),allocatable::
     &     xarr(:)
      integer ::
     &     dim_tot,ex

      call rts(getAttribute(curarg,atr_len),dim_tot)
      call rts(getAttribute(curarg,atr_kind),type)
      if (ntest.ge.100) then 
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,'(" dim_tot:",i3)')dim_tot 
         write(lulog,'(" type:",i3)')type
         write(lulog,'("unconverted input:",a,":")')
     &        getAttribute(curarg,atr_val)
      end if 

         succ=.false.
      if(.not.hasAttribute(curarg,atr_val))then 
         num=0
         return
      else
         num=1
      end if 

      select case(type)
      case (vtyp_log)
         allocate(larr(dim_tot))
         call rts(getAttribute(curarg,atr_val),larr
     &        ,iostat=ex,num=num) 
         if (ex.le. 0) succ = .true.
      case (vtyp_int)
         allocate(iarr(dim_tot))
         call rts(getAttribute(curarg,atr_val),iarr(1:dim_tot)
     &        ,iostat=ex,num=num) 
         if (ex.le. 0) succ = .true.
      case (vtyp_rl8)
         allocate(xarr(dim_tot))
         call rts(getAttribute(curarg,atr_val),xarr
     &        ,iostat=ex,num=num) 
         if (ex.le. 0) succ = .true.
      case (vtyp_str)
         num=dim_tot
         succ = .true.
      end select
      return 
      end subroutine 
























*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!     Iterating search functions
!     if called repeatedly these functions will iterate over the remaining tree in a depths first search
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*



*----------------------------------------------------------------------*
!>    keyword_specific version of the dsearch function
!!
!!
!!    @param curkey pointer to current keyword
!!    @param[inout] level on in old level on out new level
*----------------------------------------------------------------------*
      function key_dsearch(curkey, level) result(nxtkey) 
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="key_dsearch"

      type(Node), pointer::
     &     nxtkey

      type(Node), pointer,intent(in)::
     &     curkey

      integer,intent(inout),optional ::
     &     level

      integer::
     &     dlevel !levelchange
      
      nxtkey=>curkey

      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) "starting at:",
     &        getAttribute(nxtkey,atr_name)
      end if 

      dlevel=0
      main_loop: do
         nxtkey=>elem_dsearch(nxtkey,dlevel)

         if(associated(nxtkey).and. ntest .ge.100) then 
            write(lulog,*) "looking at:",
     &           getAttribute(nxtkey,atr_name)
            write(lulog,*) " "
         else if (ntest.ge.100) then
            write(lulog,*) "not associated"
         end if

         if (.not.associated(nxtkey))exit main_loop

         if (getNodeName(nxtkey).eq. trim(key_tag))then
            exit main_loop
         else if (getNodeName(nxtkey).eq.key_root_tag) then
            nxtkey=>null()
            exit main_loop
         end if
      end do main_loop

      if (present(level))level=level+dlevel
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
         write (lulog,*) " found node"
      else if (ntest.ge.100)then
         write (lulog,*) " no node found"
      end if 
      end function


*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
!     guided search functions 
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
! external subroutines should only be concerned with 
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
     &     i_am="inp_key_from_context"
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

      reversed=.false.
      if (present(latest)) reversed = latest
      
      finnode=> inp_fetch_root()
      finnode=> arg_from_context(context,name,finnode,reversed,
     &     ikeycount,iargcount)

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

      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (latest) write (lulog,*) "reverse search"
         if (.not.latest) write (lulog,*) "normal search"
         write (lulog,'(x,i3,"th occurrence requested")') 
      end if

      finnode=> inp_fetch_root()
      finnode=> key_from_context(context,finnode,reversed,ikeycount)

      if (ntest.gt.100.and. associated(finnode))
     &     write (lulog,*)"success"
      if (ntest.gt.100.and. .not. associated(finnode))
     &     write (lulog,*)"no success"
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
     &     max_stack=16             ! 
      
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
     &     keyc_stack(max_stack), ikeyc_stack,
     &     ierr
      logical ::
     &     ctxt_end
      
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

      call stack_push(st_stack,ist_stack,0,ierr)
      call stack_push(keyc_stack,ikeyc_stack,0,ierr)
      if (ierr.gt.0) call quit(1,i_am,"no stack?")
      do
         ipnd = index(context(ipst:),context_sep)+ipst-2
         if (ipnd.lt.ipst) ctxt_end=.True.
         if(ctxt_end) ipnd=len



         if (ntest.ge.100) then
            write(lulog,*) ' current subkeyword: "',
     &           context(ipst:ipnd),'"'
            write(lulog,*) ' in final level',ctxt_end
            write(lulog,*) ' keys to be found:',keycount
         end if
        
         if (ctxt_end)then
            ! 
            ! key_getSubkey looks for the keycount'th subnode
            !  AND reduces keycount by the amount of found nodes!!!
            curnode=>key_getSubkey(finnode,context(ipst:ipnd),latest,
     &           keycount)
            if (ntest.ge.100) then
               write(lulog,*) ' currently below: "',
     &              getAttribute(finnode,atr_name),'"'
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



         ! which direction to go
         if (associated(curnode))then 
            ! go one level down
            if (ntest .gt. 100) then
               write(lulog,*) "found",getAttribute(curnode,atr_name)
               write(lulog,*) "going one level up"
            end if 

            call stack_push(st_stack,ist_stack,ipst,ierr)
            call stack_push(keyc_stack,ikeyc_stack,sikeycount,ierr)
            if (ierr.gt.0) call quit (1,i_am,"stack exceeded:")
      
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
     &           "OK, WTF?? trying to access negative stack")
            finnode=> getParentNode(finnode)
            sikeycount=sikeycount+1
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

      subroutine stack_pop(stack,stack_pointer,i,ierr)
      implicit none
      integer,intent(inout)::
     &     stack(max_stack),stack_pointer
      integer,intent(out)::
     &     i
      integer, intent(out)::
     &     ierr
  
      stack_pointer=stack_pointer-1
      if (stack_pointer.le.0) then 
         ierr=ierr+1
         return
      end if 
      i=stack(stack_pointer)
      end subroutine

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
