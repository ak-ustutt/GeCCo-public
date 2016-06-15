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
      include 'write_styles.h'
      include 'par_vtypes.h'
      integer,parameter::
     &     file_loc_len=23
      character(len=file_loc_len),parameter ::
     &     rel_file_loc="/data/keyword_registry2"
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
      

      character,parameter::      !tags 
     &     key_root_tag*8="key_root",
     &     key_tag*7="keyword",
     &     arg_tag*8="argument"



      type(Node), pointer :: 
     &     registry_doc,         !root of the tree that represents the keyword registry
     &     key_root,            !root element for keywords in registry
     &     input_doc,           !document element of input
     &     input_root,          !root element of input; corresponds to key_root
     &     history_pointer      ! to  this "calculate" block the calculation has advanced

      

      contains 

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
! Subroutines coupled to the module variables
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
!>   wrapper to uncouple set_input_status_ from module variables
!!
!!  @param  one_more logical to show if there are unprocessed blocks remaining.
!!          input status not used
*----------------------------------------------------------------------*
      subroutine process_input_(one_more)
*----------------------------------------------------------------------*
      implicit none
      logical, intent(inout)::
     &     one_more
      type(Node),pointer::
     &     history

      call set_input_status_(input_root, history_pointer,one_more)
      end subroutine 





*----------------------------------------------------------------------*
!!    initializes the module variables relating to the registry
*----------------------------------------------------------------------*
      subroutine keyword_init_()
*----------------------------------------------------------------------*
      implicit none
      registry_doc =>null()
      key_root=> null()
      end subroutine

*----------------------------------------------------------------------*
!!    initializes the module variables relating to the input
*----------------------------------------------------------------------*
      subroutine input_init_()
*----------------------------------------------------------------------*
      implicit none
      input_doc=> null()
      input_root=>null()
      history_pointer=> null()
      end subroutine

*----------------------------------------------------------------------*
!!    creates the document to hold the input variables
*----------------------------------------------------------------------*

      subroutine input_start_()
      call input_init_()

      !easiest way to create a document
      input_doc=> parseString(
     &     "<"//key_root_tag//"></"//key_root_tag//">")
      input_root=> getFirstChild(input_doc)
      call setAttribute(input_root,atr_name,"input")
      end subroutine input_start_

*----------------------------------------------------------------------*
!!    parses the keyword_file and initializes the input
!!  
*----------------------------------------------------------------------*
      subroutine set_keywords()
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=13),parameter ::
     &     i_am="set_keywords"
      integer, parameter ::
     &     ntest=00
      type(DOMException)::
     &     ex
      type(Node),pointer ::
     &     tmpnode
      type(NodeList),pointer::
     &     nodes_list

      call keyword_init_()
      call input_start_()


      registry_doc => parseFile(trim(get_keyword_file()), ex=ex)
      !! @TODO more special error handling 
      

      if (getExceptionCode(ex) .ne. 0)
     &     call quit(0,i_am,"could not read keyword registry")

      if (ntest .ge. 100) write (lulog,*) "file parsed"
      nodes_list=>getElementsByTagName(registry_doc, key_root_tag)    
      !  assumes there is only on of these elements
      if (getLength(nodes_list) .ne. 1)
     &     call quit(0,i_am,"more than one key_root element")

      if (ntest .ge. 100) then
         write (lulog,'("expected: 1 element")') 
         write (lulog,'("found: ",i2,"element")')getLength(nodes_list)
      end if
      key_root=> item(nodes_list, 0) ! 0-indexed ... that's no true Fortran
      
      call setAttribute(key_root,atr_name,"registry")
      contains

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

      end subroutine 

*----------------------------------------------------------------------*
!!   wrapper to uncouple keyword_parse_ from module variables
!!
!!  @param  luin 
*----------------------------------------------------------------------*
      subroutine keyword_parse(luin)
*----------------------------------------------------------------------*
      implicit none
      integer, intent(in) ::
     &     luin

      call keyword_parse_(luin,input_doc,key_root)
      end subroutine keyword_parse

      subroutine show_input()
      include "stdunit.h"
      call keyword_list(lulog,input_root," ",show_args=.True.)
      end subroutine


      subroutine show_registry()
      include "stdunit.h"
      call keyword_list(lulog,key_root," ",show_args=.True., 
     &     show_status=.False.)
      end subroutine












*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
! self contained subroutines
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     set all previous keywords with same key as present keyword to
*     inactive (+ all sub-levels)
*----------------------------------------------------------------------*
      subroutine unset_previous_keywords(keywd)
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
     &     call set_keyword_status(current,-1)
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
!!   @param one_more logical if unprocessed blocks remain
*----------------------------------------------------------------------*
      subroutine set_input_status_(in_root, history,one_more)
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
      end if 
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
         call set_keyword_status(current,+1)
         call unset_previous_keywords(current)



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
      
      end subroutine 

*----------------------------------------------------------------------*
!>    sets the given keyword and all subkeywords to a given status
!!
!!   @param keywd pointer to a keyword
!!   @param status integer with the status we set it to 
*----------------------------------------------------------------------*
      subroutine set_keyword_status(keywd,status)
*----------------------------------------------------------------------*
      implicit none
      type(Node), pointer, intent(in)::
     &     keywd
      integer,intent(in)::
     &     status
      type(Node),pointer::
     &     current 
      current => keywd
      
      do while(associated(current))
         call setAttribute(current,atr_stat,str(status))
         call dsearch_next_key(current,key_tag)
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
     &     jcount, status
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

        if (status.gt.0) then
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
        call dsearch_next_key(current, key_tag)
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
!>   creates a node (keyword or argument)
!!
!!   copies all attributes from template but does not copy child nodes and 
!!     inserts the new node in document doc with the context of the original node.
!!     context must already exist in doc.
!!     the node is set inactive
!!   overrides the argument value with value. for keywords value is irrelevant
!!   @param[inout] doc document that will own the newly created node. 
!!   @param[in] template template for the newly created node
!!   @param[in] value string which can be converted to input
!!   @TODO error checking
*----------------------------------------------------------------------*
      subroutine create_node(doc,template,value)
*----------------------------------------------------------------------*
      implicit none 
      include "stdunit.h"
      character(len=11),parameter::
     &     i_am="create_node"
      integer,parameter ::
     &     ntest=00
      integer,parameter::
     &     ERR_TO_MANY_ELEMENTS=1,
     &     ERR_UNCONVERTIBLE=2

      integer, parameter ::
     &     maxlen  = 256

      type(node), pointer ::
     &     doc,template
      character(len=*),intent(in)::
     &     value
      type(node),pointer ::
     &     new_elem,new_parent,root
      character(len=maxlen) ::
     &     context
      logical::
     &     full_context
      integer::
     &     kind,dim,ierr 

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of doc:"
     &        ,associated(doc)
         write (lulog,*) "testing association of template:"
     &        ,associated(template)
         write (lulog,*) "tag, name of template:",
     &        getNodeName(template),getAttribute(template,atr_name)
      end if

      if (.not. associated(doc)) 
     &     call quit(1,i_am,"doc not set")
      
      if (.not. associated(template)) 
     &     call quit(1,i_am,"template not set")
      

      new_elem=>createElement(doc,getNodeName(template))

      if (ntest.ge.100) then
         write (lulog,*) "tag of new element:",
     &        getNodeName(new_elem)
      end if 
      context=" "
      call  keyword_get_context(context,template)
      if (ntest.ge.100) then
         write (lulog,*) "creating a new keyword/argument on context:",
     &   trim(context)
      end if


!> @TODO this assumes that key_root is the first child of doc
      root => getFirstChild(doc)
      call find_node(root,new_parent,trim(context),latest=.True.)
      if (ntest.ge.100) then
         write (lulog,*) "tag, name of parent:",
     &        getNodeName(new_parent),getAttribute(new_parent,atr_name)
      end if

      new_elem=>appendChild(new_parent, new_elem)
      new_parent=> getParentNode(new_elem)
      call setAttribute(new_elem,atr_name,
     &     getAttribute(template,atr_name))

      if (hasChildNodes(new_parent))then 
         root=> getFirstChild(new_parent)
         do while(associated (root))
            root => getNextSibling(root)
         end do 
      end if 
      if (getNodeName(new_elem) .eq. arg_tag)then 

         call setAttribute(new_elem,atr_len
     &        ,getAttribute(template,atr_len))
         call setAttribute(new_elem,atr_kind,
     &        getAttribute(template,atr_kind))
         call rts(getAttribute(template,atr_kind),kind)
         call rts(getAttribute(template,atr_len),dim)
         if (kind.eq.vtyp_log)then
            call setAttribute(new_elem, atr_val,
     &           trim(conv_logical_inp(value,dim,ierr)))
         else 
            call setAttribute(new_elem, atr_val,value)
         end if 
         select case (ierr)
            case (ERR_TO_MANY_ELEMENTS)
               call quit(0,i_am,
     &              "trying to set to many elements for:"//
     &              getAttribute(new_elem,atr_name)//" "//
     &              value)
            case (ERR_UNCONVERTIBLE)
               call quit(0,i_am,
     &              "value for "//getAttribute(new_elem,atr_name)//
     &              " not convertible:"//value)
         end select
      end if
      call setAttribute(new_elem,atr_stat,str(-1))



      end subroutine

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
     &     status, type, dim,
     &     ii

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
            call rts(getAttribute(curkey,atr_stat),status)
            if (status.gt.0) then
               write(fmtstr,'("(""A"",",i3,"x,a)")') 2*level+1
            else
               write(fmtstr,'("(""I"",",i3,"x,a)")') 2*level+1
            end if
         else 
             write(fmtstr,'("("">"",",i3,"x,a)")') 2*level+1
         end if
            
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
         call dsearch_next_key(curkey,key_tag,level)
         if (.not.associated(curkey)) exit key_loop
         if (present(n_descent))then
            if (level.gt.n_descent) exit key_loop
         end if
         
      end do key_loop
      end subroutine
*----------------------------------------------------------------------*
*     parse the keywords on unit luin
*     the unit should be a formatted, sequential file, positioned
*     at the place where the parser should start
*----------------------------------------------------------------------*
      subroutine keyword_parse_(luin,in_doc,k_root)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest=00
      character(len=17),parameter ::
     &     i_am="keyword_parse"


      integer, intent(in) ::
     &     luin

      type(Node), pointer ::
     &     in_doc,               ! document element of the input
     &     k_root                ! root element of the preset keyword tree
      integer, parameter ::
     &     maxlen  = 256,
     &     n_delim = 8
      logical, parameter ::
     &     new_line_is_delim = .true.
      character, parameter ::
     &     delimiter(n_delim) = 
     &     (/' ', ';', ',', '(', ')', '"', '!', '=' /)

      integer, parameter ::
     &     ispace =   1,
     &     isemicolon = 2,
     &     icomma =   3,
     &     iparen_o = 4,
     &     iparen_c = 5,
     &     iquote   = 6,
     &     icomment = 7,
     &     iequal   = 8

      integer, parameter ::
     &     n_allowed_start = 1,
     &     allowed_start(n_allowed_start) = (/icomment/),
     &     n_allowed_after_key = 6,
     &     allowed_after_key(n_allowed_after_key) = 
     &     (/ispace,iequal,isemicolon,
     &     iparen_o,iquote,icomment/),
     &     n_allowed_after_arg = 4,
     &     allowed_after_arg(n_allowed_after_arg) = 
     &     (/isemicolon,iequal,icomma,icomment/)

      character ::
     &     line*(maxlen)


      type(Node), pointer ::
     &     curkey,nxtkey
      type(Node), pointer ::
     &     curarg
      character(len=maxlen) ::
     &     context
      integer ::
     &     allowed_delim(n_delim), n_allowed_delim
      integer ::
     &     ipst, ipnd, itest, lenline, ierr

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of key_root:"
     &        ,associated(k_root)
      end if

      context = " "
      curkey => k_root

      allowed_delim(1:n_allowed_start)=allowed_start(1:n_allowed_start)
      n_allowed_delim = n_allowed_start

      ierr = 0
      file_loop: do
         read(luin, '(a)', end=100, err=200) line
         ipst = 1
         if (ntest .ge. 100)then 
            write (lulog,*) "reading line"
            write (lulog,*) "'",line,"'"
         end if 
         call clean_line(line,delimiter,n_delim)
         
         lenline=len_trim(line)

         !empty line?

         if (lenline.le.0) cycle file_loop 

         !ensure that the line does not begin with a delimiter
         itest = next_delim(line(ipst:ipst),
     &        allowed_delim,n_allowed_delim)

         if (itest.le.0)then
            ierr = ierr+1
            call error_delim(line,ipst)
            ipst = ipst+1
         end if

         ! comment?       
         if (line(ipst:ipst).eq.delimiter(icomment)) cycle file_loop 
         
         ! any delim that can end an argument or keyword name
         allowed_delim(1:n_allowed_after_key) =
     &        allowed_after_key(1:n_allowed_after_key)
         n_allowed_delim = n_allowed_after_key        


         line_loop: do while(ipst.le.lenline)
            itest = next_delim(line(ipst:),
     &           allowed_delim,n_allowed_delim)
            
            ipnd = abs(itest)+ipst-2
            if (ntest .gt. 100)then 
               write (lulog, *) "disassembling line"
               write (lulog, *) "start_index: ",ipst
               write (lulog, *) "current word: ",line(ipst:ipnd)
            end if 
            
            if (itest.le.0) then
               ierr = ierr+1
               call error_delim(line,ipnd+1)
            end if
            
! is it an argument key?
            call arg_node(curarg,curkey,line(ipst:ipnd))
            if (ntest.ge.100) 
     &           write (lulog,*) "This is an argument?",
     &           associated(curarg)
            if (associated(curarg)) then 
!     check that a value is assigned
               if(  next_delim(line(ipnd+1:ipnd+1),
     &              [iequal],1)
     &              .le. 0 ) then 
                  ierr = ierr+1
                  call error_misseq(line,ipnd+1)
               end if

               
               ipst = ipnd+2
               if (ipst.le.lenline) then
                  if (line(ipst:ipst).eq.delimiter(iparen_o)) then
                     ipnd=index(line(ipst:),delimiter(iparen_c))+ipst-1
                     if (itest.le.0) then
                        ierr = ierr + 1
                        call error_misspc(line,ipnd)
                        exit line_loop
                     end if
                     
                     allowed_delim(1:n_allowed_after_arg) =
     &                    allowed_after_arg(1:n_allowed_after_arg)
                     n_allowed_delim = n_allowed_after_arg
                     call create_node(in_doc,curarg,line(ipst+1:ipnd-1))
                  else
                     allowed_delim(1:n_allowed_after_arg) =
     &                    allowed_after_arg(1:n_allowed_after_arg)
                     n_allowed_delim = n_allowed_after_arg
                     itest = next_delim(line(ipst:),
     &                    allowed_delim,n_allowed_delim)
                     ipnd = abs(itest)+ipst-1
                     if (itest.le.0) then
                        ierr = ierr+1
                        call error_delim(line,ipnd)
                     end if
                     ipnd = ipnd-1
                     call create_node(in_doc,curarg,line(ipst:ipnd))
                  end if 

               else
                  ierr = ierr+1
                  call error_eol(line,ipst)
                  allowed_delim(1:n_allowed_after_key) =
     &                 allowed_after_key(1:n_allowed_after_key)
                  n_allowed_delim = n_allowed_after_key
                  exit line_loop
               end if
            else                ! keyword
               call next_node(curkey,nxtkey,line(ipst:ipnd))
               if (ntest.ge.100) 
     &              write (lulog,*) "Is it a keyword?",
     &              associated(nxtkey)
               if (.not.associated(nxtkey)) then
                  ierr = ierr+1
                  call error_keywd(line,ipst,curkey)
               else
                  curkey => nxtkey
                  call create_node(in_doc,curkey," ")
                  
!     add node to keyword history
               end if
            end if

            ipst = ipnd+2
         end do line_loop
      end do file_loop

 100  continue                  ! read EOF reached
      if (ierr.gt.0)
     &     call quit(0,i_am
     &                ,'input errors detected, see above')
      return 

 200  continue                  !  read err
      call quit(0,i_am,'I/O error on reading input file')
      return
*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*
*     local functions
*     variables of main function can be accessed 
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
!>     find next delimiter
!!    
!!    accesses global par delimiter and n_delim
!!    @param str string to be searched
!!    @param delim(n_delim) array of allowed delimiters
!!    @param n_delim
!!    return value:
!!      position of next delimiter
!!      len+1 if line end is delimiter
!!      negative value, if delimiter is not allowed in context
!!    
*----------------------------------------------------------------------*
      pure integer function next_delim(str,delim,nrdelim)
*----------------------------------------------------------------------*
      implicit none
      character, intent(in) ::
     &     str*(*)
      integer,intent(in)::
     &     nrdelim
      integer,dimension(nrdelim),intent(in)::
     &     delim

      logical ::
     &     ok
      integer ::
     &     ipos, jpos, len, idelim, jdelim
c dbg
CC remove the pure keyword if you use this
c       write (lulog,*)"debug: next_delim"
c       write (lulog,*) "called with: str='",trim(str),"'"
c       write (lulog,*) "called with: n_delim=",nrdelim,"'"
c       write (lulog,*) "delim:",delim(1)
c dbgend 
      len = len_trim(str)
      if (len.eq.0) then
        next_delim = 0
        return
      end if

      ipos = len+1
      idelim = 0
      do jdelim = 1, n_delim
         
        jpos = index(trim(str),delimiter(jdelim))
        if (jpos.gt.0.and.jpos.lt.ipos) then
          ipos = jpos
          idelim = jdelim
        end if
        if (ipos.eq.1) exit
      end do
c dbg
CC remove the pure keyword if you use this
c       write (lulog,*)"length:",len
c       write (lulog,*) "found delimiter at",ipos,"'"
c       if (idelim.gt.0)
c     &      write (lulog,*) "found  delimiter='",
c     &      delimiter(idelim),"'"
c dbgend 
      ok = .true.
      if (idelim.gt.0) then
        ok = .false.
        do jdelim = 1, nrdelim
          if (idelim.eq.delim(jdelim)) then
            ok = .true.
            exit
          end if
        end do
      end if

      if (ok) next_delim = ipos
      if (.not.ok) next_delim = -ipos

      return
      end function
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
!     Error notifier
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
      subroutine error_delim(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character(len=80) ::
     &     fmtstr, fmtstr2

      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(fmtstr2,'("(x,a,",i3,"(a1,x))")') n_allowed_delim
 
      write(lulog,'(x,a)') trim(line)
      write(lulog,fmtstr) 
      write(lulog,fmtstr2) 'INPUT ERROR: unexpected delimiter, '//
     &     'expected one of ',
     &     delimiter(allowed_delim(1:n_allowed_delim))
      if (lulog.ne.luout) then
        write(luout,'(x,a)') trim(line)
        write(luout,fmtstr) 
        write(luout,fmtstr2) 'INPUT ERROR: unexpected delimiter, '//
     &     'expected one of ',
     &     delimiter(allowed_delim(1:n_allowed_delim))
      end if

      return
      end subroutine
*----------------------------------------------------------------------*
      subroutine error_eol(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character(len=80) ::
     &     fmtstr

      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1

      write(lulog,'(x,a)') trim(line)
      write(lulog,fmtstr) 
      write(lulog,'(x,a)') 'INPUT ERROR: unexpected end-of-line'
      if (lulog.ne.luout) then
        write(luout,'(x,a)') trim(line)
        write(luout,fmtstr) 
        write(luout,'(x,a)') 'INPUT ERROR: unexpected end-of-line'
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine error_misspc(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character(len=80) ::
     &     fmtstr

      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1

      write(lulog,'(x,a)') trim(line)
      write(lulog,fmtstr) 
      write(lulog,'(x,a)') 'INPUT ERROR: missing )'
      if (lulog.ne.luout) then
        write(luout,'(x,a)') trim(line)
        write(luout,fmtstr) 
        write(luout,'(x,a)') 'INPUT ERROR: missing )'
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine error_misseq(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character(len=80) ::
     &     fmtstr

      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1

      write(lulog,'(x,a)') trim(line)
      write(lulog,fmtstr) 
      write(lulog,'(x,a)') 'INPUT ERROR: missing ='
      if (lulog.ne.luout) then
        write(luout,'(x,a)') trim(line)
        write(luout,fmtstr) 
        write(luout,'(x,a)') 'INPUT ERROR: missing ='
      end if

      return
      end subroutine


*----------------------------------------------------------------------*
      subroutine error_keywd(str,ipos,curkey)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      type(Node), intent(in) ::
     &     curkey
      character(len=80) ::
     &     fmtstr

      if(abs(ipos).EQ.1) then
       write(fmtstr,'("(x,""^"")")')
      else
       write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      end if
      write(lulog,'(x,a)') trim(line)
      write(lulog,fmtstr) 
      write(lulog,'(x,a)') 'INPUT ERROR: unexpected keyword '
      if (lulog.ne.luout) then
        write(luout,'(x,a)') trim(line)
        write(luout,fmtstr) 
        write(luout,'(x,a)') 'INPUT ERROR: unexpected keyword '
      end if
  
      return
      end subroutine

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
      
      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' context = "',trim(context),'"'
        if (forward) write(lulog,*) ' forward search'
        if (.not.forward) write(lulog,*) ' backward search'
      end if 

      finnode => null()
      current => tree_root

      ipst = 1
      len = len_trim(context)
      if (len.eq.0) finnode => current

      subk_loop: do while(ipst.le.len)

        found=.False.
        ipnd = index(context(ipst:),".")+ipst-2
        if (ipnd.lt.ipst) ipnd = len

        if (ntest.ge.100) then
          write(lulog,*) ' current subkeyword: "',context(ipst:ipnd),'"'
        end if


        ! get all keywords one level down
        if (.not. hasChildNodes(current)) exit subk_loop
        nodes_list=>getChildNodes(current)
        if (forward)then
           node_loop: do ii=0,getLength(nodes_list)-1 !DAMN YOU C!!!
           current=> item(nodes_list,ii)
           ! compare tags to filter keywords and name for specific keyword

           if (getNodeName(current) 
     &          .eq. key_tag .and. 
     &          getAttribute(current, atr_name)
     &          .eq. context(ipst:ipnd)) then 
              current=>item(nodes_list,ii)
              found=.True.
              exit node_loop
           end if
           
           end do node_loop
        else
           nodes_loop: do ii=getLength(nodes_list)-1,0,-1
           current=> item(nodes_list,ii)
           if (getNodeName(current) 
     &          .eq. key_tag .and. 
     &          getAttribute(current, atr_name)
     &          .eq. context(ipst:ipnd)) then 
              current=>item(nodes_list,ii)
              found=.True.
              exit nodes_loop
           end if

           end do nodes_loop
        end if

        if (.not.found) then
           exit subk_loop
        else if (found .and. ipnd.eq.len )then
           finnode => current
           exit subk_loop
        end if

        ipst=ipnd+2
      end do subk_loop

      if (ntest.ge.100) then
        if (associated(finnode)) write(lulog,*) 'success'
        if (.not.associated(finnode)) write(lulog,*) 'no success'
      end if
      end subroutine

*----------------------------------------------------------------------*
!>    delivers the next keyword in a depth first search
!!
!!    prefilters with a specific tag.
!!    @param[inout] nxtnode on entry: node where the search is continued
!!                          on exit : nextnode  (or null if key_root would be next node)
!!    @param[in] key string that contains the tag we filter for.
*----------------------------------------------------------------------*
      subroutine dsearch_next_key(nxtnode, tag,level)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      integer,parameter::
     &     ntest= 00
      character(len=16),parameter ::
     &     i_am="dsearch_next_key"

      character,intent(in) ::
     &     tag*(*)
      integer,intent(inout),optional ::
     &     level

      integer::
     &     dlevel !levelchange
      type(Node), pointer, intent(inout)::
     &     nxtnode
      type(Node), pointer::
     &     current,siblnode
      

      current => nxtnode

      if (ntest .ge. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' looking for tag = "',trim(tag),'"'
         write(lulog,*) ' association status of nxtnode:',
     &        associated(nxtnode)
      end if 

      dlevel=0
      
      main_loop: do
         siblnode=>getNextSibling(current)

         if (hasChildNodes(current))then
            current=> getFirstChild(current)
            dlevel=dlevel+1

         else if (associated(siblnode))then
            current=> siblnode
         else
 
            up_loop: do while(getNodeName(current).ne.key_root_tag)
               current => getParentNode(current)
               dlevel=dlevel-1

               siblnode=>getNextSibling(current)
               if (associated(siblnode))then 
                  current=> siblnode
                  exit up_loop
               end if 
            end do up_loop
         end if 

         if (getNodeName(current).eq. trim(tag))then
            nxtnode=>current
            exit main_loop
         else if (getNodeName(current).eq.key_root_tag) then
            nxtnode=>null()
            exit main_loop
         end if
      end do main_loop
      if (present(level))level=level+dlevel
      end subroutine 
*----------------------------------------------------------------------*
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


*----------------------------------------------------------------------*
!>    function to convert any input for a logical argument to a format readable by rts     
!!    
!!    it converts the string to logical(s) and back
!!    @param valstr string to be converted
!!    @param dim maximum dimension of the input
!!    @param[out] ierr error code
*----------------------------------------------------------------------*
      function conv_logical_inp(valstr,dim,ierr) 
*----------------------------------------------------------------------*
      implicit none
      integer,intent(in)::
     &     dim
      character(len=dim*6):: ! worst case dim*'false,' (=dim*6-1) 
     &     conv_logical_inp
      character(len=16),parameter::
     &     i_am="input_to_logical"
      integer,parameter ::
     &     ntest=00
      integer,parameter::
     &     ERR_TO_MANY_ELEMENTS=1,
     &     ERR_UNCONVERTIBLE=2

      character(len=*),intent(in)::
     &     valstr
      integer,intent(out)::
     &     ierr

      logical,Dimension(dim)::
     &     larr
      
      integer ::
     &     ii, 
     &     ipst, ipnd
      
      ierr=0
      conv_logical_inp=" "
      ipst=1
      ipnd=index(valstr(ipst:),",")+ipst-1

      if (ipnd.gt.ipst)then
         ii=0
         do while (ipnd.gt.ipst)
            ii=ii+1
            if (ii.gt.dim)then 
               ierr=ERR_TO_MANY_ELEMENTS
               return
            end if 

            larr(ii)=str_to_logical(valstr(ipst:ipnd),ierr)

            if (ierr.gt.0) return 

            ipst=ipnd+2
            ipnd=index(valstr(ipst:),",")+ipst-1
         end do
         conv_logical_inp=str(larr(1:ii))
      else
         conv_logical_inp=str(str_to_logical(valstr,ierr))
      end if
      end function

*----------------------------------------------------------------------*
!> converts a single element of user input to a logical      
*----------------------------------------------------------------------*
      logical function str_to_logical(instr,ierr) 
*----------------------------------------------------------------------*
      implicit none
      integer,parameter::
     &     ERR_UNCONVERTIBLE=2

      character,intent(in)::
     &     instr*(*)
      integer,intent(out)::
     &     ierr
      ierr=0

      select case (trim(instr))
         case("T","t")
            str_to_logical=.true.
         case("F","f")
            str_to_logical=.false.
         case default
            str_to_logical=.false.
            ierr=ERR_UNCONVERTIBLE
      end select
      end function 
      end module
