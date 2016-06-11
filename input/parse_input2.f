      module parse_input2
      use FoX_dom,only: Node,Nodelist,DOMException,
     &     parseFile,parseString,getExceptionCode,getElementsByTagName,
     &     hasChildNodes,getChildNodes,getParentNode,getFirstChild,
     &     getLength,getNodeName,getAttribute,importNode,item,
     &     appendChild,setAttribute
      implicit none
      include 'write_styles.h'

      integer,parameter::
     &     file_loc_len=16
      character(len=file_loc_len),parameter ::
     &     rel_file_loc="/keyword_registry"
      !> @TODO make this usable on systems where / is not directory separator
      character,parameter::     !attribute names
     &     atr_name*3="key",
     &     atr_kind*4="kind",
     &     atr_len*3="len",
     &     atr_val*3="val"

      

      character,parameter::      !tags 
     &     key_root_tag*8="key_root",
     &     key_tag*7="keyword",
     &     arg_tag*8="argument"



      type(Node), pointer :: 
     &     preset_doc,         !root of the tree that represents the keyword registry
     &     key_root,            !root element for keywords in registry
     &     input_doc           !root element of input; corresponds to key_root
      



      contains 
*----------------------------------------------------------------------*
!>    returns name of the keyword_file
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
!>    parses the keyword_file 
!!    @TODO think if that should be done here or in set_keywords
*----------------------------------------------------------------------*
      subroutine keyword_init_()
*----------------------------------------------------------------------*
      implicit none
      character(len=13),parameter ::
     &     i_am="keyword_init_"
      type(DOMException)::
     &     ex
      type(Node),pointer::     tmpnode
      type(NodeList),pointer::
     &     nodes_list
      integer ::
     &     i
      
      preset_doc => parseFile(trim(get_keyword_file()), ex=ex)
      !> @TODO more special error handling 
      if (getExceptionCode(ex) .ne. 0)
     &     call quit(1,i_am,"could not read keyword registry")
      nodes_list=>getElementsByTagName(preset_doc, key_root_tag)
      i=1
      tmpnode=> item(nodes_list, i)!  assumes there is only on of these elements
      input_doc=> parseString(
     &     "<"//key_root_tag//"></"//key_root_tag//">")
      end subroutine 


*----------------------------------------------------------------------*
      subroutine find_active_node(tree,finnode,context,icount)
*----------------------------------------------------------------------*

      implicit none

      type(Node), pointer ::
     &     tree
      type(Node), pointer ::
     &     finnode
      character ::
     &     context*(*)
      integer ::
     &     icount

      character ::
     &     curcontext*1024     
      type(Node), pointer ::
     &     current
      integer ::
     &     jcount, status

      if (.not.associated(tree%down_h))
     &     call quit(1,'find_active_node','invalid keyword tree')

      current => tree%down_h

      finnode => null()

      jcount = 0
      key_loop: do 
        call rts(getAttribute(current,atr_stat),status)
        if (status.gt.0) then
          call keyword_get_context(curcontext,current)
          if (trim(context).eq.trim(curcontext)) jcount = jcount+1 
          if (icount.eq.jcount) then
            finnode => current
            exit key_loop
         end if
        end if

        if (status.gt.0.and.hasChildNodes(current)) then
          ! go down
          current => getFirstChild(current)
        else if (associated(getNextSibling(current))) then
          ! else stay within level
          current => getNextSibling(current)
        else
          ! else find an upper level, where a next
          ! node exists:
          up_loop: do
             if (getNodeName(getParentNode(current))
     &         .ne. key_root_tag) then
                current => getParentNode(current)
                if (associated(getNextSibling(current))) then
                   current => getNextSibling(current)
                   exit up_loop
                end if
             else
                exit key_loop
             end if
          end do up_loop
       end if        
      end do key_loop

      return
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

      
      if (ntest .gt. 100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) ' context = "',trim(context),'"'
        if (forward) write(lulog,*) 'forward search'
        if (.not.forward) write(lulog,*) 'backward search'
      end if 
      forward = .true.
      if (present(latest)) forward = .not.latest

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
           node_loop: do ii=1,getLength(nodes_list) 

           ! compare tags to filter keywords and name for specific keyword

           if (getNodeName(item(nodes_list,ii)) 
     &          .eq. key_tag .and. 
     &          getAttribute(item(nodes_list,ii), atr_name)
     &          .eq. context(ipst:ipnd)) then 
              current=>item(nodes_list,ii)
              found=.True.
              exit node_loop
           end if
           
           end do node_loop
        else
           nodes_loop: do ii=getLength(nodes_list),1,-1

           if (getNodeName(item(nodes_list,ii)) 
     &          .eq. key_tag .and. 
     &          getAttribute(item(nodes_list,ii), atr_name)
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
!>    
*----------------------------------------------------------------------*
*     look whether keyword cur_key hosts an argument with key "key"
*     and return the corresponding node
*----------------------------------------------------------------------*
      subroutine arg_node(arg,cur_key,key)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'

      character(len=9),parameter::
     &     i_am="next_node"
      integer, parameter ::
     &     ntest = 00

      type(Node), pointer ::
     &     cur_key
      character, intent(in) ::
     &     key*(*)
      type(Node), pointer ::
     &     arg
      type(NodeList),pointer ::
     &     nodes_list
      integer ::
     &     ii

      arg => null()
      if (hasChildNodes(cur_key))then
         nodes_list=> getChildNodes(cur_key)
         do ii=1,getLength(nodes_list)
            if (getNodeName(item(nodes_list,ii))
     &           .eq. arg_tag .and. 
     &           getAttribute(item(nodes_list,ii), atr_name)
     &           .eq. key ) then
               arg => item(nodes_list,ii)
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
     &     ntest = 00

      type(node), pointer ::
     &     cur_node
      type(node), pointer ::
     &     nxt_node
      character, intent(in) ::
     &     key*(*)

      type(node), pointer ::
     &     current
      type(NodeList), pointer ::
     &     keylist
      integer ::
     &     ii
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
!         write(lulog,*) ' start = "',trim(get),'"'
!         write(lulog,*) ' search for "',trim(key),'"'
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
         do ii=1,getLength(keylist)

            if (getNodeName(item(keylist,ii)) 
     &           .eq. key_tag .and. 
     &           getAttribute(item(keylist,ii), atr_name)
     &           .eq. key ) then 
               current=>item(keylist,ii)
               nxt_node=> item(keylist,ii)
               exit node_loop
            end if

         end do

         if (getNodeName(current).ne. getNodeName(key_root))then
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
!!   overrides the argument value with value. for keywords value is irrelevant
!!   @param[inout] doc document that will own the newly created node. 
!!   @param[in] template template for the newly created node
!!   @param[in] value string which can be converted to input
!!   @TODO error checking
*----------------------------------------------------------------------*
      subroutine create_node(doc,template,value)
*----------------------------------------------------------------------*
      implicit none 
      character(len=11),parameter::
     &     i_am="create_node"
      integer,parameter ::
     &     ntest=00

      integer, parameter ::
     &     maxlen  = 256

      type(node), pointer ::
     &     doc,template
      character(len=*),intent(in)::
     &     value



      type(node),pointer ::
     &     new_elem,new_parent
      character(len=maxlen) ::
     &     context

      if (.not. associated(doc)) 
     &     call quit(1,i_am,"doc not set")
      
      if (.not. associated(template)) 
     &     call quit(1,i_am,"template not set")
      
      new_elem=>importNode(doc, template, .False.)


      call  keyword_get_context(context,template)

!> @TODO this assumes that key_root is the first child of doc      
      call find_node(getFirstChild(doc),new_parent,context)

      new_elem=>appendChild(new_parent, new_elem)
      if (getNodeName(new_elem) .eq. arg_tag) 
     &     call setAttribute(new_elem,atr_val, value)
         

      end subroutine

*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine keyword_get_context(curcontext,current)
*----------------------------------------------------------------------*
      character(len=19),parameter::
     &     i_am="keyword_get_context"
      integer,parameter::
     &     ntest=00
      character ::
     &     curcontext*(*)
      type(node), pointer ::
     &     current 
      
      type(node), pointer ::
     &     internal
 
      
      internal=> getParentNode(current)
      do while (getNodeName(internal).ne. key_root_tag)
         curcontext=getAttribute(internal,atr_name)//
     &   "."//trim(curcontext)
         internal=> getParentNode(internal)
      end do

      end subroutine

*----------------------------------------------------------------------*
!Output subroutine; to be implemented later
*----------------------------------------------------------------------*
!      subroutine keyword_list(lulog,tree_root,
!     &     context,n_descent,show_args)


*----------------------------------------------------------------------*
*     parse the keywords on unit luin
*     the unit should be a formatted, sequential file, positioned
*     at the place where the parser should start
*----------------------------------------------------------------------*
      subroutine keyword_parse(luin)
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



      context = " "
      curkey => key_root

      allowed_delim(1:n_allowed_start)=allowed_start(1:n_allowed_start)
      n_allowed_delim = n_allowed_start


      ierr = 0
      file_loop: do
         read(luin, '(a)', end=100, err=200) line
         ipst = 1
         
         call clean_line(line,delimiter,n_delim)
         
         lenline=len_trim(line)

         !empty line?

         if (lenline.le.0) cycle file_loop 

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
            itest = next_delim(line(ipst:ipst),
     &           allowed_delim,n_allowed_delim)
            
            ipnd = abs(itest)+ipst-2
            
            if (itest.le.0) then
               ierr = ierr+1
               call error_delim(line,ipnd+1)
            end if
            
! is it an argument key?
            call arg_node(curarg,curkey,line(ipst:ipnd))
            
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
                  end if 
                  call create_node(input_doc,curarg,line(ipst:ipnd))
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
               if (.not.associated(nxtkey)) then
                  ierr = ierr+1
                  call error_keywd(line,ipst,curkey)
               else
                  curkey => nxtkey
                   call create_node(input_doc,curkey," ")
!     add node to keyword history
               end if
            end if
            
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
!!    accesses global var delimiter and n_delim
!!    @param str string to be searched
!!    @param delim(n_delim) array of allowed delimiters
!!    @param n_delim
!!    return value:
!!      position of next delimiter
!!      len+1 if line end is delimiter
!!      negative value, if delimiter is not allowed in context
!!    
*----------------------------------------------------------------------*
      integer function next_delim(str,delim,n_delim)
*----------------------------------------------------------------------*
      implicit none

      character, intent(in) ::
     &     str*(*)
      integer,intent(in)::
     &     n_delim
      integer,dimension(n_allowed_delim),intent(in)::
     &     delim

      logical ::
     &     ok
      integer ::
     &     ipos, jpos, len, idelim, jdelim

      len = len_trim(str)

      if (len.eq.0) then
        next_delim = 0
        return
      end if

      ipos = len+1
      idelim = 0
      do jdelim = 1, n_delim
        jpos = index(str,delimiter(jdelim))
        if (jpos.gt.0.and.jpos.lt.ipos) then
          ipos = jpos
          idelim = jdelim
        end if
        if (ipos.eq.1) exit
      end do

      ok = .true.
      if (idelim.gt.0) then
        ok = .false.
        do jdelim = 1, n_allowed_delim
          if (idelim.eq.allowed_delim(jdelim)) then
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


      end subroutine keyword_parse



      end module
