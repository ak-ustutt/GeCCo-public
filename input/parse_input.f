      module parse_input
      use keyword_trees
      use FoX_common, only: rts,str
      include 'par_vtypes.h'
      integer,parameter::
     &     maxvalstr_len=256


      character(len=*),parameter::
     &     calculate_name="calculate"   ! important in input postprocessing
      contains

!=======================================================================
! output subroutines
!=======================================================================

*----------------------------------------------------------------------*
!>    wrapper for keyword_list for the input
*----------------------------------------------------------------------*
      subroutine inp_show(unit)
*----------------------------------------------------------------------*
      integer,intent(in)::
     &     unit

      type(tree_t)::
     &     input
      input = fetch_input_keyword_tree()
      call keyword_list(unit,input,show_args=.True.)
      end subroutine

*----------------------------------------------------------------------*
!>    wrapper for keyword_list for the input
*----------------------------------------------------------------------*
      subroutine reg_show(unit)
*----------------------------------------------------------------------*
      integer,intent(in)::
     &     unit
      type(tree_t),pointer::
     &     keytree
      keytree = fetch_input_keyword_tree()
      call keyword_list(unit,keytree,show_args=.True.)
      end subroutine



*----------------------------------------------------------------------*
!>    prints the keyword_tree below a certain keyword:
!!
!!
!!    @param luwrt ouput unit
!!    @param tree_root starting keyword
!!    @param n_descent maximum number of levels to descent
!!    @param show_args if arguments should be shown
!!    @param show_status if status should be shown
*----------------------------------------------------------------------*
      subroutine keyword_list(luwrt,tree,
     &     n_descent,show_args)
*----------------------------------------------------------------------*
      use FoX_common, only:str,rts
      implicit none
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
      type(tree_t), intent(inout)::
     &     tree

      integer,optional,intent(in)::  ! not implemented
     &     n_descent
      logical, optional, intent(in):: 
     &     show_args

      logical:: 
     &     args_vis,status_vis
      integer::
     &     level, 
     &      type, dim,
     &     ii, dummy
      character::
     &     status

      type(Node),pointer::
     &     curkey,curarg
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)

      end if

      args_vis=.false.
      if (present(show_args)) args_vis=show_args

      if (ntest.ge.100) then
         write (lulog,*) "show_status?",status_vis
         write (lulog,*) "show_arguments?",args_vis
      end if


      curarg=>null()
      curkey=> tree_iterate(tree)
      
      key_loop: do 
         call print_keyword(luwrt,curkey,getLevel(tree))
         if ( args_vis)then
            curarg=> elem_getFirstChild(curkey,arg_tag)
            arg_loop: do while (associated(curarg))
               call print_argument(luwrt,curarg, getLevel(tree))
               curarg=> elem_getNextSibling(curarg, arg_tag)
            end do arg_loop 
         end if
         curkey=> tree_iterate(tree)
         if (.not.associated(curkey)) exit key_loop
         do while (getNodeName(curkey).ne. key_tag)
            curkey=> tree_iterate(tree)
            if (.not.associated(curkey)) exit key_loop
         end do 
      end do key_loop
      return
      contains
*----------------------------------------------------------------------*
!>    prints a keyword 
!!
!!   @param luwrt output unit
!!   @param arg argument node
!!   @param status_vis if status should be printed as well
*----------------------------------------------------------------------*
      subroutine print_keyword(luwrt,keywd,level)
*----------------------------------------------------------------------*
      integer, parameter ::
     &     ntest=00
      character(len=*),parameter ::
     &     i_am="print_keyword"
      integer, intent(in)::
     &     luwrt
      type(Node),pointer, intent(in)::
     &     keywd
      character(len=64)::
     &     fmtstr
      integer,intent(in)::
     &     level

      if (hasAttribute(keywd,atr_stat))then 
         status=getAttribute(keywd,atr_stat)
         if (ntest.ge.100)then 
            write (lulog,*) "status:",getAttribute(keywd,atr_stat)
         end if 
         if (status.eq.status_active) then
            write(fmtstr,'("(""A"",",i3,"x,a)")') 2*level+1
         else
            write(fmtstr,'("(""I"",",i3,"x,a)")') 2*level+1
         end if
      else 
         write(fmtstr,'("("">"",",i3,"x,a)")') 2*level+1
      end if
      write (luwrt,fmtstr)getAttribute(curkey,atr_name)
      end subroutine print_keyword
*----------------------------------------------------------------------*
!>    prints a keyword 
!!
!!   @param luwrt output unit
!!   @param arg argument node
!!   @param status_vis if status should be printed as well
*----------------------------------------------------------------------*
      subroutine print_argument(luwrt,arg,level)
*----------------------------------------------------------------------*
      integer, parameter ::
     &     ntest=00
      character(len=*),parameter ::
     &     i_am="print_argument"
      integer, intent(in)::
     &     luwrt
      integer,intent(in)::
     &     level

      type(Node),pointer, intent(in)::
     &     arg
      character(len=64)::
     &     fmtstr
      integer::
     &     type, dim
      call rts(getAttribute(arg,atr_kind),type)
      call rts(getAttribute(arg,atr_len),dim)
      write(fmtstr,'("(x,",i3,"x,a,x,i2,x,a)")') 2*level+4
      write(luwrt,fmtstr) getAttribute(curarg,atr_name)//" "
     &     //trim(type_array(type))//" of len",dim,": "// 
     &     getAttribute(curarg,atr_val)
      end subroutine
      end subroutine










!======================================================================!
!   postprocess routines
!======================================================================!



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
      type(tree_t)::
     &     inputtree

      inputtree=fetch_input_keyword_tree()
      call tree_toggle_status_(inputtree, one_more)
      end subroutine 






*----------------------------------------------------------------------*
!!    toggles the status of the toplevel keywords (keywords below the root of tree) 
!!
!!   looks for the first "calculate" block with neither active nor inactive status
!!    if none is found, advances over all context
!!   for any given toplevel keyword(including calculate) only the last block 
!!       before and including this "calculate" is set active
!!   all other previous contexts are set inactive.
!!   @param tree the tree object
!!   @param history pointer to the last active history file (may be null())
!!   @param one_more true if active blocks were found
*----------------------------------------------------------------------*
      subroutine tree_toggle_status_(tree,one_more)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      integer,parameter::
     &     ntest= 00
      character(len=*),parameter ::
     &     i_am="tree_toggle_status"
      type(tree_t),intent(in)::
     &     tree
      logical, intent(inout)::
     &     one_more
      type(Node),pointer::
     &     calculate_ptr
      type(Node),pointer::
     &     current,nxtnode
      integer::
     &     i
      one_more=.true.
      current=> getRoot(tree)
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) " root: ",getAttribute(current,atr_name)
      end if 
      current=> elem_getFirstChild(tree%root)
      
      if (.not.associated(current))
     &     call quit(1,i_am,'not a single keyword given?')
      
      do while (associated(nxtnode))
         current=>nxtnode
         if ( ( getAttribute(current, atr_name) .eq. calculate_name) 
     &        .and. (.not. hasAttribute(current, atr_stat) ) )
     &        exit
         nxtnode=> elem_getNextSibling(current)
      end do

      if (.not. associated(nxtnode) .and. 
     &     hasAttribute(current, atr_stat))then
         one_more=.false.
         return
      end if

      calculate_ptr=> current !points to a calculate or the last toplevel keyword in input

      set_active_loop: do while(associated(current))
         call set_status(current,status_active)
         call unset_previous_keywords(current)
         current=> elem_getPreviousSibling(current)
         backtrack: do while(associated(current))
            if (hasAttribute(current, atr_stat))then
               if (getAttribute(current,atr_stat).eq.status_active)
     &              exit set_active_loop
            else
               exit backtrack
            end if
            current=> elem_getPreviousSibling(current)
         end do backtrack
      end do set_active_loop
      
      end subroutine 

*----------------------------------------------------------------------*
!>     set all previous keywords with same key as present keyword to
!!     inactive (+ all sub-levels)
*----------------------------------------------------------------------*
      subroutine unset_previous_keywords(keywd)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=*),parameter ::
     &     i_am="unset_previous_keywords"
      integer, parameter ::
     &     ntest=00
      type(Node), pointer,intent(in) ::
     &     keywd     
      type(Node), pointer ::
     &     current
      character(len=max_name_len) ::
     &     name
      integer :: i
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of keyword:"
     &        ,associated(keywd)
      end if 

      current=>elem_getPreviousSibling(keywd, key_tag)
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
     &        call set_status(current,status_inactive)
         current =>elem_getPreviousSibling(current,key_tag)
         if (associated(current) )then
            continue
         else
            exit 
         end if
      end do 
      return 
      end subroutine

*----------------------------------------------------------------------*
!>    sets the given keyword  to a given status
!!
!!   @param keywd pointer to a keyword
!!   @param status integer with the status we set it to 
*----------------------------------------------------------------------*
      subroutine set_status(keywd,status)
*----------------------------------------------------------------------*
      implicit none
      type(Node), pointer, intent(in)::
     &     keywd
      character,intent(in)::
     &     status

      call setAttribute(keywd,atr_stat,status)

      end subroutine











*----------------------------------------------------------------------*
!>    for a given argument determines the actual length (and type)
!!    
!!    @TODO improve string handling
*----------------------------------------------------------------------*
      subroutine get_argument_dimension_core(curarg,num,type,succ)
*----------------------------------------------------------------------*
      use Fox_common, only: rts
      include 'par_vtypes.h'
      include 'stdunit.h'
      integer,parameter::
     &     ntest=00
      character(len=*),parameter ::
     &     i_am="get_argument_dimension_core"

      type(Node), pointer,intent(inout) ::
     &     curarg
      integer , intent(out)::
     &     num,type
      logical , intent(out)::
     &     succ
     
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
         call test_conv_argument(getAttribute(curarg,atr_val),dim_tot,
     &       type,succ=succ,dim=num)
      end if 


      if (ntest.ge.100) then 
        write(lulog,'("length:",i3)')num 
      end if
      return 
      end subroutine 



*----------------------------------------------------------------------*
!>    tests if a string can successfully be converted into arguments of a given type
!!
!!    as a side effect the subroutine also determines the number of actually existing elements
!!    @param[in] valstr the string to be converted
!!    @param[in] dim_tot maximal number of elements to convert
!!    @param[in] type type as vtype 
!!    @param[out] succ indicates a successful conversion
!!    @param[out] dim actual number of converted elements
!!    Note: there is a FoXdom utility function, which does about the same
*----------------------------------------------------------------------*
      subroutine test_conv_argument(valstr,dim_tot,type,succ,dim)
*----------------------------------------------------------------------*
      integer, parameter ::
     &     ntest=00
      character(len=17),parameter ::
     &     i_am="test_argument_conversion"

      character(len=*),intent(in)::
     &     valstr
      integer,intent(in)::
     &     dim_tot,type
      logical,intent(out),optional::
     &     succ
      integer,intent(out),optional::
     &     dim

      logical::
     &     insucc
      integer::
     &     indim, ex

      logical,allocatable::
     &     larr(:)
      integer,allocatable::
     &     iarr(:)
      real(8),allocatable::
     &     xarr(:)
      
      insucc=.false.
      indim=1


      select case(type)
      case (vtyp_log)
         allocate(larr(dim_tot))
         call rts(valstr,larr,iostat=ex,num=indim) 
         if (ex.le. 0) insucc = .true.
      case (vtyp_int)
         allocate(iarr(dim_tot))
         call rts(valstr,iarr
     &        ,iostat=ex,num=indim) 
         if (ex.le. 0) insucc = .true.
      case (vtyp_rl8)
         allocate(xarr(dim_tot))
         call rts(valstr,xarr
     &        ,iostat=ex,num=indim) 
         if (ex.le. 0) insucc = .true.
      case (vtyp_str)
         num=dim_tot
         insucc = .true.
      end select
      if(present(succ)) succ=insucc
      if(present(dim)) dim=indim
      end subroutine
*----------------------------------------------------------------------*
!!   wrapper to uncouple keyword_parse_ from module variables
!!
!!  @param  luin input unit
*----------------------------------------------------------------------*
      subroutine inp_parse(luin)
*----------------------------------------------------------------------*
      implicit none
      integer, intent(in) ::
     &     luin
      type(tree_t)::
     &     keytree, input

      keytree= fetch_registry_keyword_tree()
      input= fetch_input_keyword_tree()
      call keyword_parse_(luin,keytree,input)
      end subroutine inp_parse

*----------------------------------------------------------------------*
!>     parse the keywords on unit luin
!!     the unit should be a formatted, sequential file, positioned
!!     at the place where the parser should start
*----------------------------------------------------------------------*
      subroutine keyword_parse_(luin,keytree,input)
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

      type(tree_t),intent(inout)::
     &     keytree,input               ! root element of the preset keyword tree
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
     &     curkey,nxtkey,dummykey
      type(Node), pointer ::
     &     curarg
      character(len=maxlen) ::
     &     context
      integer ::
     &     allowed_delim(n_delim), n_allowed_delim
      integer ::
     &     ipst, ipnd, itest, lenline, ierr, dummy

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if

      context = " "
      curkey => getRoot(keytree)

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
            
!     is it an argument key?

            
            

            dummy=1
            curarg=>getSubNode(curkey,line(ipst:ipnd),
     &           arg_tag,latest=.false.,icount=dummy )


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
                     call create_node(input,curarg,line(ipst+1:ipnd-1))
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
                     call create_node( input, curarg,line(ipst:ipnd))
                  end if 

               else
                  ierr = ierr+1
                  call error_eol(line,ipst)
                  allowed_delim(1:n_allowed_after_key) =
     &                 allowed_after_key(1:n_allowed_after_key)
                  n_allowed_delim = n_allowed_after_key
                  exit line_loop
               end if
            else      
               nxtkey=> tree_goback_to_element(keytree,line(ipst:ipnd),
     &           key_tag )

               if (.not.associated(nxtkey)) then
                  ierr = ierr+1
                  call error_keywd(line,ipst)
               else
                  curkey => nxtkey
                  call create_node(input, curkey," ")
                  
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
!>    writes the line and generates pointer to the position of the error
*----------------------------------------------------------------------*
      subroutine error_pointer(ipos,line,msg)
*----------------------------------------------------------------------*
      implicit none
      include "stdunit.h"
      character(len=*),parameter::
     &     unit="O"

      character(len=*), intent(in) ::
     &     line,msg
      integer, intent(in) ::
     &     ipos

      character(len=80) ::
     &     fmtstr, outstr
      call print_out(' ',unit)
      write(outstr,'(x,a)') line
      call print_out(outstr,unit)
      if (ipos.gt.1)then
         write(fmtstr,'("(x,""",a,""",""^"")")') repeat("-",abs(ipos)-1)
      else
         write(fmtstr,'("(x,""^"")")')
      end if
      write(luout,*)fmtstr
      write(outstr,fmtstr)
      call print_out(outstr,unit)
      write(outstr,'(x,"INPUT ERROR: ",a)') msg
      call print_out(outstr,unit)
      end subroutine


*----------------------------------------------------------------------*
      subroutine error_delim(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character(len=80) ::
     &     fmtstr, msg_str

      write(fmtstr,'("(x,a,",i3,"(a1,x))")') n_allowed_delim

      write(msg_str,fmtstr) 'unexpected delimiter, '//
     &     'expected one of ',
     &     delimiter(allowed_delim(1:n_allowed_delim))

      call error_pointer(ipos,trim(line), trim(msg_str))

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

      call error_pointer(ipos,trim(line), 'unexpected EOL')
      

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

      call error_pointer(ipos,trim(line), 'missing )')

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

      call error_pointer(ipos,trim(line),'missing =')
      return
      end subroutine


*----------------------------------------------------------------------*
      subroutine error_keywd(line,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character(len=*), intent(in) ::
     &     line
      integer, intent(in) ::
     &     ipos

      call error_pointer(ipos,trim(line), 'unexpected keyword')
  
      return
      end subroutine

      end subroutine 

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
      subroutine create_node(input, template,value)
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
      type(tree_t),intent(inout)::
     &     input
      type(node), pointer ::
     &     template
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
      character(len=maxvalstr_len)::
     &     invalue
      logical ::
     &     convertible

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write (lulog,*) "testing association of template:"
     &        ,associated(template)
         write (lulog,*) "tag, name of template:",
     &        getNodeName(template),getAttribute(template,atr_name)
      end if

      
      if (.not. associated(template)) 
     &     call quit(1,i_am,"template not set")
      
      new_elem=>tree_create_new_element( input, template)

      call setAttribute(new_elem,atr_name,
     &     getAttribute(template,atr_name))

      if (getNodeName(new_elem) .eq. arg_tag)then 

         call setAttribute(new_elem,atr_len
     &        ,getAttribute(template,atr_len))
         call setAttribute(new_elem,atr_kind,
     &        getAttribute(template,atr_kind))
         call rts(getAttribute(template,atr_kind),kind)
         call rts(getAttribute(template,atr_len),dim)

         invalue=trf_in_val(value,dim,kind,ierr)
         call test_conv_argument(invalue,dim,kind,succ=convertible)
         if (.not. convertible)then 
            call quit(0,i_am,"cannot convert given values:"//value
     &           //"for argument"//getAttribute(template,atr_name))
         end if
!! @TODO: implement a check for invalid input
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
         call setAttribute(new_elem, atr_val,trim(invalue))
      end if
      return 
      contains 
      function trf_in_val(value,dim,kind,ierr) result(ret)

      character(len=maxvalstr_len)::
     &     ret
      character(len=*),intent(in)::
     &     value
      integer,intent(in)::
     &     dim, kind
      integer,intent(out)::
     &     ierr
      ierr=0
      select case(kind)
      case (vtyp_log)
         ret=conv_logical_inp(value,dim,ierr)
      case default
         ret=value
      end select
      return
      end function
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
      return
      end function 

      end subroutine

      end module  
