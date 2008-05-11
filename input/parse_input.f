*----------------------------------------------------------------------*
      module parse_input
*----------------------------------------------------------------------*
*     routines to set up a keyword-list and for parsing
*     user input into this list
*----------------------------------------------------------------------*

      implicit none

      ! parameters
      integer, parameter ::
     &     vtyp_log = 1,    ! input variable types
     &     vtyp_int = 2,
     &     vtyp_rl8 = 4,
     &     vtyp_str = 8,
     &     sttyp_undef = 0,
     &     sttyp_set   = 1,
     &     sttyp_auto  = 2
      integer, parameter ::
     &     lenkey = 16

      type value
        integer ::
     &     type, len
        logical, allocatable ::
     &     lval(:)
        integer, allocatable ::
     &     ival(:)
        real(8), allocatable ::
     &     xval(:)
        character, allocatable ::
     &     cval(:)
      end type value

      type argument
        character ::
     &     key*(lenkey)
        type(argument), pointer ::
     &       prev, next
        type(value) ::
     &       val
        integer ::
     &       status
      end type argument

      type keyword
        character ::
     &     key*(lenkey)
        type(keyword), pointer ::
     &       prev, next, up, down_h, down_t
        type(argument), pointer ::
     &       arg_h, arg_t
        logical ::
     &       required
        integer ::
     &       status
      end type keyword

      type(keyword), pointer ::
     &     keyword_root
      type(keyword), pointer ::
     &     keyword_history
      type(keyword), pointer ::
     &     history_pointer
      integer ::
     &     keyword_status = -1

      contains

*----------------------------------------------------------------------*
      subroutine keyword_init()
*----------------------------------------------------------------------*
*     inititialization routine
*----------------------------------------------------------------------*

      implicit none

      allocate(keyword_root,keyword_history)
      keyword_root%key="gecco keywords"
      nullify(keyword_root%prev)
      nullify(keyword_root%next)
      nullify(keyword_root%up)
      nullify(keyword_root%down_h)
      nullify(keyword_root%down_t)
      nullify(keyword_root%arg_h)
      nullify(keyword_root%arg_t)
      keyword_history%key="keyword hist"
      nullify(keyword_history%prev)
      nullify(keyword_history%next)
      nullify(keyword_history%up)
      nullify(keyword_history%down_h)
      nullify(keyword_history%down_t)
      nullify(keyword_history%arg_h)
      nullify(keyword_history%arg_t)

      keyword_status = 0

      return
      end subroutine

c gives internal compiler error:
c      pure function find_node(context)
c      implicit none
c      type(keyword), pointer ::
c     &     find_node
c      character, intent(in) ::
c     &     context*(*)
c
c      type(keyword), pointer ::
c     &     current
c
c      find_node => null()
c
c      return
c      end function

*----------------------------------------------------------------------*
      subroutine find_node(tree_root,node,context,latest)
*----------------------------------------------------------------------*
*     context is a string as e.g. "key.subkey.subsubkey"
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      type(keyword), pointer ::
     &     tree_root
      type(keyword), pointer ::
     &     node
      character, intent(in) ::
     &     context*(*)
      logical, intent(in), optional ::
     &     latest

      logical ::
     &     forward
      integer ::
     &     ipst, ipnd, len
      type(keyword), pointer ::
     &     current

      forward = .true.
      if (present(latest)) forward = .not.latest

      if (ntest.ge.100) then
        write(luout,*) '-----------'
        write(luout,*) ' find_node'
        write(luout,*) '-----------'
        write(luout,*) ' context = "',trim(context),'"'
        if (forward) write(luout,*) 'forward search'
        if (.not.forward) write(luout,*) 'backward search'
      end if

      node => null()

      ! we start at the root level
      current => tree_root

      ipst = 1
      len = len_trim(context)
      if (len.eq.0) node => current ! no context = root

      subk_loop: do while(ipst.le.len)
        ipnd = index(context(ipst:),".")+ipst-2
        if (ipnd.lt.ipst) ipnd = len

        if (ntest.ge.100) then
          write(luout,*) ' current subkeyword: "',context(ipst:ipnd),'"'
        end if

        if (forward) then
          ! start search at head of next sublevel
          if (.not.associated(current%down_h)) exit subk_loop
          current => current%down_h
        else
          ! start search at tail of next sublevel
          if (.not.associated(current%down_t)) exit subk_loop
          current => current%down_t
        end if

        node_loop: do

          if (trim(current%key).eq.context(ipst:ipnd)) then
            ! last part of context? return this node
            if (ipnd.eq.len) node => current
            exit node_loop
          end if
          
          if (forward) then
            ! look at next node
            if (associated(current%next)) then
              current => current%next
            else
              ! last element => search unsuccessful
              exit subk_loop
            end if
          else
            ! look at previous node
            if (associated(current%prev)) then
              current => current%prev
            else
              ! last element => search unsuccessful
              exit subk_loop
            end if
          end if

        end do node_loop
        
        ipst = ipnd+2
      end do subk_loop

      if (ntest.ge.100) then
        if (associated(node)) write(luout,*) 'success'
        if (.not.associated(node)) write(luout,*) 'no success'
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine next_node(cur_node,nxt_node,key)
*----------------------------------------------------------------------*
*     context is a string as e.g. "key.subkey.subsubkey"
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      type(keyword), pointer ::
     &     cur_node
      type(keyword), pointer ::
     &     nxt_node
      character, intent(in) ::
     &     key*(*)

      type(keyword), pointer ::
     &     current

      if (ntest.ge.100) then
        write(luout,*) '-----------'
        write(luout,*) ' next_node'
        write(luout,*) '-----------'
        write(luout,*) ' start = "',trim(cur_node%key),'"'
        write(luout,*) ' search for "',trim(key),'"'
      end if

      nxt_node => null()
      current => cur_node

      ! start search at head of next sublevel
      if (associated(current%down_h)) then
        current => current%down_h
      else
        ! rewind to first element in level
        do while(associated(current%prev))
          current => current%prev
        end do
      end if

      node_loop: do

        if (trim(current%key).eq.trim(key)) then
          ! matches key? return this node
          nxt_node => current
          exit node_loop
        end if
          
        ! find our way through tree (only same level and up)
        if (associated(current%next)) then
          current => current%next
        else if (associated(current%up)) then
          current => current%up
          ! rewind to first element in level
          do while(associated(current%prev))
            current => current%prev
          end do
        else
          ! root level reached => search unsuccessful
          exit node_loop
        end if

      end do node_loop
        
      if (ntest.ge.100) then
        if (associated(nxt_node)) write(luout,*) 'success'
        if (.not.associated(nxt_node)) write(luout,*) 'no success'
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine arg_node(node,cur_key,key)
*----------------------------------------------------------------------*
*     look whether keyword cur_key hosts an argument with key "key"
*     and return the corresponding node
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      type(keyword), pointer ::
     &     cur_key
      type(argument), pointer ::
     &     node
      character, intent(in) ::
     &     key*(*)

      type(keyword), pointer ::
     &     current
      type(argument), pointer ::
     &     curarg

      if (ntest.ge.100) then
        write(luout,*) '-----------'
        write(luout,*) ' arg_node'
        write(luout,*) '-----------'
        write(luout,*) ' start = "',trim(cur_key%key),'"'
        write(luout,*) ' search for "',trim(key),'"'
      end if

      node => null()
      current => cur_key

      ! start search at head of arguments
      if (associated(current%arg_h)) then
        curarg => current%arg_h
      else
        return
      end if

      node_loop: do

        if (trim(curarg%key).eq.trim(key)) then
          ! matches key? return this node
          node => curarg
          exit node_loop
        end if
          
        ! find our way through tree (only same level)
        if (associated(curarg%next)) then
          curarg => curarg%next
        else
          ! search unsuccessful
          exit node_loop
        end if

      end do node_loop
        
      if (ntest.ge.100) then
        if (associated(node)) write(luout,*) 'success'
        if (.not.associated(node)) write(luout,*) 'no success'
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine keyword_add(key,context,required,status)
*----------------------------------------------------------------------*
*     add a new keyword to level below context (default: keyword_root)
*     context is a string like "<key>.<key>.<key>"
*     where <key> is a valid key refering to an existing keyword node
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     key*(*)
      character, intent(in), optional ::
     &     context*(*)
      logical, intent(in), optional ::
     &     required
      integer, intent(in), optional ::
     &     status

      type(keyword), pointer ::
     &     current_up, current_prev, current

      current_up => keyword_root

      if (present(context)) then
        call find_node(keyword_root,current_up,context)
        if (.not.associated(current_up))
     &      call quit(1,'keyword_add','unknown context "'//context//'"')
      end if

      ! add new keyword at tail of list
      if (associated(current_up%down_t)) then
        ! tail of list exists
        current_prev => current_up%down_t
        allocate(current_prev%next)
        current => current_prev%next
        current_up%down_t => current
        current%up => current_up
        current%prev => current_prev
      else
        ! new list
        allocate(current_up%down_h)
        current => current_up%down_h
        current_up%down_t => current
        current%up => current_up
        nullify(current%prev)
      end if

      ! initialize node
      nullify(current%next)
      nullify(current%down_h)
      nullify(current%down_t)
      current%key = key
      if (present(status)) then
        current%status = status
      else
        current%status = -1
      end if
      nullify(current%arg_h)
      nullify(current%arg_t)
      if (present(required)) then
        current%required=required
      else
        current%required=.false.
      end if

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine argument_add(argkey,context,type,len,
     &     idef,xdef,ldef,cdef)
*----------------------------------------------------------------------*
*     add argument to keyword given by context
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     argkey*(*), context*(*)
      integer, intent(in), optional ::
     &     type,len
      integer, intent(in), optional ::
     &     idef(*)
      real(8), intent(in), optional ::
     &     xdef(*)
      logical, intent(in), optional ::
     &     ldef(*)
      character, intent(in), optional ::
     &     cdef(*)

      type(keyword), pointer ::
     &     curkey
      type(argument), pointer ::
     &     current, current_prev
      integer ::
     &     len_arg, type_arg

      call find_node(keyword_root,curkey,context)
      
      if (associated(curkey%arg_t)) then
        current_prev => curkey%arg_t
        allocate(current_prev%next)
        current => current_prev%next
        current%prev => current_prev
        curkey%arg_t => current
        nullify(current%next)
      else
        allocate(curkey%arg_h)
        current => curkey%arg_h
        curkey%arg_t => current
        nullify(current%prev)
        nullify(current%next)
      end if


      current%key = argkey
      current%status = -1
      type_arg = vtyp_log
      len_arg = 1
      if (present(type)) type_arg = type
      if (present(len))  len_arg = len      

      current%val%len = len_arg
      current%val%type = type_arg
      
      if (iand(type_arg,vtyp_log).gt.0) then
        if (present(ldef)) then
          allocate(current%val%lval(len_arg))
          current%val%lval(1:len_arg) = ldef(1:len_arg)
        end if
      end if
      if (iand(type_arg,vtyp_int).gt.0) then
        if (present(idef)) then
          allocate(current%val%ival(len_arg))
          current%val%ival(1:len_arg) = idef(1:len_arg)
        end if
      end if
      if (iand(type_arg,vtyp_rl8).gt.0) then
        if (present(xdef)) then
          allocate(current%val%xval(len_arg))
          current%val%xval(1:len_arg) = xdef(1:len_arg)
        end if
      end if
      if (iand(type_arg,vtyp_str).gt.0) then
        if (present(cdef)) then
          allocate(current%val%cval(len_arg))
          current%val%cval(1:len_arg) = cdef(1:len_arg)
        end if
      end if

      return
      end subroutine


*----------------------------------------------------------------------*
      subroutine keyword_add_history(context)
*----------------------------------------------------------------------*
*     add a new keyword to history
*     context is a string like "<key>.<key>.<key>"
*     where <key> is a valid key refering to an existing keyword node
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     context*(*)

      type(keyword), pointer ::
     &     current_up, current_prev, current
      integer ::
     &     idx_split

      idx_split = index(context,'.',back=.true.)
      if (idx_split.le.0) then
        call find_node(keyword_history,current_up,'',
     &       latest=.true.)
      else
        call find_node(keyword_history,current_up,
     &       context(1:idx_split-1),
     &       latest=.true.)
      end if
      if (.not.associated(current_up))
     &      call quit(1,'keyword_add_history',
     &     'unknown context "'//trim(context)//'"')

      ! add new keyword at tail of list
      if (associated(current_up%down_t)) then
        ! tail of list exists
        current_prev => current_up%down_t
        allocate(current_prev%next)
        current => current_prev%next
        current_up%down_t => current
        current%up => current_up
        current%prev => current_prev
      else
        ! new list
        allocate(current_up%down_h)
        current => current_up%down_h
        current_up%down_t => current
        current%up => current_up
        nullify(current%prev)
      end if

      ! initialize node
      nullify(current%next)
      nullify(current%down_h)
      nullify(current%down_t)
      current%key = context(idx_split+1:)
      current%status = -1
      nullify(current%arg_h)
      nullify(current%arg_t)

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine argument_add_history(context,arg_node,str)
*----------------------------------------------------------------------*
*     add argument to keyword history as given by context
*     the argument is given as string and parsed according to type
*     given in arg_node
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     context*(*)
      character, intent(in) ::
     &     str*(*)
      type(argument), pointer ::
     &     arg_node

      type(keyword), pointer ::
     &     curkey
      type(argument), pointer ::
     &     current, current_prev
      integer ::
     &     len_arg, type_arg

      call find_node(keyword_history,curkey,context,latest=.true.)
      
      if (associated(curkey%arg_t)) then
        current_prev => curkey%arg_t
        allocate(current_prev%next)
        current => current_prev%next
        current%prev => current_prev
        curkey%arg_t => current
        nullify(current%next)
      else
        allocate(curkey%arg_h)
        current => curkey%arg_h
        curkey%arg_t => current
        nullify(current%prev)
        nullify(current%next)
      end if

      ! get information on argument type from arg_node
      current%key      = arg_node%key
      current%val%len  = arg_node%val%len
      current%val%type = arg_node%val%type
      type_arg = current%val%type

      len_arg = check_array(str)
      current%val%len = len_arg

      ! allocate space for taking arguments
      if (iand(type_arg,vtyp_log).gt.0) then
        allocate(current%val%lval(len_arg))
        call read_log(current%val%lval,str)
      end if
      if (iand(type_arg,vtyp_int).gt.0) then
        allocate(current%val%ival(len_arg))
        call read_int(current%val%ival,str)
      end if
      if (iand(type_arg,vtyp_rl8).gt.0) then
        allocate(current%val%xval(len_arg))
        call read_real(current%val%xval,str)
      end if
      if (iand(type_arg,vtyp_str).gt.0) then
        allocate(current%val%cval(len_arg))
        call read_str(current%val%cval,str)
      end if

      return

      contains

      ! quick-and-dirty version:
      integer function check_array(str)
      implicit none

      character, intent(in) ::
     &     str*(*)
      integer ::
     &     ipos
      
      if (str(1:1).eq.'(') then
        ipos = 1
        check_array = 1
        do while(index(str(ipos:),',').gt.0)
          check_array = check_array+1
          ipos = index(str(ipos:),',')+ipos
        end do
      else if ! one number (very preliminary)
     &   (str(1:1).eq.'0'.or.str(1:1).eq.'1'.or.str(1:1).eq.'2'.or.
     &    str(1:1).eq.'3'.or.str(1:1).eq.'4'.or.str(1:1).eq.'5'.or.
     &    str(1:1).eq.'6'.or.str(1:1).eq.'7'.or.str(1:1).eq.'8'.or.
     &    str(1:1).eq.'9'.or.str(1:1).eq.'.'.or.str(1:1).eq.'-') then
        check_array = 1
      else ! a string
        check_array = index(str,',')
        if (check_array.eq.0) check_array = len_trim(str)
      end if

      end function

      subroutine read_log(larr,str)
      implicit none

      character, intent(in) ::
     &     str*(*)
      logical, intent(out) ::
     &     larr(*)
      integer ::
     &     ipst, ipnd, len, idx

      if (str(1:1).eq.'(') then
        ipst = 2
        len = len_trim(str)
        idx = 1
        do while(ipst.le.len)
          ipnd = index(str(ipst:),',')+ipst-2
          if (ipnd.lt.ipst) ipnd=len-1
          read(str(ipst:ipnd),*) larr(idx)
          ipst = ipnd+2
          idx = idx+1
        end do
      else
        read(str,*) larr(1)
      end if

      end subroutine

      subroutine read_int(iarr,str)
      implicit none

      character, intent(in) ::
     &     str*(*)
      integer, intent(out) ::
     &     iarr(*)
      integer ::
     &     ipst, ipnd, len, idx
      
      if (str(1:1).eq.'(') then
        ipst = 2
        len = len_trim(str)
        idx = 1
        do while(ipst.le.len)
          ipnd = index(str(ipst:),',')+ipst-2
          if (ipnd.lt.ipst) ipnd=len-1
          read(str(ipst:ipnd),*) iarr(idx)
          ipst = ipnd+2
          idx = idx+1
        end do
      else
        read(str,*) iarr(1)
      end if

      end subroutine

      subroutine read_real(xarr,str)
      implicit none

      character, intent(in) ::
     &     str*(*)
      real(8), intent(out) ::
     &     xarr(*)
      integer ::
     &     ipst, ipnd, len, idx
      
      if (str(1:1).eq.'(') then
        ipst = 2
        len = len_trim(str)
        idx = 1
        do while(ipst.le.len)
          ipnd = index(str(ipst:),',')+ipst-2
          if (ipnd.lt.ipst) ipnd=len-1
          read(str(ipst:ipnd),*) xarr(idx)
          ipst = ipnd+2
          idx = idx+1
        end do
      else
        read(str,*) xarr(1)
      end if

      end subroutine

      subroutine read_str(carr,str)
      implicit none

      character, intent(in) ::
     &     str*(*)
      character, intent(out) ::
     &     carr(*)
      integer ::
     &     ipst, ipnd, len, idx
      
      do idx = 1, len_trim(str)
        carr(idx) = str(idx:idx)
      end do

      end subroutine

      end subroutine

*----------------------------------------------------------------------*
      subroutine keyword_list(luout,tree_root,
     &     context,n_descent,show_args)
*----------------------------------------------------------------------*
*     print keyword list to unit luout
*     if context given, start here
*     if n_descent given, descent at most n_descent levels(default: all)
*     if show_args given, show possible argument keys to keyword
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     luout
      type(keyword), pointer ::
     &     tree_root
      character, intent(in), optional ::
     &     context*(*)
      integer, intent(in), optional ::
     &     n_descent
      logical, intent(in), optional ::
     &     show_args

      type(keyword), pointer ::
     &     current
      type(argument), pointer ::
     &     curarg
      logical ::
     &     shw_arg
      integer ::
     &     level, max_level
      character ::
     &     fmtstr*80
      
      level = 0
      current => tree_root
      max_level = 1000  ! means: show all levels
      shw_arg = .false.
      if (present(context))
     &     call find_node(tree_root,current,context)
      if (present(n_descent))
     &     max_level = n_descent
      if (present(show_args))
     &     shw_arg = show_args

      key_loop: do

        ! show keyword
        if (current%status.gt.0) then
          write(fmtstr,'("(""A"",",i3,"x,a)")') 2*level+1
        else
          write(fmtstr,'("(""I"",",i3,"x,a)")') 2*level+1
        end if
        write(luout,fmtstr) trim(current%key)
        ! show arguments to keyword (if applicable)
        if (shw_arg.and.associated(current%arg_h)) then
          curarg => current%arg_h
          write(fmtstr,'("(x,",i3,"x,a,i2,x,i2)")') 2*level+4

          arg_loop: do
            
            write(luout,fmtstr) trim(curarg%key),
     &         curarg%val%type,curarg%val%len
            if (allocated(curarg%val%ival))
     &           write(luout,*) ' > ',curarg%val%ival
            if (allocated(curarg%val%xval))
     &           write(luout,*) ' > ',curarg%val%xval
            if (allocated(curarg%val%cval))
     &           write(luout,*) ' > ',curarg%val%cval
            if (allocated(curarg%val%lval))
     &           write(luout,*) ' > ',curarg%val%lval

            if (.not.associated(curarg%next)) exit arg_loop
            curarg => curarg%next

          end do arg_loop

        end if

        ! find a way through the tree:
        if (associated(current%down_h) .and. level.lt.max_level) then
          ! go down
          current => current%down_h
          level = level+1
        else if (associated(current%next)) then
          ! else stay within level
          current => current%next
        else if (level.gt.0) then
          ! else find an upper level, where a next
          ! node exists:
          up_loop: do
            if (associated(current%up)) then
              current => current%up
              level = level-1
              if (associated(current%next)) then
                current => current%next
                exit up_loop
              end if
            else
              exit key_loop
            end if
          end do up_loop
        else
          exit key_loop
        end if

      end do key_loop

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine keyword_get_context(context,cur_node)
*----------------------------------------------------------------------*
      implicit none

      character, intent(out) ::
     &     context*(*)
      type(keyword), pointer ::
     &     cur_node

      integer, parameter ::
     &     maxscr = 1024

      type(keyword), pointer ::
     &     current
      character ::
     &     cscr*(maxscr)

      current => cur_node
      
      context = ""
      if (associated(current%up)) context = trim(current%key)

      do while (associated(current%up)) 
        current => current%up
        cscr = context
        if (len_trim(cscr)+len_trim(current%key)+1.gt.maxscr)
     &       call quit(1,'keyword_get_context','cscr too small')
        if (associated(current%up))
     &       context = trim(current%key)//'.'//trim(cscr)
      end do

      return
      end subroutine

*----------------------------------------------------------------------*
      subroutine find_active_node(tree,node,context,icount)
*----------------------------------------------------------------------*

      implicit none

      type(keyword), pointer ::
     &     tree
      type(keyword), pointer ::
     &     node
      character ::
     &     context*(*)
      integer ::
     &     icount

      character ::
     &     curcontext*1024     
      type(keyword), pointer ::
     &     current
      integer ::
     &     jcount

      if (.not.associated(tree%down_h))
     &     call quit(1,'find_active_node','invalid keyword tree')

      current => tree%down_h

      node => null()

      jcount = 0
      key_loop: do 
        if (current%status.gt.0) then
          call keyword_get_context(curcontext,current)
          if (trim(context).eq.trim(curcontext)) jcount = jcount+1 
          if (icount.eq.jcount) then
            node => current
            exit key_loop
         end if
        end if

        if (current%status.gt.0.and.associated(current%down_h)) then
          ! go down
          current => current%down_h
        else if (associated(current%next)) then
          ! else stay within level
          current => current%next
        else
          ! else find an upper level, where a next
          ! node exists:
          up_loop: do
            if (associated(current%up)) then
              current => current%up
              if (associated(current%next)) then
                current => current%next
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
      subroutine keyword_parse(luin)
*----------------------------------------------------------------------*
*     parse the keywords on unit luin
*     the unit should be a formatted, sequential file, positioned
*     at the place where the parser should start
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'ifc_baserout.h'

      integer, intent(in) ::
     &     luin

      integer, parameter ::
     &     maxlen  = 256,
     &     n_delim = 9
      logical, parameter ::
     &     new_line_is_delim = .true.
      character, parameter ::
     &     delimiter(n_delim) = 
     &     (/' ', ';', ',', '(', ')', '\\', '"', '!', '=' /)
      integer, parameter ::
     &     ispace =   1,
     &     isemicolon = 2,
     &     icomma =   3,
     &     iparen_o = 4,
     &     iparen_c = 5,
     &     ibacksl  = 6,
     &     iquote   = 7,
     &     icomment = 8,
     &     iequal   = 9
      integer, parameter ::
     &     n_allowed_start = 1,
     &     allowed_start(n_allowed_start) = (/icomment/),
     &     n_allowed_after_key = 7,
     &     allowed_after_key(n_allowed_after_key) = 
     &     (/ispace,iequal,isemicolon,
     &     iparen_o,ibacksl,iquote,icomment/),
     &     n_allowed_after_arg = 5,
     &     allowed_after_arg(n_allowed_after_arg) = 
     &     (/isemicolon,iequal,icomma,ibacksl,icomment/)
      character ::
     &     line*(maxlen)

      type(keyword), pointer ::
     &     curkey,nxtkey
      type(argument), pointer ::
     &     curarg
      character*(maxlen) ::
     &     context
      integer ::
     &     allowed_delim(n_delim), n_allowed_delim
      integer ::
     &     ipst, ipnd, itest, lenline, ierr

      context = " "
      curkey => keyword_root

      allowed_delim(1:n_allowed_start)=allowed_start(1:n_allowed_start)
      n_allowed_delim = n_allowed_start

      ierr = 0
      file_loop: do
        read(luin,'(a)',end=100,err=200) line

c        ipst = first_nonblank(line)
        call clean_line(line,delimiter,n_delim)

        lenline = len_trim(line)
        ! empty line?
        if (lenline.le.0) cycle file_loop 
        ipst = 1

        itest = next_delim(line(ipst:ipst))
        
        if (itest.le.0) then
          ierr = ierr+1
          call error_delim(line,ipst)
          ipst = ipst+1
        end if
        
        ! comment?
        if (line(ipst:ipst).eq.delimiter(icomment)) cycle file_loop 

        allowed_delim(1:n_allowed_after_key) =
     &       allowed_after_key(1:n_allowed_after_key)
        n_allowed_delim = n_allowed_after_key

        line_loop: do while(ipst.le.lenline)
          itest = next_delim(line(ipst:))
          ipnd = abs(itest)+ipst-1
          if (itest.le.0) then
            ierr = ierr+1
            call error_delim(line,ipnd)
          end if
          ipnd = ipnd-1
          ! is it an argument key?
          call arg_node(curarg,curkey,line(ipst:ipnd))

          if (associated(curarg)) then
            allowed_delim(1) = iequal
            n_allowed_delim = 1
            itest = next_delim(line(ipnd+1:ipnd+1))
            if (itest.le.0) then
              ierr = ierr+1
              call error_misseq(line,ipnd+1)
            end if

            ipst = ipnd+2
            if (ipst.le.lenline) then
              if (line(ipst:ipst).eq.delimiter(iparen_o)) then
                ipnd = index(line(ipst:),delimiter(iparen_c))+ipst-1
                if (ipnd.lt.ipst) then
                  ierr = ierr + 1
                  call error_misspc(line,ipnd)
                  exit line_loop
                end if
                
              else

                allowed_delim(1:n_allowed_after_arg) =
     &               allowed_after_arg(1:n_allowed_after_arg)
                n_allowed_delim = n_allowed_after_arg
                itest = next_delim(line(ipst:))
                ipnd = abs(itest)+ipst-1
                if (itest.le.0) then
                  ierr = ierr+1
                  call error_delim(line,ipnd)
                end if
                ipnd = ipnd-1
              end if

              call keyword_get_context(context,curkey)
              call argument_add_history(context,curarg,line(ipst:ipnd))

            else
              ierr = ierr+1
              call error_eol(line,ipst)
              allowed_delim(1:n_allowed_after_key) =
     &             allowed_after_key(1:n_allowed_after_key)
              n_allowed_delim = n_allowed_after_key
              exit line_loop
            end if

          else
            ! find keyword in tree
            call next_node(curkey,nxtkey,line(ipst:ipnd))

            if (.not.associated(nxtkey)) then
              ierr = ierr+1
              call error_keywd(line,ipst,curkey)
            else
              curkey => nxtkey
              call keyword_get_context(context,curkey)
              ! add node to keyword history
              call keyword_add_history(context)

            end if          
          end if

          ipst = ipnd+2
        end do line_loop
      end do file_loop
      
 100  continue

      if (iprlvl.ge.10)
     &    call keyword_list(luout,keyword_history,show_args=.true.)

      if (ierr.gt.0)
     &     call quit(0,'parse_input','input errors detected, see above')

      return

 200  call quit(0,'parse_input','I/O error on reading input file')
      
      return
      
*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*
*     local functions
*     variables of main function can be accessed
*----------------------------------------------------------------------*
      integer function next_delim(str)
*----------------------------------------------------------------------*
*     find next delimiter
*     return value:
*       position of next delimiter
*       len+1 if line end is delimiter
*       negative value, if delimiter is not allowed in context
*----------------------------------------------------------------------*

      implicit none

      character, intent(in) ::
     &     str*(*)

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
      subroutine error_delim(str,ipos)
*----------------------------------------------------------------------*
      implicit none
      
      character, intent(in) ::
     &     str*(*)
      integer, intent(in) ::
     &     ipos
      character ::
     &     fmtstr*80
      write(luout,'(x,a)') trim(line)
      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(luout,fmtstr) 
      write(fmtstr,'("(x,a,",i3,"(a1,x))")') n_allowed_delim
      write(luout,fmtstr) 'INPUT ERROR: unexpected delimiter, '//
     &     'expected one of ',
     &     delimiter(allowed_delim(1:n_allowed_delim))

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
      character ::
     &     fmtstr*80

      write(luout,'(x,a)') trim(line)
      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(luout,fmtstr) 
      write(luout,'(x,a)') 'INPUT ERROR: unexpected end-of-line'
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
      character ::
     &     fmtstr*80

      write(luout,'(x,a)') trim(line)
      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(luout,fmtstr) 
      write(luout,'(x,a)') 'INPUT ERROR: missing )'
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
      character ::
     &     fmtstr*80

      write(luout,'(x,a)') trim(line)
      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(luout,fmtstr) 
      write(luout,'(x,a)') 'INPUT ERROR: missing ='
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
      type(keyword), intent(in) ::
     &     curkey
      character ::
     &     fmtstr*80

      write(luout,'(x,a)') trim(line)
      write(fmtstr,'("(x,",i3,"x,""^"")")') abs(ipos)-1
      write(luout,fmtstr) 
      write(luout,'(x,a)') 'INPUT ERROR: unexpected keyword '

      return
      end subroutine

      end subroutine keyword_parse

*----------------------------------------------------------------------*

      end module parse_input

