*----------------------------------------------------------------------*
      subroutine get_argument_value(
     &     context,argkey,keycount,argcount,
     &     ival,iarr,lval,larr,xval,xarr,string)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input2
      use FoX_dom,only:Node,Nodelist,DOMException,
     &     getFirstChild,hasChildNodes,getChildNodes,getLength,item,
     &     getNodeName,getAttribute,getTextContent
      use FoX_common,only:rts
      implicit none
      integer,parameter::
     &     ntest=1000
      character(len=18),parameter ::
     &     i_am="get_argument_value"

      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      logical, intent(out), optional ::
     &     lval, larr(*)
      integer, intent(out), optional ::
     &     ival, iarr(*)
      real(8), intent(out), optional ::
     &     xval, xarr(*)
      character, intent(out), optional ::
     &     string*(*)
      type(Node), pointer ::
     &     curkey,nxtkey
      type(Node), pointer ::
     &     curarg,input_root,key_root
      type(NodeList),pointer::
     &     child_list
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target, type, dim, idx
      logical ::
     &     try_default, succ ,dummy
      integer ::
     &     ii, num,len,ex

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (present(keycount)) write (lulog,*) "keycount:",keycount
         if (present(argcount)) write (lulog,*) "argcount:",argcount
         write (lulog,*) "looking for ",argkey," in context:",context
      end if

      input_root=>inp_fetch_root()
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')
      
      key_root=>reg_fetch_root()


      icount_target = 1
      if (present(keycount)) icount_target = keycount
      
      try_default = .false.
      call find_active_node(input_root,curkey,
     &     context,icount_target)
      iargcount = 0
      iargcount_target = 1
      if (present(argcount)) iargcount_target = argcount

      succ = .false.
      ! either get 
      repetition_loop: do 
! try default instead

         if (associated(curkey))then
            if (.not.hasChildNodes(curkey)
     &           .and.iargcount_target.eq.1) then
               try_default = .true.
            end if
         else if(.not. associated(curkey))then
            try_default = .true.
         end if

         if (try_default)
     &        call find_node(key_root,curkey,context)


         if (associated(curkey).and.hasChildNodes(curkey)) then
            child_list=> getChildNodes(curkey)
            arg_loop: do ii=0,getLength(child_list)-1 ! MY LISTS START AT 1!!!!
               curarg=>item(child_list,ii)
               if (getNodeName(curarg) .ne. arg_tag) cycle arg_loop
               if (getAttribute(curarg,atr_name).eq.trim(argkey)) 
     &              iargcount = iargcount+1

               if (getAttribute(curarg,atr_name).eq.trim(argkey).and.
     &              iargcount.eq.iargcount_target) then
                  
                  call get_argument_dimension_core(curarg,len,type,
     &                 dummy)

                  call rts(getAttribute(curarg,atr_len),dim)
                  
                  if (ntest.ge.100) then 
                     write(lulog,'(" dim:",i3)')dim 
                     write(lulog,'(" type:",i3)')type
                     write(lulog,'("unconverted input:",a,":")')
     &                    getAttribute(curarg,atr_val)
                  end if 

                  if(.not.hasAttribute(curarg,atr_val))then
                     num=0
                     exit arg_loop
                  else
                     num=1
                  end if
                  select case(type)
               case (vtyp_log)
                  if (.not.(present(lval).or.present(larr)))
     &                 call quit(1,i_am,
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no l-value array present')
                  if (present(lval)) 
     &                 call rts(getAttribute(curarg,atr_val),lval,
     &                 iostat=ex)
                  if (present(larr))  
     &                 call rts(getAttribute(curarg,atr_val),
     &                 larr(1:len),iostat=ex,num=num) 
                  if (ex.le. 0) succ = .true.
               case (vtyp_int)
                  if (.not.(present(ival).or.present(iarr)))
     &                 call quit(1,i_am,
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no i-value array present')
                  if (present(ival)) then
                     call rts(getAttribute(curarg,atr_val),
     &                 ival,iostat=ex,num=num)
                  end if 
                  if (present(iarr)) 
     &                 call rts(trim(getAttribute(curarg,atr_val)),
     &                 iarr(1:len),iostat=ex,num=num)
                  if (ex.le. 0) succ = .true.
            
               case (vtyp_rl8)
                  if (.not.(present(xval).or.present(xarr)))
     &                 call quit(1,i_am,
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no r-value array present')
                  if (present(xval)) 
     &                 call rts(getAttribute(curarg,atr_val),xval
     &                 ,iostat=ex)
                  if (present(xarr)) 
     &                 call rts(getAttribute(curarg,atr_val),
     &                 xarr(1:len),iostat=ex,num=num) 
                  if (ex.le. 0) succ = .true.
               case (vtyp_str)
                  if (.not.(present(string)))
     &                 call quit(1,'get_argument_value',
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no r-value array present')
                  string = trim(getAttribute(curarg,atr_val))
                  succ = .true.
               end select
               if (ntest.ge.100) then 
                  write(lulog,'(" transferred ",i3," elements")') num
                  write(lulog,*) succ 
               end if
               end if
               if (succ) exit repetition_loop
            end do arg_loop
         end if 
         if (try_default) exit repetition_loop
!     else try default (if applicable)
         try_default = .true.
         call find_node(key_root,curkey,context)
         
         if (.not.associated(curkey).or.
     &        .not.hasChildNodes(curkey)) exit repetition_loop
      end do repetition_loop 
            
      if (.not.succ)
     &     call quit(1,i_am,
     &     'Could not provide any value for '//trim(context)//
     &     '.'//trim(argkey))

      return
      end
