*----------------------------------------------------------------------*
      subroutine get_argument_value(
     &     context,argkey,keycount,argcount,
     &     ival,iarr,lval,larr,xval,xarr,str)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input2,only : find_active_node,find_node,
     &     arg_tag,atr_name,atr_len,atr_kind,atr_val,
     &     input_root,input_doc,key_root
      use FoX_dom,only:Node,Nodelist,DOMException,
     &     getFirstChild,hasChildNodes,getChildNodes,getLength,item,
     &     getNodeName,getAttribute
      use FoX_common,only:rts
      implicit none
      include 'par_vtypes.h'
      integer,parameter::
     &     ntest=00
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
     &     str*(*)
      integer::
     &     ex
      type(Node), pointer ::
     &     curkey
      type(Node), pointer ::
     &     curarg
      type(NodeList),pointer::
     &     child_list
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target, type, dim, idx
      logical ::
     &     try_default, succ
      integer ::
     &     ii, num

      input_root=getFirstChild(input_doc)
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')


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
         if ((.not.associated(curkey).or.
     &     .not.hasChildNodes(curkey)) !careful, attributes are child nodes
     &     .and.iargcount_target.eq.1) then
            try_default = .true.
            call find_node(key_root,curkey,context)
         end if 

         if (associated(curkey).and.hasChildNodes(curkey)) then
            child_list=> getChildNodes(curkey)
            arg_loop: do ii=1,getLength(child_list) 
               curarg=>item(child_list,ii)

               if (getNodeName(curarg) .ne. arg_tag) cycle arg_loop
               
               if (getAttribute(curarg,atr_name).eq.trim(argkey)) 
     &              iargcount = iargcount+1
               
               if (getAttribute(curarg,atr_name).eq.trim(argkey).and.
     &              iargcount.eq.iargcount_target) then
                  call rts(getAttribute(curarg,atr_len),dim)
                  call rts(getAttribute(curarg,atr_kind),type)
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
     &                 call rts(getAttribute(curarg,atr_val),larr(1:dim)
     &                 ,iostat=ex) 
                  if (ex.gt. 0) succ = .true.
               case (vtyp_int)
                  if (.not.(present(ival).or.present(iarr)))
     &                 call quit(1,i_am,
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no i-value array present')
                  if (present(lval)) 
     &                 call rts(getAttribute(curarg,atr_val),ival
     &                 ,iostat=ex)
                  if (present(larr)) 
     &                 call rts(getAttribute(curarg,atr_val),iarr(1:dim)
     &                 ,iostat=ex) 
                  if (ex.gt. 0) succ = .true.
            
               case (vtyp_rl8)
                  if (.not.(present(xval).or.present(xarr)))
     &                 call quit(1,i_am,
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no r-value array present')
                  if (present(lval)) 
     &                 call rts(getAttribute(curarg,atr_val),xval
     &                 ,iostat=ex)
                  if (present(larr)) 
     &                 call rts(getAttribute(curarg,atr_val),xarr(1:dim)
     &                 ,iostat=ex) 
                  if (ex.gt. 0) succ = .true.
               case (vtyp_str)
                  if (.not.(present(str)))
     &                 call quit(1,'get_argument_value',
     &                 trim(context)//'->'//trim(argkey)//
     &                 'no r-value array present')
                  do idx = 1, dim
                     str = trim(getAttribute(curarg,atr_val))
                  end do
                  succ = .true.
                  end select
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
