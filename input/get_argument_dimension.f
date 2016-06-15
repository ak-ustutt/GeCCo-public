


*----------------------------------------------------------------------*
      subroutine get_argument_dimension(dim,
     &     context,argkey,keycount,argcount,type)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use parse_input2,only : find_active_node,find_node,
     &     arg_tag,atr_name,atr_len,atr_kind,atr_val,
     &     input_root,input_doc,key_root,hasAttribute,
     &     get_argument_dimension_core
      use FoX_dom,only:Node,Nodelist,DOMException,
     &     getFirstChild,hasChildNodes,getChildNodes,getLength,item,
     &     getNodeName,getAttribute,getTextContent
      use FoX_common,only:rts
      implicit none
      include 'par_vtypes.h'
      include 'stdunit.h'

      integer,parameter::
     &     ntest=00
      character(len=22),parameter ::
     &     i_am="get_argument_dimension"
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      integer, intent(out),optional::
     &     type




      type(Node), pointer ::
     &     curkey,nxtkey
      type(Node), pointer ::
     &     curarg
      type(NodeList),pointer::
     &     child_list
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target, itype, dim_tot, 
     &     idx
      logical ::
     &     try_default, succ
      integer ::
     &     ii, dim, ex

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (present(keycount)) write (lulog,*) "keycount:",keycount
         if (present(argcount)) write (lulog,*) "argcount:",argcount
         write (lulog,*) "looking for ",argkey," in context:",context
      end if


      input_root=getFirstChild(input_doc)
      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')
      
      if (ntest.ge.100) then
         write(lulog,*) "assoc. status of input_root",
     &        associated(input_root)
      end if

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

                  call get_argument_dimension_core(curarg,dim,
     &                 itype,succ)
               if (ntest.ge.100) then 
                  write(lulog,'(" argument has ",i3," elements")') dim 
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
      if (present(type))type=itype
      return
      end
