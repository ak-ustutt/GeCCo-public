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

      use parse_input2, only : Node,NodeList,
     &     input_root, arg_tag,atr_name,atr_len,atr_kind,
     &     getFirstChild, hasChildNodes, getChildNodes, getLength,
     &     getNodeName, getAttribute, item, 
     &     find_node,find_active_node
      use FoX_common, only : rts
      implicit none

      integer,parameter::
     &     ntest=00
      character(len=22),parameter ::
     &     i_am="get_argument_dimension"

      integer, intent(out) ::
     &     dim
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      integer, intent(out), optional ::
     &     type

      type(Node), pointer ::
     &     curkey
      type(Node), pointer ::
     &     curarg

      type(NodeList), pointer ::
     &     curargs

      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, icount_target, iargcount_target, ii 

      print *, "ntered get_argument_dimension" 

      if (.not.hasChildNodes(input_root))
     &     call quit(1,i_am,'invalid keyword history')

      icount_target = 1
      if (present(keycount)) icount_target = keycount

      print *, "active node finding"
      call find_active_node(input_root,curkey,
     &     context,icount_target)
      print *, "active node found"
      dim = -1
      if (present(type)) type = -1
      iargcount = 0
      iargcount_target = 1
      if (present(argcount)) iargcount_target = argcount

      if (associated(curkey).and.hasChildNodes(curkey)) then
        curargs => getChildNodes(curkey)
        arg_loop: do ii=1,getLength(curargs)
        curarg=>item(curargs, ii)
          if (getNodeName(curarg) .eq. arg_tag) then 
             if (trim(getAttribute(curarg,atr_name))
     &            .eq. trim(argkey) ) iargcount = iargcount+1

             if (iargcount.eq.iargcount_target) then
                call rts(getAttribute(curarg, atr_len), dim)
                if (present(type)) 
     &               call rts( getAttribute(curarg,atr_kind),
     &               type)
                exit arg_loop
             end if
          end if 
        end do arg_loop
!> TODO convert length and kind into normalized attributes
      end if

      return
      end
