*----------------------------------------------------------------------*
      subroutine argument_add(argkey,context,type,len,
     &     idef,xdef,ldef,cdef)
*----------------------------------------------------------------------*
*     add argument to keyword given by context
*----------------------------------------------------------------------*
      
      use parse_input
      implicit none
      include 'par_vtypes.h'

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
