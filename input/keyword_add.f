*----------------------------------------------------------------------*
      subroutine keyword_add(key,context,required,status)
*----------------------------------------------------------------------*
*     add a new keyword to level below context (default: keyword_root)
*     context is a string like "<key>.<key>.<key>"
*     where <key> is a valid key refering to an existing keyword node
*----------------------------------------------------------------------*

      use parse_input
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
