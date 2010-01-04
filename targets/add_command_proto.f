*----------------------------------------------------------------------*
      subroutine add_command_proto(name,tgt_info)
*----------------------------------------------------------------------*
*     allocate a new slot for an (prototype) action structure 
*     in tgt_info
*     the entries are (currently) used for setting default values
*     the list is extended and the array is updated
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     name

      type(action_list), pointer ::
     &     list_pnt

      integer ::
     &     itgts, idx, iarg

      integer, external ::
     &     idx_action

      idx = idx_action(name,tgt_info)
      if (idx.gt.0)
     &     call quit(1,'add_action',
     &     'action already exists: '//trim(name) )

      list_pnt => tgt_info%act_list
      ! advance to end of target list:
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%act)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        list_pnt%next => null()
      end if
      allocate (list_pnt%act)

      ! add user-supplied label
      if (len_trim(name).gt.len_command_name)
     &   call quit(1,'add_command_proto','name too long: "'
     &     //trim(name)//'"')
      list_pnt%act%command(1:len_command_name) = ' '
      list_pnt%act%command = trim(name)
      list_pnt%act%n_arguments = 0
      list_pnt%act%arg => null()

      ! increment counter
      tgt_info%nactions = tgt_info%nactions+1

      ! update target array
      call update_tgt_info(tgt_info)

      return
      end
