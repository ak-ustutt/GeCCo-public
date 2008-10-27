*----------------------------------------------------------------------*
      subroutine add_action(name,type,narguments,
     &                      arg_label,arg_type,tgt_info)
*----------------------------------------------------------------------*
*     allocate a new slot for an (prototype) action structure 
*     in tgt_info
*     the list is extended and the array is updated
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      integer, intent(in) ::
     &     type, narguments, arg_type(narguments)
      character(len=*), intent(in) ::
     &     name, arg_label(*)

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
     &   call quit(1,'add_action','name too long: "'
     &     //trim(name)//'"')
      list_pnt%act%command(1:len_command_name) = ' '
      list_pnt%act%type = type
      list_pnt%act%command = trim(name)
      list_pnt%act%n_update = 1 ! preliminary fix
      list_pnt%act%n_arguments = narguments
      if (narguments.eq.0) then
        list_pnt%act%arg => null()
      else
        allocate(list_pnt%act%arg(narguments))
        do iarg = 1, narguments
          list_pnt%act%arg(iarg)%type      = arg_type(iarg)
          list_pnt%act%arg(iarg)%arg_label = trim(arg_label(iarg))
          list_pnt%act%arg(iarg)%val_label => null()
          list_pnt%act%arg(iarg)%val_log => null()
          list_pnt%act%arg(iarg)%val_int => null()
          list_pnt%act%arg(iarg)%val_occ => null()
          list_pnt%act%arg(iarg)%val_restr => null()
          list_pnt%act%arg(iarg)%val_rl8  => null()
          list_pnt%act%arg(iarg)%val_str  => null()
        end do

      end if

      ! increment counter
      tgt_info%nactions = tgt_info%nactions+1

      ! update target array
      call update_tgt_info(tgt_info)

      return
      end
