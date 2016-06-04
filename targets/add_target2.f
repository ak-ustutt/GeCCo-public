*----------------------------------------------------------------------*
      subroutine add_target2(name,required,tgt_info)
*----------------------------------------------------------------------*
*     allocate a new slot for a target structure in tgt_info
*     the list is extended and the array is updated
*     new version: type is no longer required
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      logical, intent(in) ::
     &     required
      character, intent(in) ::
     &     name*(*)

      type(target_list), pointer ::
     &     list_pnt

      integer ::
     &     itgts, idx

      integer, external ::
     &     idx_target

      idx = idx_target(name,tgt_info)
      if (idx.gt.0)
     &     call quit(1,'add_target2',
     &     'target already exists: '//trim(name) )

      list_pnt => tgt_info%list
      ! advance to end of target list:
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%tgt)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        list_pnt%next => null()
      end if
      allocate (list_pnt%tgt)

      ! add user-supplied label
      if (len_trim(name).gt.len_target_name)
     &   call quit(1,'add_target2','name too long: "'
     &     //trim(name)//'"')
      list_pnt%tgt%name(1:len_target_name) = ' '
      list_pnt%tgt%type = ttype_gen ! now always generic
      list_pnt%tgt%name = name

      ! increment counter
      tgt_info%ntargets = tgt_info%ntargets+1
      ! set my_idx and last_mod
      list_pnt%tgt%my_idx = tgt_info%ntargets
      list_pnt%tgt%last_mod = -1

      ! set other defaults
      list_pnt%tgt%required = required  ! primary target?
      list_pnt%tgt%n_joined_with = 0
      list_pnt%tgt%n_depends_on = 0
      list_pnt%tgt%joined_with => null()
      list_pnt%tgt%depends_on => null()
      list_pnt%tgt%idx_joined_with => null()
      list_pnt%tgt%idx_depends_on => null()

      list_pnt%tgt%n_rules = 0
      list_pnt%tgt%rules => null()

      ! update target array
      call update_tgt_info(tgt_info)

      return
      end
