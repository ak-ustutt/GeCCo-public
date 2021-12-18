*----------------------------------------------------------------------*
      subroutine del_me_list(name,op_info)
*----------------------------------------------------------------------*
*     delete info for me_list with label "name"
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      character(*), intent(in) ::
     &     name

      type(me_list_list), pointer ::
     &     list_pnt
      type(me_list), pointer ::
     &     mel

      list_pnt => op_info%mel_list
      ! advance to respective entry on me_list list:
      do while (trim(list_pnt%mel%label).ne.trim(name)
     &     .and.associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      if (trim(list_pnt%mel%label).ne.trim(name)) then
        call quit(1,'del_me_list','unknown label: "'//trim(name)//'"')
      end if

      if (associated(list_pnt%prev)) list_pnt%prev%next => list_pnt%next
      if (associated(list_pnt%next)) list_pnt%next%prev => list_pnt%prev
      mel => list_pnt%mel
      call dealloc_me_list(mel)
      deallocate(list_pnt%mel)
      deallocate(list_pnt)

      ! decrement counter
      op_info%nmels = op_info%nmels-1

      ! update ME-list array
      call update_mel_arr(op_info)
      ! make sure that also the operator->ME assignments get updated
      call update_op_arr(op_info)

      return
      end
