*----------------------------------------------------------------------*
      subroutine act_list2arr(act_list,act_arr,nacts)
*----------------------------------------------------------------------*
*     set up an array of pointers that point to the entries of the
*     list on act_list
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'def_target.h'
      include 'def_action_list.h'
      include 'def_action_array.h'

      type(action_list), intent(in), target ::
     &     act_list
      integer, intent(in) ::
     &     nacts
      type(action_array), intent(out) ::
     &     act_arr(nacts)

      type(action_list), pointer ::
     &     current

      integer ::
     &     iact

      current => act_list

      do iact = 1, nacts
        if (.not.associated(current%act))
     &       call quit(1,'act_list2arr','unallocated action on list')
        act_arr(iact)%act => current%act
        if (iact.lt.nacts.and..not.(associated(current%next)))
     &       call quit(1,'act_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
