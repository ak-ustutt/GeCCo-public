*----------------------------------------------------------------------*
!>     set up an array of pointers that point to the entries of the
!>     list on mel_list
!>     new version with a real pointer array
!>     @param[in] mel_list linked list of ME-lists
!>     @param[in] nmels number of ME-Lists
!>     @param[out] mel_arr array to be populated
*----------------------------------------------------------------------*
      subroutine mel_list2arr(mel_list,mel_arr,nmels)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_me_list_list.h'
      include 'def_me_list_array.h'

      type(me_list_list), intent(in), target ::
     &     mel_list
      integer, intent(in) ::
     &     nmels
      type(me_list_array), intent(out) ::
     &     mel_arr(nmels)

      type(me_list_list), pointer ::
     &     current

      integer ::
     &     imel

      current => mel_list

      do imel = 1, nmels
        if (.not.associated(current%mel))
     &       call quit(1,'mel_list2arr','unallocated ME-list on list')
        mel_arr(imel)%mel => current%mel
        if (imel.lt.nmels.and..not.(associated(current%next)))
     &       call quit(1,'mel_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
