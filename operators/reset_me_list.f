*----------------------------------------------------------------------*
      subroutine reset_me_list(name,op_info)
*----------------------------------------------------------------------*
*     reset me_list with label "name"
*     i.e. set modification mark to undefined and delete file
*     associated with list
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

      mel => list_pnt%mel

      if (associated(mel%fhand%last_mod))
c     &     mel%fhand%last_mod = 0
     &     mel%fhand%last_mod = -1

      call file_delete(mel%fhand)

      return
      end
