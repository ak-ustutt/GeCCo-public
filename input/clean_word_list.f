      subroutine clean_word_list(wlist)

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt, wl_next

      ! not active?
      if (.not.associated(wlist%head)) return

      wl_pnt => wlist%head

      ! deallocate all entries
      do 
        wl_next => wl_pnt%next
        deallocate(wl_pnt)
        if (.not.associated(wl_next)) exit
        wl_pnt => wl_next
      end do

      ! reset all pointers
      wlist%head => null()
      wlist%tail => null()
      wlist%current => null()

      end
