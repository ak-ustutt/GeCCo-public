      subroutine print_word_list(luout,wlist)

      implicit none

      include 'def_word_list.h'

      integer, intent(in) ::
     &     luout
      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt

      ! not active?
      if (.not.associated(wlist%head)) then
        write(luout,*) 'word list is empty'
        return
      end if

      wl_pnt => wlist%head

      ! loop over all entries
      do 

        write(luout,*) '"',trim(wl_pnt%word),'", sep=',wl_pnt%sep

        if (.not.associated(wl_pnt%next)) exit
        wl_pnt => wl_pnt%next
      end do

      end
