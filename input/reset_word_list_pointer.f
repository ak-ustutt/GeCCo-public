      subroutine reset_word_list_pointer(wlist)
      !
      ! set position "current" to first entry
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist

      type(word_list_entry), pointer ::
     &     wl_pnt
      integer ::
     &     len_word

      if (.not.associated(wlist%head)) 
     &     call quit(1,'reset_word_list_entry','wlist is undefined')

      wlist%current => wlist%head

      end
