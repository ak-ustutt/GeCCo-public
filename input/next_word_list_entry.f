      logical function next_word_list_entry(word,wlist)
      !
      ! get entry from position "current"
      ! return false, if no further entry is available
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist
      character(len=*), intent(inout) ::
     &     word

      type(word_list_entry), pointer ::
     &     wl_pnt
      integer ::
     &     len_word

      len_word = len(word)

      word(1:len_word) = ' '
      next_word_list_entry = .false.
      if (.not.associated(wlist%head)) return

      if (len_trim(wlist%current%word).gt.len_word)
     &     call quit(1,'next_word_list_entry','output string too small')

      word = wlist%current%word

      next_word_list_entry = associated(wlist%current%next)

      if (next_word_list_entry) wlist%current => wlist%current%next

      end
