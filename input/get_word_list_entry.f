      subroutine get_word_list_entry(word,sep,wlist)
      !
      ! get entry from position "current"
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist
      character(len=*), intent(inout) ::
     &     word
      character(len=1), intent(inout) ::
     &     sep

      type(word_list_entry), pointer ::
     &     wl_pnt
      integer ::
     &     len_word

      len_word = len(word)

      word(1:len_word) = ' '
      sep = ' '
      if (.not.associated(wlist%head)) return

      if (len_trim(wlist%current%word).gt.len_word)
     &     call quit(1,'get_word_list_entry','output string too small')

      word = wlist%current%word
      sep  = wlist%current%sep

      end
