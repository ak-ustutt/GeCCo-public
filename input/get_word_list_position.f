      subroutine get_word_list_position(line,col,wlist)
      !
      ! get line and column (in original file) from "current" entry
      !

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist
      integer, intent(out) ::
     &     line, col

      type(word_list_entry), pointer ::
     &     wl_pnt

      line = 0
      col  = 0
      if (.not.associated(wlist%head)) return

      line = wlist%current%line
      col  = wlist%current%col

      end
