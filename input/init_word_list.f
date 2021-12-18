      subroutine init_word_list(wlist)

      implicit none

      include 'def_word_list.h'

      type(word_list), intent(inout) ::
     &     wlist

      wlist%head => null()
      wlist%tail => null()
      wlist%current => null()

      end
