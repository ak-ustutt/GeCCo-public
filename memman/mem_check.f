      subroutine mem_check(label)
      use memman
      implicit none
      include 'stdunit.h'

      character, intent(in) ::
     &     label*(*)

      call memman_check(lulog,label)

      return
      end
