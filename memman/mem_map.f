      subroutine mem_map(check)
      use memman
      implicit none
      include 'stdunit.h'

      logical, intent(in) ::
     &     check

      call memman_map(luout,check)

      return
      end
