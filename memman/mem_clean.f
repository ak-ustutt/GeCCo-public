      subroutine mem_clean()
      use memman
      implicit none
      include 'stdunit.h'

      call memman_stat(lulog)
      call memman_clean()

      return
      end
