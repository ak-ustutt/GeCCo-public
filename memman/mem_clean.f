      subroutine mem_clean()
      use memman
      implicit none
      include 'stdunit.h'

      call memman_stat(lulog)
      if (luout.ne.lulog) call memman_stat(luout)
      call memman_clean()

      return
      end
