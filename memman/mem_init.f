      subroutine mem_init(mem_free_init)
      use memman
      implicit none

      integer, intent(in) ::
     &     mem_free_init
      
      call memman_init(mem_free_init)

      return
      end
