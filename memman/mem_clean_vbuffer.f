*----------------------------------------------------------------------*
      subroutine mem_clean_vbuffer(name)
*----------------------------------------------------------------------*
*     reset buffer, deallocate all slots
*----------------------------------------------------------------------*
      use memman
      implicit none

      character*(*) ::
     &     name

      call memman_clean_vbuffer(name)

      return
      end
