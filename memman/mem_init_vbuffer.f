*----------------------------------------------------------------------*
      subroutine mem_init_vbuffer(ffbuf,name,max_len,max_slots)
*----------------------------------------------------------------------*
*     assign a memory buffer to file ffbuf
*----------------------------------------------------------------------*
      use memman
      implicit none
c      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     ffbuf

      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     max_len,max_slots

      call memman_init_vbuffer(ffbuf,name,max_len,max_slots)

      return
      end
