*----------------------------------------------------------------------*
      subroutine mem_print_vbuffer(luout,ffbuf)
*----------------------------------------------------------------------*
*     print info on virtual memory buffer (for debugging)
*----------------------------------------------------------------------*
      use memman
      implicit none
c      include 'def_filinf.h'

      type(filinf) ::
     &     ffbuf
      integer, intent(in) ::
     &     luout

      call print_vbuffer_int(luout,ffbuf%buf_id)

      return
      end
