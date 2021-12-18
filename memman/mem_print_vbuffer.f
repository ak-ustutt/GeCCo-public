*----------------------------------------------------------------------*
      subroutine mem_print_vbuffer(lulog,ffbuf)
*----------------------------------------------------------------------*
*     print info on virtual memory buffer (for debugging)
*----------------------------------------------------------------------*
      use memman
      implicit none
c      include 'def_filinf.h'

      type(filinf) ::
     &     ffbuf
      integer, intent(in) ::
     &     lulog

      call print_vbuffer_int(lulog,ffbuf%buf_id,.false.)

      return
      end
