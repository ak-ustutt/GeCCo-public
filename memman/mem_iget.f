*----------------------------------------------------------------------*
      subroutine mem_iget(ffbuf,ibuf,idxst,idxnd)
*----------------------------------------------------------------------*
*     get integer_array(idxst:idxnd) on ibuf(1:)
*----------------------------------------------------------------------*
      use memman
      implicit none
c      include 'def_filinf.h'

      type(filinf) ::
     &     ffbuf
      integer, intent(in) ::
     &     idxst, idxnd
      integer, intent(out) ::
     &     ibuf(*)

      call memman_vbuffer_get_int(ffbuf%buf_id,idxst,idxnd,ibuf)

      return
      end
