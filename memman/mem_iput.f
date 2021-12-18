*----------------------------------------------------------------------*
      subroutine mem_iput(ffbuf,ibuf,idxst,idxnd)
*----------------------------------------------------------------------*
*     put ibuf(1:) to integer_array(idxst:idxnd)
*----------------------------------------------------------------------*
      use memman
      implicit none
c      include 'def_filinf.h'

      type(filinf) ::
     &     ffbuf
      integer, intent(in) ::
     &     idxst, idxnd
      integer, intent(in) ::
     &     ibuf(*)

      call memman_vbuffer_put_int(ffbuf%buf_id,idxst,idxnd,ibuf)

      return
      end
