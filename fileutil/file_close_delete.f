
*----------------------------------------------------------------------*
      subroutine file_close_delete(fhand)
*----------------------------------------------------------------------*
*     close and delete file with handle fhand (unit number is reset)
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     fhand

      if (fhand%unit.le.0)
     &     call quit(1,'file_close_delete',
     &     'illegal unit number on handle')

      call relunit(fhand%unit,'delete')
      fhand%unit = -1

      return
      end
