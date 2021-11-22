*----------------------------------------------------------------------*
      subroutine file_close_keep(fhand)
*----------------------------------------------------------------------*
*     close and keep file with handle fhand (unit number is reset)
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     fhand

      if (fhand%buffered) return

      if (fhand%unit.le.0)
     &     call quit(1,'file_close_keep',
     &     'illegal unit number on handle')

      call relunit(fhand%unit,'keep')
      fhand%unit = -1

      return
      end
