*----------------------------------------------------------------------*
      subroutine file_delete(fhand)
*----------------------------------------------------------------------*
*     delete file with handle fhand
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     fhand

      ! not open? 
      if (fhand%unit.le.0) call file_open(fhand)

      call relunit(fhand%unit,'delete')
      fhand%unit = -1

      return
      end
