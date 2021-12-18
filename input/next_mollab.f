*----------------------------------------------------------------------*
      subroutine next_mollab(label,lu,ierr)
*----------------------------------------------------------------------*
*
*  adapted from mollab (DALTON)
*
*  Purpose:
*     Forward to next MOLECULE label on file LU
*
*----------------------------------------------------------------------*
      implicit none

      character*8, intent(out) ::
     &     label
      integer, intent(out) ::
     &     ierr
      integer, intent(in) ::
     &     lu
      character*8 ::
     &     mark(4)
      character*8, parameter ::
     &     magic = '********'

      ierr = 0 
      do
        read (lu,end=3,err=6) mark
        if (mark(1).eq.magic) then
          label = mark(4)
          exit
        end if
      end do

      return

 3    ierr = -1
      return

 6    ierr = -2
      return

      end
