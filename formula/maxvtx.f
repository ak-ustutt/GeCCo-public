
*----------------------------------------------------------------------*
      integer function maxvtx(iconn,narc)
*----------------------------------------------------------------------*
*     return maximum vertex number currently in use 
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     narc, iconn(3,narc)

      integer ::
     &     idx

      maxvtx = -1
      do idx = 1, narc
        maxvtx = max(maxvtx,iconn(2,idx),iconn(3,idx))
      end do

      return
      end
