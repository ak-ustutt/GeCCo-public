
*----------------------------------------------------------------------*

      integer function ifndmax(ivec,idxoff,lvec,inc)

      implicit none

      integer, intent(in) ::
     &     ivec(*), lvec, inc, idxoff

      integer ::
     &     i, imx, idx

      imx = -huge(imx)
      idx = idxoff
      do i = 1, lvec
        imx = max(imx,ivec(idx))
        idx = idx + inc
      end do
      
      ifndmax = imx

      return
      end
