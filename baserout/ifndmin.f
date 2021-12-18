
*----------------------------------------------------------------------*

      integer function ifndmin(ivec,idxoff,lvec,inc)

      implicit none

      integer, intent(in) ::
     &     ivec(*), lvec, inc, idxoff

      integer ::
     &     i, imn, idx

      imn = huge(imn)
      idx = idxoff
      do i = 1, lvec
        imn = min(imn,ivec(idx))
        idx = idx + inc
      end do
      
      ifndmin = imn

      return
      end
