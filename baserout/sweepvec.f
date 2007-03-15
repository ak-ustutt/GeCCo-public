
*----------------------------------------------------------------------*

      subroutine sweepvec(vec,ndim)

* purpose: replace numerical zeroes by real zeroes 
*          (convenient for debugging)

      implicit none

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     vec(ndim)
      
      integer ::
     &     i
      real(8) ::
     &     thr

      thr = 100d0*epsilon(1d0)
      do i = 1, ndim
        if (abs(vec(i)).lt.thr) vec(i) = 0d0
      end do

      return
      end
