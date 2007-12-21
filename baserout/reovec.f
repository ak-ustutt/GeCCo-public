      subroutine reovec(vecout,vecin,ireo,ndim)

      implicit none

      integer, intent(in) ::
     &     ndim, ireo(ndim)
      real(8), intent(in) ::
     &     vecin(ndim)
      real(8), intent(out) ::
     &     vecout(ndim)

      integer ::
     &     idx

      do idx = 1, ndim
        vecout(idx) = vecin(ireo(idx))
      end do

      return
      end
