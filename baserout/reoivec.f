      subroutine reoivec(vec,ireo,ndim)

      implicit none

      integer, intent(in) ::
     &     ndim, ireo(ndim)
      integer, intent(inout) ::
     &     vec(ndim)

      integer ::
     &     idx, tmp(ndim)

      tmp = vec

      do idx = 1, ndim
        vec(idx) = tmp(ireo(idx))
      end do

      return
      end
