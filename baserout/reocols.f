      subroutine reocols(matout,matin,ireo,nrow,ncol)

      implicit none

      integer, intent(in) ::
     &     ncol, nrow, ireo(ncol)
      real(8), intent(in) ::
     &     matin(ncol,nrow)
      real(8), intent(out) ::
     &     matout(ncol,nrow)

      integer ::
     &     idx

      do idx = 1, nrow
        matout(1:ncol,idx) = matin(1:ncol,ireo(idx))
      end do

      return
      end
