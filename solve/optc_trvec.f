
      subroutine optc_trvec(imet,xvec,ndim,xvec0)

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 0

      integer, intent(in) ::
     &     imet, ndim
      real(8), intent(in) ::
     &     xvec0(*)
      real(8), intent(out) ::
     &     xvec(*)

      integer ::
     &     ii

      if (ntest.ge.100) then
        write(luout,*) 'modus: ',imet
        write(luout,*) 'input vector: '
        call wrtmat2(xvec0,1,ndim+1,1,ndim+1)
      end if

      if (imet.eq.0) then

        do ii = 1, ndim
          xvec(ii) = xvec0(ii) - xvec0(ii+1)
        end do

      else

        do ii = 1, ndim
          xvec(ii) = xvec0(ii) - xvec0(ndim+1)
        end do

      end if

      if (ntest.ge.100) then
        write(luout,*) 'output vector: '
        call wrtmat2(xvec,1,ndim,1,ndim)
      end if

      return
      end
