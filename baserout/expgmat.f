*----------------------------------------------------------------------*
*     follows: exp and log of matrices
*----------------------------------------------------------------------*
      subroutine expgmat(ndim,expx,xmat,xscr1,xscr2,thrsh)
*----------------------------------------------------------------------*
*     calculate exp(X), returned on expx, of (ndim,ndim)-matrix X,
*     input on xmat, by Taylor expansion (threshold thrsh)
*     xscr is a scratch matrix of the same dimensions as xmat, expx
*
*     any quadratic matrix may be supplied
*
*     andreas, aug 2004
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter :: ntest = 00, maxn = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     thrsh
      real(8), intent(inout) ::
     &     expx(ndim,ndim), xmat(ndim,ndim),
     &     xscr1(ndim,ndim),  xscr2(ndim,ndim)

      logical ::
     &     conv
      integer ::
     &     n, ndim2, ii
      real(8) ::
     &     xnrm, fac

      real(8), external ::
     &     inprod

      expx(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)
      xscr2(1:ndim,1:ndim) = xmat(1:ndim,1:ndim)

      do ii = 1, ndim
        expx(ii,ii) = expx(ii,ii) + 1d0
      end do

      ndim2 = ndim*ndim
      n = 1
      conv = .false.

      do while (.not.conv)
        n = n+1
        if (n.gt.maxn) exit

        fac = 1d0/dble(n)

        ! Xscr = 1/N Xscr * X
c        call matml7(xscr1,xscr2,xmat,
c     &              ndim,ndim,
c     &              ndim,ndim,
c     &              ndim,ndim,0d0,fac,0)
        call dgemm('n','n',ndim,ndim,ndim,
     &             fac,xscr2,ndim,
     &                 xmat,ndim,
     &             0d0,xscr1,ndim)


        xnrm = sqrt(inprod(xscr1,xscr1,ndim2))
        if (xnrm.lt.thrsh) conv = .true.

        if (ntest.ge.10)
     &       write(luout,*) ' N = ',n,'  |1/N! X^N| = ',xnrm

        expx(1:ndim,1:ndim) = expx(1:ndim,1:ndim) + xscr1(1:ndim,1:ndim)

        xscr2(1:ndim,1:ndim) = xscr1(1:ndim,1:ndim)

      end do

      if (.not.conv) then
        write(luout,*) ' Taylor expansion of exp(X) did not converge!'
        call quit(1,'expgmat','no convergence')
      end if

      return
      end
