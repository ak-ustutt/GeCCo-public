*----------------------------------------------------------------------*
      subroutine mgs(xmat,smat,ndim,xscr)
*----------------------------------------------------------------------*
*     modified Gram-Schmidt orthonormalization
*     adapted from jeppe olsen's mgs3() routine (LUCIA)
*
*     input:  smat(ndim,ndim) - overlap matrix of vectors <v_i|v_j>
*     output: xmat(ndim,ndim) - orthogonalization matrix:
*                               |w_i> = |v_j> x(j,i)
*     scratch: xscr(ndim)
*
*     zero column vectors may be returned in case of linear dependency
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     smat(ndim,ndim)
      real(8), intent(out) ::
     &     xmat(ndim,ndim)
      real(8), intent(inout) ::
     &     xscr(ndim)

      integer ::
     &     ii, jj
      real(8) ::
     &     xnorm, factor, xsx

      real(8), external ::
     &     ddot

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'mgs')
        write(luout,*) 'overlap-matrix:'
        call wrtmat2(smat,ndim,ndim,ndim,ndim)
      end if

      ! set unit matrix 
      xmat(1:ndim,1:ndim) = 0d0
      do ii = 1, ndim
        xmat(ii,ii) = 1d0
      end do

      ! loop over vectors
      do ii = 1, ndim
        call dgemv('n',ndim,ndim,1d0,smat,ndim,xmat(1,ii),1,
     &                 0d0,xscr,1)
        xnorm = ddot(ndim,xmat(1,ii),1,xscr,1)
        if (xnorm.lt.epsilon(1d-6)) then
          xmat(1:ndim,ii) = 0d0
          cycle
        else
          factor = 1d0/sqrt(xnorm)
        end if
        xmat(1:ndim,ii) = factor*xmat(1:ndim,ii)
        xscr(1:ndim) = factor * xscr(1:ndim)

        ! subtract xmat(1:ndim,ii) from all remaining vectors
        do jj = ii+1, ndim
          xsx = ddot(ndim,xscr,1,xmat(1,jj),1)
          xmat(1:ndim,jj) = xmat(1:ndim,jj) - xsx*xmat(1:ndim,ii)
        end do
c dbg
c        write(luout,*) 'updated X-matrix:'
c        call wrtmat2(xmat,ndim,ndim,ndim,ndim)
c dbg

      end do

      if (ntest.ge.100) then
        write(luout,*) 'trafo-matrix:'
        call wrtmat2(xmat,ndim,ndim,ndim,ndim)
      end if

      return
      end
