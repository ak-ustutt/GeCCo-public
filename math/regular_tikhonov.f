*----------------------------------------------------------------------*
      subroutine regular_tikhonov(ncol,nlin,mat,sv,omega2)
*----------------------------------------------------------------------*
*     multiplies the columns of the input matrix with a regularization
*     factor based on the corresponding singular value,
*     in order to achieve Tikhonov regularization.
*     We assume that the input matrix already contains sv^(-1/2),
*     such that the needed factor is sv/sqrt(sv^2+omega2)
*
*     matthias, jul 2012
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ncol, nlin
      real(8), intent(inout) ::
     &     mat(nlin,ncol)
      real(8), intent(in) ::
     &     sv(ncol), omega2

      real(8) ::
     &     fac
      integer ::
     &     icol
      real(8), external ::
     &     ddot

      if (nlin.eq.0) return

c dbg
c          write(*,*) 'before Tikhonov:'
c          call wrtmat2(mat,nlin,ncol,nlin,ncol)
c          write(*,*) 'singular values:'
c          call wrtmat2(sv,1,ncol,1,ncol)
c dbgend

      do icol = 1, ncol
        fac = sv(icol)/sqrt(sv(icol)**2+omega2)
c        ! modified regularization: s^-1/2 --> 1/(s^2+omega2)^0.25
c        fac = sqrt(sv(icol))/((sv(icol)**2+omega2)**0.25d0)
c        ! mod: just unitary matrix
c        fac = sqrt(sv(icol))
        mat(1:nlin,icol) = fac * mat(1:nlin,icol)
      end do
c dbg
c          write(*,*) 'after Tikhonov:'
c          call wrtmat2(mat,nlin,ncol,nlin,ncol)
c          write(*,*) 'norms:'
c          do icol = 1, ncol
c            fac = ddot(nlin,mat(1,icol),1,mat(1,icol),1)
c            if (fac.gt.1d-12) print *,icol,fac
c          end do
c dbgend

      return
      end
