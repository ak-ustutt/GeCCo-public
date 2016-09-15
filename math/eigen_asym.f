*----------------------------------------------------------------------*
      subroutine eigen_asym(ndim,xmat,eigr,eigi,vecs,xscr,ierr)
*----------------------------------------------------------------------*
*     Solve the eigenvalue problem for an asymmetric matrix
*     Wrapper for eispack call
*     covers some ordering and renormalization issues as needed
*     in GeCCo
*
*     xmat -> matrix to be diagonalised
*     eigr -> real part of the resulting eigenvalues
*     eigi -> imaginary part of the resulting eigenvalues
*     vecs -> resulting eigenvectors
*     xscr -> scratch variable (?)
*     ierr -> error
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     xmat(ndim,ndim)
      real(8), intent(inout) ::
     &     eigr(ndim), eigi(ndim), vecs(ndim,ndim),
     &     xscr(ndim,ndim)
      integer, intent(inout) ::
     &     ierr

      integer ::
     &     irg, idx,
     &     ivec(ndim)
      real(8) ::
     &     rnorm,
     &     xvec(ndim)

      real(8), external ::
     &     ddot, dnrm2

      irg = 1
      call rg(ndim,ndim,xmat,eigr,eigi,irg,vecs,ivec,xvec,ierr)
      if (ierr.ne.0) return

      ! normalize and re-orthogonalize degenerate pairs
      idx = 1
      do while(idx.le.ndim)
        if (eigi(idx).eq.0d0) then
          rnorm = dnrm2(ndim,vecs(1,idx),1)
          vecs(1:ndim,idx) = (1d0/rnorm)*vecs(1:ndim,idx)
          idx = idx+1
        else
          rnorm = dnrm2(2*ndim,vecs(1,idx),1)
          vecs(1:ndim,idx) = (1d0/rnorm)*vecs(1:ndim,idx)
          vecs(1:ndim,idx+1) = (1d0/rnorm)*vecs(1:ndim,idx+1)
          ! re-ortho. part still missing (cf. lucia_ccrsp, l3110ff)
          idx = idx+2
        end if
      end do

      ! init reordering array
      do idx = 1, ndim
        ivec(idx) = idx
      end do
      ! sort solutions according to real part (ascending)
      ! eigr is resorted on the way
      call idxsortx(eigr,ivec,ndim,+1)

      ! resort eigi as well
      call reovec(xvec,eigi,ivec,ndim)
      eigi(1:ndim) = xvec(1:ndim)

      ! resort vectors
      call reocols(xscr,vecs,ivec,ndim,ndim)
      vecs = xscr

      return
      end
