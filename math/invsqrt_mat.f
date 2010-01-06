*----------------------------------------------------------------------*
      subroutine invsqrt_mat(ndim,mat)
*----------------------------------------------------------------------*
*     calculates U*mat^(-0.5) using MAT = U*mat*U^+
*     mat must be quadratic and symmetric (not checked so far)
*
*     matthias, dec 2009
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00
      real(8), parameter ::
     &     min_sv = 1d-10 ! singular value threshold for calc. of pseudo-inv.
      integer, intent(in) ::
     &     ndim
      real(8), intent(inout), target ::
     &     mat(ndim,ndim)

      integer ::
     &     nrot, idx
      real(8), pointer ::
     &     eigen_val(:), eigen_vec(:,:)

      allocate(eigen_vec(ndim,ndim),eigen_val(ndim))

      ! Calculate eigenvalues and eigenvectors
      ! NOTE: May/Should be replaced by svd-routine in the long run!
      call jacobi(mat,ndim,ndim,eigen_val,eigen_vec,nrot)

      if (ntest.ge.100) then
        write(luout,*) 'eigen_val: ',eigen_val
        write(luout,*) 'eigen_vec:'
        call wrtmat2(eigen_vec,ndim,ndim,ndim,ndim)
      end if

      ! square root of (pseudo) inverse:
      ! A^(-1/2) = U D^(-1/2) U^+ = U D^(-1/4) [U D^(-1/4)]^+
      do idx = 1, ndim
        if (eigen_val(idx).gt.min_sv) then
c          eigen_vec(1:ndim,idx) = eigen_vec(1:ndim,idx)
c     &                         * (eigen_val(idx)**(-0.25d0))
          mat(1:ndim,idx) = eigen_vec(1:ndim,idx)
     &                         * (eigen_val(idx)**(-0.5d0))
        else
c          eigen_vec(1:ndim,idx) = 0d0
          mat(1:ndim,idx) = 0d0
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'U*s^(-0.5):'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

c      call dgemm('n','t',ndim,ndim,ndim,
c     &           1d0,eigen_vec,ndim,
c     &               eigen_vec,ndim,
c     &           0d0,mat,ndim)

      deallocate(eigen_vec,eigen_val)

      return
      end
