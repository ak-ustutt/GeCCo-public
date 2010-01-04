*----------------------------------------------------------------------*
      subroutine invsqrt_mat(ndim,mat)
*----------------------------------------------------------------------*
*     calculates mat^(-0.5)
*     mat must be quadratic and symmetric (not checked so far)
*
*     matthias, dec 2009
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00
      real(8), parameter ::
     &     min_sv = 1d-15 ! singular value threshold for calc. of pseudo-inv.
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
          eigen_vec(1:ndim,idx) = eigen_vec(1:ndim,idx)
     &                         * (eigen_val(idx)**(-0.25d0))
        else
          eigen_vec(1:ndim,idx) = 0d0
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'U*s^(-0.25):'
        call wrtmat2(eigen_vec,ndim,ndim,ndim,ndim)
      end if

      call dgemm('n','t',ndim,ndim,ndim,
     &           1d0,eigen_vec,ndim,
     &               eigen_vec,ndim,
     &           0d0,mat,ndim)

      deallocate(eigen_vec,eigen_val)

      return
      end
