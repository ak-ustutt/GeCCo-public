*----------------------------------------------------------------------*
      subroutine check_inv(ndim,mat,inv)
*----------------------------------------------------------------------*
*     Check if inv is the inverse matrix of mat, giving a report
*
*     ndim - dimension of the matrix
*
*     yuri, set 2015
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     mat(ndim,ndim), inv(ndim,ndim)

      integer ::
     &     i,j,k
      real(8) ::
     &     K_delta(ndim,ndim)

      do i=1,ndim
       do j=1,ndim
        K_delta(i,j) = 0d0
        do k=1,ndim
         K_delta(i,j) = K_delta(i,j) + mat(i,k)*inv(k,j)
        end do
       end do
      end do

      write(lulog,*) "====================="
      write(lulog,*) "   check_inv"
      write(lulog,*) "Initial matrix, A:"
      do i=1,ndim
       write(lulog,*) mat(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "Proposed inverse matrix, A^-1:"
      do i=1,ndim
       write(lulog,*) inv(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "Their product:"
      do i=1,ndim
       write(lulog,*) K_delta(i,:)
      end do
      write(lulog,*) "====================="

      end
