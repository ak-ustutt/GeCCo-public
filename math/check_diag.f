*----------------------------------------------------------------------*
      subroutine check_diag(ndim,nvec,mat,eigenvec,eigenval)
*----------------------------------------------------------------------*
*     Check if eigenvec are really eigenvectors of mat associated
*     to the eigenvalues eigenval, giving a report
*
*     ndim - dimension of the matrix
*     nvec - number of vectors
*     mat - the ndim x ndim matrix to be diagonalised
*     eigenvec - ndim x nvec matrix with eigenvectors
*     eigenval - nvec vector with eigenvalues
*
*     yuri, oct 2014
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim, nvec
      real(8), intent(in) ::
     &     mat(ndim,ndim)
      real(8), intent(out) ::
     &     eigenvec(ndim,nvec), eigenval(nvec)

      integer ::
     &     i,j,k
      real(8) ::
     &     l_vec(ndim,nvec), m_vec(ndim,nvec), rmsd(nvec)

      do i=1,ndim
       do j=1,nvec
        m_vec(i,j) = 0d0
        do k=1,ndim
         m_vec(i,j) = m_vec(i,j) + mat(i,k)*eigenvec(k,j)
        end do
       end do
      end do

      do i=1,ndim
       rmsd(nvec) = 0d0
       do j=1,nvec
        l_vec(i,j) = eigenval(j)*eigenvec(i,j)
        rmsd(j) = rmsd(j) + (l_vec(i,j)-m_vec(i,j))**2
       end do
       rmsd(j) = sqrt(rmsd(j)/nvec)
      end do
      
      write(lulog,*) "====================="
      write(lulog,*) "   check_diag"
      write(lulog,*) "Initial matrix, A:"
      do i=1,ndim
       write(lulog,*) mat(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "Proposed eigenvalues, v:"
      write(lulog,*) eigenval(:)
      write(lulog,*) "---"
      write(lulog,*) "Proposed eigenvectors. l:"
      do i=1,ndim
       write(lulog,*) eigenvec(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "A x v:"
      do i=1,ndim
       write(lulog,*) m_vec(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "l x v:"
      do i=1,ndim
       write(lulog,*) l_vec(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "Root mean square deviations:"
      write(lulog,*) rmsd(:)
      write(lulog,*) "---"      
      write(lulog,*) "====================="

      end
