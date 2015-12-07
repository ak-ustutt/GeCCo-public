*----------------------------------------------------------------------*
      subroutine check_diag(ndim,nvec,mat,eigenvec,eigenval,s_mat)
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
      real(8), intent(in) ::
     &     eigenvec(ndim,nvec), eigenval(nvec)
      real(8), intent(in), optional ::
     &     S_mat(ndim,ndim)

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
        if (present(S_mat)) then
         l_vec(i,j) = 0.0d0
         do k=1,ndim
          l_vec(i,j) = l_vec(i,j) + S_mat(i,k)*eigenvec(k,j)
         end do
        else
         l_vec(i,j) = eigenvec(i,j)
        end if
        l_vec(i,j) = eigenval(j)*l_vec(i,j)
        rmsd(j) = rmsd(j) + (l_vec(i,j)-m_vec(i,j))**2
       end do
       rmsd(j) = sqrt(rmsd(j)/nvec)
      end do
      
      write(lulog,*) "====================="
      write(lulog,*) "   check_diag"
      write(lulog,*) "Initial matrix, H:"
      do i=1,ndim
       write(lulog,*) mat(i,:)
      end do
      write(lulog,*) "---"
      if (present(S_mat)) then
       write(lulog,*) "Overlap matrix, S:"
       do i=1,ndim
        write(lulog,*) S_mat(i,:)
       end do
       write(lulog,*) "---"
      end if
      write(lulog,*) "Proposed eigenvalues, E:"
      write(lulog,*) eigenval(:)
      write(lulog,*) "---"
      write(lulog,*) "Proposed eigenvectors, v:"
      do i=1,ndim
       write(lulog,*) eigenvec(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "H x v:"
      do i=1,ndim
       write(lulog,*) m_vec(i,:)
      end do
      write(lulog,*) "---"
      if (present(S_mat)) then
       write(lulog,*) "E x S x v:"
      else
       write(lulog,*) "E x v:"
      end if
      do i=1,ndim
       write(lulog,*) l_vec(i,:)
      end do
      write(lulog,*) "---"
      write(lulog,*) "Root mean square deviations:"
      write(lulog,*) rmsd(:)
      write(lulog,*) "---"      
      write(lulog,*) "====================="

      end
