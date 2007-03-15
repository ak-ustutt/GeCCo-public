
*----------------------------------------------------------------------*
      subroutine nrmvec(ndim,eigvec,eigvi)
*----------------------------------------------------------------------*
*     normalize the eigenvectors in array eigvec(ndim,ndim)
*     imaginary pairs are handled as described in rg(), eispack.f
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     eigvi(ndim)
      real(8), intent(inout) ::
     &     eigvec(ndim,ndim)

      integer ::
     &     ivec
      real(8) ::
     &     xnrm

      real(8), external ::
     &     inprod

      ivec = 0
      do while (ivec.lt.ndim)
        ivec = ivec+1
        
        xnrm = inprod(eigvec(1,ivec),eigvec(1,ivec),ndim)

        if (eigvi(ivec).ne.0d0) then
          if (ivec+1.gt.ndim) then
            write(6,*) 'inconsistency in eigenvalue structure'
            stop 'nrmvec'
          end if

          xnrm = xnrm + inprod(eigvec(1,ivec+1),eigvec(1,ivec+1),ndim)
        end if

        xnrm = sqrt(xnrm)

        call scalve(eigvec(1,ivec),1d0/xnrm,ndim)

        if (eigvi(ivec).ne.0d0) then
          call scalve(eigvec(1,ivec+1),1d0/xnrm,ndim)
          ivec = ivec+1
        end if

      end do

      return

      end
