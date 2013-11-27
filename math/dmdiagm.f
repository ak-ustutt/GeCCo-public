*----------------------------------------------------------------------*
      subroutine dmdiagm(ndim,mat,diag,right,inverted)
*----------------------------------------------------------------------*
*     performs multiplication of a matrix with a diagonal matrix
*     right    : multiply diagonal matrix from right (T) or left (F)
*     inverted : use diag^-1 (T) or just diag (F)
*
*     matthias, dec 2010
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00
      real(8), parameter ::
     &     thresh = 1d-7 ! because in current use we have the
                         ! square root of diagonal elements on diag
      integer, intent(in) ::
     &     ndim
      real(8), intent(inout) ::
     &     mat(ndim,ndim)
      real(8), intent(in) ::
     &     diag(ndim)
      logical, intent(in) ::
     &     right, inverted

      integer ::
     &     idx

      if (ndim.eq.0) return

      if (ntest.ge.100) then
        write(lulog,'(x,a)') '---------------'
        write(lulog,'(x,a)') 'dmdiagm at work'
        write(lulog,'(x,a)') '---------------'
        write(lulog,*) 'input matrix:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
        write(lulog,*) 'normalization vector:'
        call wrtmat2(diag,1,ndim,1,ndim)
      end if

      if (inverted) then
        if (right) then
          do idx = 1, ndim
            if (diag(idx).lt.thresh) then
c              mat(1:ndim,idx) = 0d0
              cycle
            end if
            mat(1:ndim,idx) = mat(1:ndim,idx)/diag(idx)
          end do
        else
          do idx = 1, ndim
            if (diag(idx).lt.thresh) then
c              mat(idx,1:ndim) = 0d0
              cycle
            end if
            mat(idx,1:ndim) = mat(idx,1:ndim)/diag(idx)
          end do
        end if
      else
        if (right) then
          do idx = 1, ndim
            if (diag(idx).lt.thresh) cycle ! debug
            mat(1:ndim,idx) = mat(1:ndim,idx)*diag(idx)
          end do
        else
          do idx = 1, ndim
            if (diag(idx).lt.thresh) cycle ! debug
            mat(idx,1:ndim) = mat(idx,1:ndim)*diag(idx)
          end do
        end if
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'output matrix:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      return
      end
