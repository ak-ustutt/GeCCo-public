*----------------------------------------------------------------------*
      subroutine inv_svd(vx,amat,vb,ndim,inc_b)
*----------------------------------------------------------------------*
*     solve A*x + b = 0 by singular value decomposition
*     A(ndim,ndim), vx(ndim), vb(inc_b,ndim)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim, inc_b
      real(8), intent(in) ::
     &     amat(ndim,ndim)
      real(8), intent(inout) ::
     &     vx(ndim)
      real(8), intent(in) ::
     &     vb(inc_b,ndim)

      real(8), parameter ::
     &     thrsh = 1d-12

      real(8), pointer ::
     &     umat(:,:), vtmat(:,:), wrk(:), singval(:)

      integer ::
     &     idx, kdx, info, lwrk

      lwrk=max(1024,ndim*ndim)
      allocate(umat(ndim,ndim),vtmat(ndim,ndim),
     &     wrk(lwrk),singval(ndim))

      call dgesvd('A','A',ndim,ndim,
     &     amat,ndim,singval,
     &     umat,ndim,vtmat,ndim,
     &     wrk,lwrk,info)

      if (info.ne.0) then
        write(luout,*) 'WARNING: SVD seems to be in trouble'
      end if

      do idx = 1, ndim
        wrk(idx) = 0d0
        if (abs(singval(idx)).lt.thrsh) cycle
        do kdx = 1, ndim
          wrk(idx) = wrk(idx)+
     &         umat(kdx,idx)*vb(1,kdx)
        end do
        wrk(idx) = wrk(idx)/singval(idx)
      end do
        
      do idx = 1, ndim
        vx(idx) = 0d0
        do kdx = 1, ndim
          vx(idx) = vx(idx) +
     &         vtmat(kdx,idx)*wrk(kdx)
        end do
      end do
      
      return
      end

