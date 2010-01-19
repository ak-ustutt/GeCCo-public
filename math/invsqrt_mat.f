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
     &     ntest = 100
      real(8), parameter ::
     &     min_sv = 1d-10, ! singular value threshold for calc. of pseudo-inv.
     &     warn_sv = 1d-5  ! give a warning for small singular values
      integer, intent(in) ::
     &     ndim
      real(8), intent(inout), target ::
     &     mat(ndim,ndim)
      real(8) ::
     &     singval(ndim),wrk(max(1024,ndim**2)),dum1,dum2

    

      integer ::
     &     nrot, idx, lwrk, info

      lwrk=max(1024,ndim**2)
      info = 0

      ! calculate U and singular values:
      call dgesvd('O','N',ndim,ndim,
     &     mat,ndim,singval,
     &     dum1,1,dum2,1,
     &     wrk,lwrk,info)

      if (info.ne.0) then
        write(luout,*) 'WARNING in invsqrt_mat: SVD in trouble'
      end if

      if (ntest.ge.100) then
        write(luout,*) 'singular values: ',singval
        write(luout,*) 'eigenvectors:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      ! square root of (pseudo) inverse:
      ! A^(-1/2) = U D^(-1/2) U^+ = U D^(-1/4) [U D^(-1/4)]^+
      do idx = 1, ndim
        if (singval(idx).gt.min_sv) then
          if (singval(idx).lt.warn_sv)
     &         call warn('invsqrt_mat','small singular values!')
c          eigen_vec(1:ndim,idx) = eigen_vec(1:ndim,idx)
c     &                         * (singval(idx)**(-0.25d0))
          mat(1:ndim,idx) = mat(1:ndim,idx)
     &                         * (singval(idx)**(-0.5d0))
        else
c          eigen_vec(1:ndim,idx) = 0d0
          mat(1:ndim,idx) = 0d0
        end if
c dbg   can be (mis)used to set a unit operator
c        mat(1:ndim,idx) = 0d0
c        mat(idx,idx)  =1d0
c dbgend
      end do

      if (ntest.ge.100) then
        write(luout,*) 'U*s^(-0.5):'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

c      call dgemm('n','t',ndim,ndim,ndim,
c     &           1d0,eigen_vec,ndim,
c     &               eigen_vec,ndim,
c     &           0d0,mat,ndim)


      return
      end
