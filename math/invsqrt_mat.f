*----------------------------------------------------------------------*
      subroutine invsqrt_mat(ndim,mat,half)
*----------------------------------------------------------------------*
*     half = true: calculates U*mat^(-0.5) using MAT = U*mat*U^+
*     half = false: calculates U*mat^(-0.5)*U^+
*     mat must be quadratic and symmetric (not checked so far)
*
*     matthias, dec 2009
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 0
      real(8), parameter ::
     &     min_sv = 1d-10, ! singular value threshold for calc. of pseudo-inv.
     &     warn_sv = 1d-5  ! give a warning for small singular values
      integer, intent(in) ::
     &     ndim
      real(8), intent(inout), target ::
     &     mat(ndim,ndim)
      logical, intent(in) ::
     &     half
      real(8) ::
     &     singval(ndim),wrk(max(1024,ndim**2)),dum1,dum2,expo
      real(8), pointer ::
     &     mat_tmp(:,:)

      integer ::
     &     nrot, idx, lwrk, info

      if (ndim.eq.0) return

      lwrk=max(1024,ndim**2)
      info = 0

      if (ntest.ge.100) then
        write(luout,'(x,a)') '-------------------'
        write(luout,'(x,a)') 'invsqrt_mat at work'
        write(luout,'(x,a)') '-------------------'
        write(luout,*) 'input S matrix:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      ! calculate U and singular values:
      call dgesvd('O','N',ndim,ndim,
     &     mat,ndim,singval,
     &     dum1,1,dum2,1,
     &     wrk,lwrk,info)

      if (info.ne.0) then
        write(luout,*) 'WARNING in invsqrt_mat: SVD in trouble'
      end if

      if (ntest.ge.10) write(luout,*) 'singular values s: ',singval
      if (ntest.ge.100) then
        write(luout,*) 'eigenvector matrix U:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      if (half) then
        mat_tmp => mat
        expo = -0.5d0
      else
        allocate(mat_tmp(ndim,ndim))
        expo = -0.25d0
      end if

      ! square root of (pseudo) inverse:
      ! A^(-1/2) = U D^(-1/2) U^+ = U D^(-1/4) [U D^(-1/4)]^+
      do idx = 1, ndim
        if (singval(idx).gt.min_sv) then
          if (singval(idx).lt.warn_sv)
     &         call warn('invsqrt_mat','small singular values!')
          mat_tmp(1:ndim,idx) = mat(1:ndim,idx)
     &                         * (singval(idx)**expo)
        else
          mat_tmp(1:ndim,idx) = 0d0
        end if
c dbg   can be (mis)used to set a unit operator
c        mat(1:ndim,idx) = 0d0
c        mat(idx,idx)  =1d0
c dbgend
      end do

      if (ntest.ge.100) then
        write(luout,'(a,f5.2,a)') 'U*s^(',expo,')'
        call wrtmat2(mat_tmp,ndim,ndim,ndim,ndim)
      end if

      if (.not.half) then
        call dgemm('n','t',ndim,ndim,ndim,
     &             1d0,mat_tmp,ndim,
     &                 mat_tmp,ndim,
     &             0d0,mat,ndim)
        deallocate(mat_tmp)
        if (ntest.ge.100) then
          write(luout,*) 'U*s^(-0.5)*U^+'
          call wrtmat2(mat,ndim,ndim,ndim,ndim)
        end if
      end if


      return
      end
