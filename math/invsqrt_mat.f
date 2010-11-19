*----------------------------------------------------------------------*
      subroutine invsqrt_mat(ndim,mat,mat2,half,icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins)
*----------------------------------------------------------------------*
*     half = true: calculates U*mat^(-0.5) using MAT = U*mat*U^+
*     half = false: calculates both U*mat^(-0.5) and U*1s*U^+
*                   which is a projector matrix for eliminating
*                   linear redundancies
*
*     matthias, dec 2009
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'routes.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 0
      real(8), parameter ::
     &     warn_sv = 1d-12!5  ! give a warning for small singular values
      integer, intent(in) ::
     &     ndim
      integer, intent(inout) ::
     &     icnt_sv, icnt_sv0, bins(17)
      real(8), intent(inout), target ::
     &     mat(ndim,ndim), mat2(ndim,ndim), xmax, xmin
      logical, intent(in) ::
     &     half
      real(8) ::
     &     singval(ndim),wrk(max(1024,ndim**2)),dum1,dum2,expo
      real(8), pointer ::
     &     mat_tmp(:,:)

      type(filinf) ::
     &     ffsv
      integer ::
     &     nrot, idx, lwrk, info, idx2, luinp, iostatus
      logical ::
     &     l_exist, sv_above, read_file

      if (ndim.eq.0) return

      lwrk=max(1024,ndim**2)
      info = 0

      if (sv_fix) then
        inquire(file='SINGVALS',exist=l_exist)
        call file_init(ffsv,'SINGVALS',ftyp_sq_frm,0)
        call file_open(ffsv)
        luinp = ffsv%unit
        read_file = l_exist
        if (l_exist) then
          do idx = 1, icnt_sv + 1
            read(luinp,*,iostat=iostatus) idx2, sv_above
          end do
          read_file = iostatus.ge.0
          rewind luinp
          do idx = 1, icnt_sv
            read(luinp,*) idx2, sv_above
          end do
        end if
      end if

      if (ntest.ge.100) then
        write(luout,'(x,a)') '-------------------'
        write(luout,'(x,a)') 'invsqrt_mat at work'
        write(luout,'(x,a)') '-------------------'
        write(luout,*) 'input S matrix:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      ! check if S is symmetric
      do idx = 2, ndim
        do idx2 = 1,idx-1
          if (abs(mat(idx,idx2)-mat(idx2,idx)).gt.1d-12)
     &        call quit(1,'invsqrt_mat','S must be symmetric!')
        end do
      end do

      ! calculate U and singular values:
      call dgesvd('O','N',ndim,ndim,
     &     mat,ndim,singval,
     &     dum1,1,dum2,1,
     &     wrk,lwrk,info)

      if (info.ne.0) then
        write(luout,*) 'WARNING in invsqrt_mat: SVD in trouble'
      end if

      if (ntest.ge.10) then
        if (ntest.ge.100) write(luout,*) 'singular values:'
        do idx = 1, ndim
          write(luout,'(x,a,i8,x,f24.12)') 'sv #',icnt_sv+idx,
     &                                     singval(idx)
        end do
      end if
      if (ntest.ge.100) then
        write(luout,*) 'eigenvector matrix U:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

c      icnt_sv = icnt_sv + ndim

c      if (half) then
c        mat_tmp => mat
       expo = -0.5d0
c      else
       if (.not.half) allocate(mat_tmp(ndim,ndim))
cmh        expo = -0.25d0
c        expo = 0d0  ! to obtain matrix to project out lin. dep.
c      end if

      ! square root of (pseudo) inverse:
      ! A^(-1/2) = U D^(-1/2) U^+ = U D^(-1/4) [U D^(-1/4)]^+
      do idx = 1, ndim
        icnt_sv = icnt_sv + 1
        sv_above = .false.
        if (sv_fix.and.read_file) read(luinp,*) idx2, sv_above

        ! binning of singular value
        dum1 = singval(idx)
        do idx2 = 1, 17
          if (dum1.ge.1d0) then
            bins(idx2) = bins(idx2) + 1
            exit
          end if
          dum1 = dum1*10d0
          if (idx2.eq.17) bins(idx2) = bins(idx2) + 1
        end do

        if (.not.(sv_fix.and.read_file)
     &      .and.singval(idx).gt.sv_thresh.or.
     &      sv_above) then
          sv_above = .true.
          if (abs(singval(idx)).lt.abs(xmin)) xmin = singval(idx)
          if (singval(idx).lt.warn_sv)
     &         call warn('invsqrt_mat','small singular values!')
          if (.not.half) mat_tmp(1:ndim,idx) = mat(1:ndim,idx)
          mat(1:ndim,idx) = mat(1:ndim,idx)
     &                         * (singval(idx)**expo)
        else
          sv_above = .false.
          icnt_sv0 = icnt_sv0 + 1
          if (abs(singval(idx)).gt.abs(xmax)) xmax = singval(idx)
          mat(1:ndim,idx) = 0d0
          if (.not.half) mat_tmp(1:ndim,idx) = 0d0
        end if
        if (sv_fix.and..not.read_file) write(luinp,*) icnt_sv, sv_above
      end do

      if (ntest.ge.100) then
        write(luout,'(a,f5.2,a)') 'U*s^(',expo,')'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      if (.not.half) then
        call dgemm('n','t',ndim,ndim,ndim,
     &             1d0,mat_tmp,ndim,
     &                 mat_tmp,ndim,
     &             0d0,mat2,ndim)
        deallocate(mat_tmp)
        if (ntest.ge.100) then
          write(luout,*) 'U*1s*U^+ projector:'
          call wrtmat2(mat2,ndim,ndim,ndim,ndim)
        end if
      end if
c dbg  can be used to set 1
c      mat(1:ndim,1:ndim) = 0d0
c      do idx = 1, ndim
c        mat(idx,idx) = 1d0
c      end do
c dbgend

      if (sv_fix) call file_close_keep(ffsv)

      return
      end
