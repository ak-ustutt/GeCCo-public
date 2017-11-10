*----------------------------------------------------------------------*
      subroutine mat_svd_traf(ndim,mat,mat2,
     &                       icnt_sv,icnt_sv0,
     &                       xmax,xmin,bins)
*----------------------------------------------------------------------*
*     Perform an SVD and similarity transform the matrix
*     by the right-hand unitary matrix.
*     On output, mat2 contains the transformed matrix
*     and mat the unitary matrix
*
*     matthias, nov 2013
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'routes.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00
      integer, intent(in) ::
     &     ndim
      integer, intent(inout) ::
     &     icnt_sv, icnt_sv0, bins(17)
      real(8), intent(inout), target ::
     &     mat(ndim,ndim), mat2(ndim,ndim), xmax, xmin
      real(8) ::
     &     dum1, cur_sv, singval(ndim), vt(ndim,ndim)
      real(8), pointer::
     &     wrk(:)

      type(filinf) ::
     &     ffsv
      integer ::
     &     idx, info, idx2, luinp, iostatus, lwrk, ifree
      logical ::
     &     l_exist, sv_above, read_file, symmetrize
      
      if (ndim.eq.0) return
      
      ifree = mem_setmark('mat_svd_traf')

      info = 0
      lwrk = max(1024,ndim**2)
      ifree = mem_alloc_real(wrk, lwrk,'wrk')
      
      
      
      if (jac_fix) then
        inquire(file='SINGVALS2',exist=l_exist)
        call file_init(ffsv,'SINGVALS2',ftyp_sq_frm,0)
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
        write(lulog,'(x,a)') '--------------------'
        write(lulog,'(x,a)') 'mat_svd_traf at work'
        write(lulog,'(x,a)') '--------------------'
        write(lulog,*) 'input matrix:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

      ! for testing: but as already noted by MH, this does not
      !   lead to improvements, on the contrary actually ...  AK
      symmetrize = .false.  ! so leave .false.
      if (symmetrize) then
        ! symmetrize
        do idx = 1, ndim
          do idx2 = 1,idx-1
            mat(idx2,idx) = 0.5d0 *
     &                     (mat(idx2,idx) + mat(idx,idx2))
            mat(idx,idx2) = mat(idx2,idx)
          end do
        end do
        if (ntest.ge.100) then
          write(lulog,*) 'symmetrized matrix:'
          call wrtmat2(mat,ndim,ndim,ndim,ndim)
        end if
      end if

      ! copy matrix and perform SVD
      mat2 = mat
      call dgesvd('N','S',ndim,ndim,
     &     mat,ndim,singval,
     &     dum1,1,vt,ndim,
     &     wrk,lwrk,info)
      if (info.ne.0) call quit(1,'mat_svd_traf','SVD in trouble.')
      ! transpose V^T:
      do idx = 1, ndim
        mat(1:ndim,idx) = vt(idx,1:ndim)
      end do

c      ! symmetrize
c      do idx = 1, ndim
c        do idx2 = 1,ndim
c          mat2(idx2,idx) = mat(idx,idx2)
c        end do
c      end do
c      do idx = 1, ndim
c        do idx2 = 1,ndim
c          mat2(idx2,idx) = 0.5d0 *
c     &                    (mat(idx2,idx) + mat2(idx2,idx))
c        end do
c      end do
c      if (ntest.ge.100) then
c        write(lulog,*) 'symmetrized matrix:'
c        call wrtmat2(mat2,ndim,ndim,ndim,ndim)
c      end if
c
c      ! calculate U and eigenvalues:
c      call rs(ndim,ndim,mat2,eigr,1,mat,scr1,scr2,info)
c
c      if (info.ne.0) then
c        call warn('mat_svd_traf','Some problem in rs.f')
c      end if

      if (ntest.ge.100) then
        write(lulog,*) 'unitary matrix V:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
        write(lulog,*) 'singular values:'
      end if

      ! exclude eigenvalues below threshold or as given by SINGVALS2
      do idx = 1, ndim
        icnt_sv = icnt_sv + 1
        sv_above = .false.
        if (jac_fix.and.read_file) read(luinp,*) idx2, sv_above

        ! binning of eigenvalue
        cur_sv = singval(idx)
        dum1 = cur_sv
        do idx2 = 1, 17 !9
          if (dum1.ge.1d0) then
            bins(idx2) = bins(idx2) + 1
            exit
c          else if (dum1.le.-1d0) then
c            bins(17-idx2) = bins(17-idx2) + 1
c            exit
          end if
          dum1 = dum1*10d0
          if (idx2.eq.17) bins(idx2) = bins(idx2) + 1
c          if (idx2.eq.9) bins(idx2) = bins(idx2) + 1
        end do

        if (.not.(jac_fix.and.read_file).and.
     &      singval(idx).gt.jac_thresh
     &      .or.sv_above) then
          sv_above = .true.
          if (abs(singval(idx)).lt.abs(xmin)) xmin = singval(idx)
        else
          sv_above = .false.
          icnt_sv0 = icnt_sv0 + 1
          if (abs(singval(idx)).gt.abs(xmax)) xmax = singval(idx)
          mat(1:ndim,idx) = 0d0
          singval(idx) = 0d0
        end if
        if (ntest.ge.10) write(lulog,'(x,a,i8,x,f24.12,3x,L1)') 'ev #',
     &                 icnt_sv,cur_sv,sv_above
        if (jac_fix.and..not.read_file) write(luinp,*) icnt_sv, sv_above
      end do

      if (ntest.ge.100) then
        write(lulog,'(a,f5.2,a)') 'truncated matrix V:'
        call wrtmat2(mat,ndim,ndim,ndim,ndim)
      end if

c      ! write diagonalized matrix
c      mat2(1:ndim,1:ndim) = 0d0
c      do idx = 1, ndim
c        mat2(idx,idx) = eigr(idx)
c      end do

      ! explicitly compute similarity transformed matrix
      call dgemm('t','n',ndim,ndim,ndim,
     &            1d0,mat,ndim,
     &            mat2,ndim,
     &            0d0,vt,ndim)
      call dgemm('n','n',ndim,ndim,ndim,
     &            1d0,vt,ndim,
     &            mat,ndim,
     &            0d0,mat2,ndim)

      if (ntest.ge.100) then
        write(lulog,'(a,f5.2,a)') 'transformed matrix:'
        call wrtmat2(mat2,ndim,ndim,ndim,ndim)
      end if

      if (jac_fix) call file_close_keep(ffsv)
      
      ifree = mem_flushmark()

      return
      end
