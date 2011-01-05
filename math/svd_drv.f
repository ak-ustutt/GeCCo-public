*----------------------------------------------------------------------*
      subroutine svd_drv(ndim,mat,singval)
*----------------------------------------------------------------------*
*     split matrix into non-coupling blocks and perform the svd
*     on these individual blocks
*     matrix should be symmetric
*
*     matthias, dec 2010
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'routes.h'
      include 'def_filinf.h'

      integer, parameter ::
     &     ntest = 00
      real(8), parameter ::
     &     coupl_thresh = 1d-14 ! threshold for coupling elements
      integer, intent(in) ::
     &     ndim
      real(8), intent(out) ::
     &     singval(ndim)
      real(8), intent(inout) ::
     &     mat(ndim,ndim)
      real(8) ::
     &     wrk(max(1024,ndim**2)),dum1,dum2,umat(ndim,ndim)
      real(8), pointer ::
     &     mat_tmp(:,:), sing_tmp(:)
      integer ::
     &     cur_els(ndim), nels, iels, iidim
      logical ::
     &     decomposed(ndim)

      integer ::
     &     nrot, idx, lwrk, info, idx2, luinp, iostatus
      integer, external ::
     &     idxlist

      umat = 0d0
      singval = 0d0
      decomposed = .false.

      do iidim = 1, ndim
        if (decomposed(iidim)) cycle

        ! find all coupling elements
        nels = 1
        cur_els(1) = iidim
        iels = 1
        ! find all elements which couple with elements on list
        do while(iels.le.nels)
          do idx2 = 1, ndim
            ! skip if already on list
            if (idxlist(idx2,cur_els,nels,1).gt.0) cycle
            if (abs(mat(cur_els(iels),idx2)).gt.coupl_thresh) then
              nels = nels + 1
              cur_els(nels) = idx2
            end if
          end do
          decomposed(cur_els(iels)) = .true.
          iels = iels + 1
        end do
        ! zero element? do nothing
        if (nels.eq.1.and.mat(cur_els(1),cur_els(1)).lt.coupl_thresh)
     &       cycle
c dbg
c        print *,'cur_els: ',cur_els(1:nels)
c dbgend

        ! assemble block
        allocate(mat_tmp(nels,nels),sing_tmp(nels))
        do idx = 1, nels
          do idx2 = 1, nels
            mat_tmp(idx,idx2) = mat(cur_els(idx),cur_els(idx2))
          end do
        end do

        ! call singular value decomposer
        lwrk=max(1024,nels**2)
        info = 0
        call dgesvd('O','N',nels,nels,
     &       mat_tmp,nels,sing_tmp,
     &       dum1,1,dum2,1,
     &       wrk,lwrk,info)
        if (info.ne.0) then
          write(luout,*) 'WARNING in svd_drv: SVD in trouble'
        end if

        ! copy singular values and eigenvector
        do idx = 1, nels
          singval(cur_els(idx)) = sing_tmp(idx)
          do idx2 = 1, nels
            umat(cur_els(idx),cur_els(idx2)) = mat_tmp(idx,idx2)
          end do
        end do

        deallocate(mat_tmp,sing_tmp)
      end do

      mat = umat

      return
      end
