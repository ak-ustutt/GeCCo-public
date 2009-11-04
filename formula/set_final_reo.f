*----------------------------------------------------------------------*
      subroutine set_final_reo(reo_info,xlines,xlines_tmp,
     &     op1op2,idxop,nj)
*----------------------------------------------------------------------*
*     After a final contraction, xlines should be a diagonal matrix.
*     if this is not the case, additional reordering steps are required
*     (except when the only distortion from a diagonal matrix is caused
*     by null-vertices).
*     Given the result operator in packed form, the index of the
*     required intermediate operator (before reordering), and xlines,
*     this routine sets reo_info accordingly and updates xlines.
*     In xlines_tmp, the reordering is done in a way to match the
*     occupation of the temporary intermediate that will be created
*     in course of the contraction.
*
*     matthias, fall 2009
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'

      integer, parameter ::
     &     ntest = 00

      type(reorder_info), intent(inout) ::
     &     reo_info
      integer, intent(in) ::
     &     idxop,nj
      integer(8), intent(in) ::
     &     op1op2(nj)
      integer(8), intent(inout), target ::
     &     xlines(nj,nj), xlines_tmp(nj,nj)

      integer(8), pointer ::
     &     int12(:), xlines_pnt(:,:), xlines_tmp_pnt(:,:)

      integer ::
     &     ivtx, idx, ij, occ_sh(ngastp,2), last,
     &     iocc_int(ngastp,2,nj), ivtxij(nj), nullij(0:nj+1)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_final_reo')
      end if
      xlines_pnt => xlines
      xlines_tmp_pnt => xlines_tmp

      allocate(int12(nj))
      ivtxij(1:nj) = 0
      nullij(0:nj+1) = 0

      ! sum lines: intermediate before reordering
      do ivtx = 1, nj
        int12(ivtx) = sum(xlines_pnt(ivtx,1:nj))
      end do
      call unpack_occ(iocc_int,int12,nj)

      ! if there are null-vertices, there is some flexibility in the final order
      ! we determine an order ivtxij with the smallest (?) number of reorderings
      ! at least the null-vertices are left alone and not shifted around ...
      do ij = 1, nj
        do ivtx = 1, nj
          if (xlines_pnt(ivtx,ij).ne.0) then
            ivtxij(ij) = ivtx
            exit
          end if
        end do
      end do
      do ij = 1, nj
        if (ivtxij(ij).eq.0) nullij(ij) = 1
      end do
      if (ntest.ge.100) then
        write(luout,'(a,5i2)') 'initial order: ',ivtxij(1:nj)
        write(luout,'(a,5i2)') 'null-vector  : ',nullij(1:nj)
      end if

      last = 0
      do ij = 1, nj
        if (ivtxij(ij).eq.0) cycle
        if (ivtxij(ij).le.last.or.
     &      ivtxij(ij).gt.ij+sum(nullij(ij+1:nj)).or.
     &      ivtxij(ij).lt.ij-sum(nullij(1:ij-1))) then
          ivtxij(ij) = last + 1
          do ivtx = ivtxij(ij), nj
            if (xlines_pnt(ivtx,ij).ne.0.and.
     &          ivtx.gt.last.and.
     &          ivtx.le.ij+sum(nullij(ivtx+1:nj)).and.
     &          ivtx.ge.ij-sum(nullij(1:ivtx-1))) then
              ivtxij(ij) = ivtx
              exit
            end if
          end do
        end if
        last = ivtxij(ij)
      end do
      if (ntest.ge.100) then
        write(luout,'(a,5i2)') 'aim at order: ',ivtxij(1:nj)
      end if

      ! set reo_info if reordering is necessary
      do ivtx = 1, nj
        do ij = 1, nj
          if (ivtx.eq.ivtxij(ij)) cycle
          if (xlines_pnt(ivtx,ij).ne.0) then
            if (.not.associated(reo_info%reo)) then
              reo_info%nreo = 0
              allocate(reo_info%reo(2*nj),
     &                 reo_info%nca_vtx(nj))
              reo_info%nvtx_contr = nj
              call set_nca_vtx(reo_info%nca_vtx,iocc_int,nj)
            end if
            reo_info%nreo = reo_info%nreo+1
            idx = reo_info%nreo
            reo_info%reo(idx)%idxsuper = 1
            reo_info%reo(idx)%idxop_ori = idxop
            reo_info%reo(idx)%iblkop_ori = 1
            reo_info%reo(idx)%is_bc_result = .true.
            reo_info%reo(idx)%reo_before = .false.
            reo_info%reo(idx)%to = ivtxij(ij)
            reo_info%reo(idx)%from = ivtx
            reo_info%reo(idx)%to_vtx   = ivtxij(ij)
            reo_info%reo(idx)%from_vtx = ivtx
            call unpack_occ(occ_sh,xlines_pnt(ivtx,ij),1)
            reo_info%reo(idx)%occ_shift = occ_sh
            if (ntest.ge.100) then
              write(luout,'(a,i2,a,i2)') 'Reorder from ',ivtx,
     &                                   ' to ',ivtxij(ij)
              call wrt_occ(luout,occ_sh)
            end if

            ! xlines: sort corresponding to end result
            xlines_pnt(ivtxij(ij),ij) = xlines_pnt(ivtxij(ij),ij) 
     &                            + xlines_pnt(ivtx,ij)
            xlines_pnt(ivtx,ij) = 0
          end if
        end do
      end do

      ! sort xlines_tmp corresponding to temporary intermediate
      if (reo_info%nreo.eq.0) then
        xlines_tmp_pnt = xlines_pnt
      else
        do ivtx = 1, nj
          do ij = 1, nj
            if (ivtx.eq.ij) cycle
            if (xlines_tmp_pnt(ivtx,ij).ne.0) then
              xlines_tmp_pnt(ivtx,ivtx) = xlines_tmp_pnt(ivtx,ivtx)
     &                              + xlines_tmp_pnt(ivtx,ij)
              xlines_tmp_pnt(ivtx,ij) = 0
            end if
          end do
        end do
      end if

      ! check if final result matches given final operator
      do ij = 1, nj
        if (.not.(op1op2(ij).eq.0.and.ivtxij(ij).eq.0).and.
     &      xlines_pnt(ivtxij(ij),ij).ne.op1op2(ij)) then
          write(luout,'(a,5i9.8)') 'xlines_pnt(i,i): ',
     &               (xlines_pnt(ivtxij(ivtx),ivtx), ivtx=1,nj)
          write(luout,'(a,5i9.8)') 'op1op2(i)      : ',
     &               op1op2(1:nj)
          call quit(1,'set_final_reo','should be final contraction')
        end if
      end do
      deallocate(int12)

      return
      end
