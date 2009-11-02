*----------------------------------------------------------------------*
      subroutine set_final_reo(reo_info,xlines,xlines_tmp,
     &     op1op2,idxop,nj)
*----------------------------------------------------------------------*
*     After a final contraction, xlines should be a diagonal matrix.
*     if this is not the case, additional reordering steps are required.
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
     &     ivtx, idx, ij, occ_sh(ngastp,2),
     &     iocc_int(ngastp,2,nj)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_final_reo')
      end if
      xlines_pnt => xlines
      xlines_tmp_pnt => xlines_tmp

      allocate(int12(nj))

      ! sum lines: intermediate before reordering
      do ivtx = 1, nj
        int12(ivtx) = sum(xlines_pnt(ivtx,1:nj))
      end do
      call unpack_occ(iocc_int,int12,nj)

      ! set reo_info if reordering is necessary
      do ivtx = 1, nj
        do ij = 1, nj
          if (ivtx.eq.ij) cycle
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
            reo_info%reo(idx)%to = ij
            reo_info%reo(idx)%from = ivtx
            reo_info%reo(idx)%to_vtx   = ij
            reo_info%reo(idx)%from_vtx = ivtx
            call unpack_occ(occ_sh,xlines_pnt(ivtx,ij),1)
            reo_info%reo(idx)%occ_shift = occ_sh
            if (ntest.ge.100) then
              write(luout,'(a,i2,a,i2)') 'Reorder from ',ivtx,' to ',ij
              call wrt_occ(luout,occ_sh)
            end if

            ! xlines: sort corresponding to end result
            xlines_pnt(ij,ij) = xlines_pnt(ij,ij) 
     &                            + xlines_pnt(ivtx,ij)
            ! xlines_tmp: sort corresponding to temporary intermediate
            xlines_tmp_pnt(ivtx,ivtx) = xlines_tmp_pnt(ivtx,ivtx)
     &                            + xlines_tmp_pnt(ivtx,ij)
            xlines_pnt(ivtx,ij) = 0
            xlines_tmp_pnt(ivtx,ij) = 0
          end if
        end do
      end do

      ! check if final result matches given final operator
      do ij = 1, nj
        if (xlines_pnt(ij,ij).ne.op1op2(ij)) then
          write(luout,'(a,5i9.8)') 'xlines_pnt(i,i): ',
     &               (xlines_pnt(ivtx,ivtx), ivtx=1,nj)
          write(luout,'(a,5i9.8)') 'op1op2(i)      : ',
     &               op1op2(1:nj)
          call quit(1,'set_final_reo','should be final contraction')
        end if
      end do
      deallocate(int12)

      return
      end
