*----------------------------------------------------------------------*
      subroutine set_final_reo(reo_info,xlines,xlines_tmp,
     &     op1op2,idxop,nvtx,nj)
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
     &     idxop,nj,nvtx
      integer(8), intent(in) ::
     &     op1op2(nj)
      integer(8), intent(inout), target ::
     &     xlines(max(nvtx,nj),nj),
     &     xlines_tmp(max(nvtx,nj),max(nvtx,nj))

      integer(8), pointer ::
     &     int12(:), xlines_pnt(:,:), xlines_tmp_pnt(:,:),
     &     scr(:), last_chance(:)

      logical, pointer ::
     &     assign_ok(:,:), scr_ok(:), merge_ok(:)
      integer, pointer ::
     &     must_assign(:)

      logical ::
     &     simple
      integer ::
     &     ivtx, idx, ij, occ_sh(ngastp,2), last,
     &     iocc_int(ngastp,2,max(nvtx,nj)), ivtxij(nj), nullij(0:nj+1),
     &     idx2, icnt, merge_ij, merge_vtx

      integer(8), external ::
     &     occ_overlap_p
      integer, external ::
     &     idxlist

      ! if xlines is already diagonal, do nothing
      ij_loop: do ij = 1, nj
        do ivtx = 1, max(nvtx,nj)
          if (ivtx.ne.ij.and.xlines(ivtx,ij).ne.0) exit ij_loop
        end do
        if (ij.eq.nj) return
      end do ij_loop

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_final_reo')
      end if
      xlines_pnt => xlines
      xlines_tmp_pnt => xlines_tmp

      allocate(int12(nvtx))
      ivtxij(1:nj) = 0
      nullij(0:nj+1) = 0

      ! sum lines: intermediate before reordering
      do ivtx = 1, nvtx
        int12(ivtx) = sum(xlines_pnt(ivtx,1:nj))
      end do
      iocc_int = 0
      call unpack_occ(iocc_int(1:ngastp,1:2,1:nvtx),int12,nvtx)

      ! if there are null-vertices, there is some flexibility in the final order
      ! we determine an order ivtxij with the smallest (?) number of reorderings
      ! at least the null-vertices are left alone and not shifted around ...
      do ij = 1, nj
        do ivtx = 1, nvtx
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
c        write(luout,'(a,5i2)') 'null-vector  : ',nullij(1:nj)
      end if

      ! get allowed assignments ivtx - ij 
      allocate(scr(max(nvtx,nj)),assign_ok(max(nvtx,nj),nj),
     &         must_assign(nvtx),last_chance(nvtx+1),
     &         scr_ok(max(nvtx,nj)),merge_ok(max(nvtx,nj)))
      scr = 0
      assign_ok = .true.
      scr_ok = .true.
      merge_ij = 0
      merge_vtx = 0
      must_assign = 0
      do ij = 1, nj
        do ivtx = 1, nvtx-1
          icnt = 0
          do idx = ivtx+1, nvtx
            scr(idx) = occ_overlap_p(xlines_pnt(ivtx,ij),
     &                                    xlines_pnt(idx,ij))
            if (scr(idx).ne.0) then
              icnt = icnt + 1
              do idx2 = ivtx+1, idx-1
                if (occ_overlap_p(scr(idx2),scr(idx)).ne.0)
     &             call quit(1,'set_final_reo','3-string resolution!')
              end do
            end if
          end do
          if (icnt.eq.0) cycle
          do idx = 1, max(nvtx,nj)
            if (icnt.eq.1) then
              if (ivtx.eq.idx) cycle
              scr_ok(idx) = scr(idx).ne.0
            else
              scr_ok(idx) = ivtx.eq.idx
            end if
          end do
          ! special merge required?
          simple = .false.
          do idx = 1, max(nvtx,nj)
            simple = simple.or.scr_ok(idx).and.assign_ok(idx,ij)
          end do
          if (simple) then
            do idx = 1, max(nvtx,nj)
              assign_ok(idx,ij) = assign_ok(idx,ij).and.scr_ok(idx)
            end do
          else
            if (merge_ij.eq.ij) then
              do idx = 1, max(nvtx,nj)
                merge_ok(idx) = merge_ok(idx).and.scr_ok(idx)
              end do
              if (.not.any(merge_ok(1:max(nvtx,nj))))
     &            call quit(1,'set_final_reo',
     &            'need additional special merge (1)')
            else if (merge_ij.ne.0) then
              call quit(1,'set_final_reo',
     &             'need additional special merge (2)')
            end if
            merge_ok(1:max(nvtx,nj)) = scr_ok(1:max(nvtx,nj))
            merge_ij = ij
          end if
          scr = 0
          scr_ok = .true.
        end do
      end do
      do ivtx = 1, nvtx
        do ij = 1, nj-1
          icnt = 0
          do idx = ij+1, nj
            scr(idx) = occ_overlap_p(xlines_pnt(ivtx,ij),
     &                                    xlines_pnt(ivtx,idx))
            if (scr(idx).ne.0) then
              icnt = icnt + 1
              do idx2 = ij+1, idx-1
                if (occ_overlap_p(scr(idx2),scr(idx)).ne.0)
     &             call quit(1,'set_final_reo','3-string resolution!')
              end do
            end if
          end do
          if (icnt.eq.0) cycle
          must_assign(ivtx) = 1
          do idx = 1, nj
            if (icnt.eq.1) then
              if (ij.eq.idx) cycle
              assign_ok(ivtx,idx) =assign_ok(ivtx,idx).and.scr(idx).ne.0
            else
              assign_ok(ivtx,idx) = assign_ok(ivtx,idx).and.ij.eq.idx
            end if
          end do
          scr = 0
        end do
      end do

      ! if special merge is required, a merged vtx may not be assigned
      ! to a target vtx. In case of only one merge vtx we can prevent it:
      if (merge_ij.ne.0) then
        icnt = 0
        do idx = 1, max(nvtx,nj)
          if (merge_ok(idx)) icnt = icnt + 1
        end do
        if (icnt.eq.1) then
          do idx = 1, max(nvtx,nj)
            if (merge_ok(idx)) then
              assign_ok(idx,1:nj) = .false.
              exit
            end if
          end do
        end if
      end if

      last_chance = 0
      last_chance(1) = nj+1
      idx = 0
      last = nj
      do ivtx = nvtx, 1, -1
        if (.not.must_assign(ivtx)) cycle
        idx = idx + 1
        do ij = last, 1, -1
          if (assign_ok(ivtx,ij)) then
            last_chance(idx+1) = ij
            last = ij
            exit
          end if
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'xlines:'
        do ivtx = 1,nvtx
          write(luout,'(6i10)') xlines_pnt(ivtx,1:nj)
        end do
        write(luout,*) 'assign_ok:'
        do ivtx = 1, max(nvtx,nj)
          write(luout,*) assign_ok(ivtx,1:nj)
        end do
        write(luout,'(a,6i2)') 'must_assign: ',must_assign
        if (merge_ij.ne.0) then
          write(luout,*) 'merge_ij: ',merge_ij
          write(luout,*) 'merge_ok: ',merge_ok
        end if
c        write(luout,'(a,6i2)') 'last_chance: ',last_chance(2:nvtx+1)
      end if

      ! check if possible reordering can exist
      do ij = 1, nj
        if (nullij(ij).eq.1) cycle
        if (.not.any(assign_ok(1:max(nvtx,nj),ij)))
     &      call quit(1,'set_final_reo','no possible reordering! (1)')
      end do
      do ivtx = 1, nvtx
        if (.not.must_assign(ivtx)) cycle
        if (.not.any(assign_ok(ivtx,1:nj)))
     &      call quit(1,'set_final_reo','no possible reordering! (2)')
      end do

      last = 0
      do ij = 1, nj
        if (ivtxij(ij).eq.0) cycle
        if (ivtxij(ij).le.last.or.ivtxij(ij).ne.ij.or.
     &      ivtxij(ij).gt.ij+sum(nullij(ij+1:nj))+max(nvtx,nj)-nj.or.
     &      ivtxij(ij).lt.ij-sum(nullij(1:ij-1)).or.
     &      (ij.ge.last_chance(sum(must_assign(ivtxij(ij):nvtx))+1).and.
     &       must_assign(ivtxij(ij)).eq.0)) then
          ivtxij(ij) = last + 1
          do ivtx = ivtxij(ij), nvtx
            if (xlines_pnt(ivtx,ij).ne.0.and.
     &          ivtx.gt.last.and.assign_ok(ivtx,ij).and.
     &          ivtx.le.ij+sum(nullij(ij+1:nj))+max(nvtx,nj)-nj.and.
     &          ivtx.ge.ij-sum(nullij(1:ij-1)).and.
     &          (ij.lt.last_chance(sum(must_assign(ivtx:nvtx))+1).or.
     &           must_assign(ivtx).eq.1)) then
              ivtxij(ij) = ivtx
c              exit
              if (ivtx.ge.ij) exit
            end if
          end do
        end if
        last = ivtxij(ij)
      end do

      ! find vtx which has to be merged (if necessary)
      if (merge_ij.ne.0) then
        do idx = 1, max(nvtx,nj)
          if (idxlist(idx,ivtxij,nj,1).le.0.and.merge_ok(idx)) then
            merge_vtx = idx
            if (must_assign(idx).eq.1) exit
          end if
        end do
      end if

      if (ntest.ge.100) then
        write(luout,'(a,5i2)') 'aim at order: ',ivtxij(1:nj)
        if (merge_ij.ne.0) write(luout,*) 'merge_vtx:',merge_vtx
      end if

      ! check if assignment has worked
      do ij = 1, nj
        if (nullij(ij).eq.1) cycle
        if (ivtxij(ij).le.0.or.ivtxij(ij).gt.max(nvtx,nj).or.
     &      .not.assign_ok(ivtxij(ij),ij))
     &      call quit(1,'set_final_reo','oops! messed up? (1)')
      end do
      do ivtx = 1, nvtx
        if (must_assign(ivtx).eq.0) cycle
        if (idxlist(ivtx,ivtxij,nj,1).le.0.and.ivtx.ne.merge_vtx)
     &      call quit(1,'set_final_reo','oops! messed up? (2)')
      end do

      ! set reo_info if reordering is necessary
      do ivtx = 1, nvtx
        do ij = 1, nj
          if (ivtx.eq.ivtxij(ij)) cycle
          if (xlines_pnt(ivtx,ij).ne.0) then
            if (.not.associated(reo_info%reo)) then
              reo_info%nreo = 0
              allocate(reo_info%reo(2*nj),
     &                 reo_info%nca_vtx(max(nvtx,nj)))
              reo_info%nvtx_contr = nvtx
              call set_nca_vtx(reo_info%nca_vtx,iocc_int,max(nvtx,nj))
            end if
            reo_info%nreo = reo_info%nreo+1
            idx = reo_info%nreo
            reo_info%reo(idx)%idxsuper = 1
            reo_info%reo(idx)%idxop_ori = idxop
            reo_info%reo(idx)%iblkop_ori = 1
            reo_info%reo(idx)%is_bc_result = .true.
            reo_info%reo(idx)%reo_before = .false.
            reo_info%reo(idx)%shift_i0 = ij.eq.merge_ij.and.
     &                                   ivtx.eq.merge_vtx
            reo_info%reo(idx)%to = ivtxij(ij)
            reo_info%reo(idx)%from = ivtx
            reo_info%reo(idx)%to_vtx   = ivtxij(ij)
            reo_info%reo(idx)%from_vtx = ivtx
            call unpack_occ(occ_sh,xlines_pnt(ivtx,ij),1)
            reo_info%reo(idx)%occ_shift = occ_sh
            if (ntest.ge.100) then
              if (reo_info%reo(idx)%shift_i0) then
                write(luout,'(a,i2,a,i2)') 'Reorder i0 from ',ivtx,
     &                                     ' to ',ivtxij(ij)
              else
                write(luout,'(a,i2,a,i2)') 'Reorder from ',ivtx,
     &                                     ' to ',ivtxij(ij)
              end if
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
        xlines_tmp_pnt(1:nvtx,1:nj) = xlines_pnt(1:nvtx,1:nj)
      else
        do ivtx = 1, nvtx
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
      deallocate(int12,assign_ok,scr,must_assign,last_chance,
     &           scr_ok,merge_ok)

      return
      end
