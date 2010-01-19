      subroutine topo_merge_vtxs2(ireo,nvtx_new,nvtx_bcres,
     &                           merge_sign,
     &                           topo,xlines,nvtx,nj,
     &                           svertex,isvtx1,isvtx2,
     &                           vtx_list,nlist)
      
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_operators.h'

      integer, intent(in) ::
     &     nvtx, nj, nlist,
     &     vtx_list(nlist),
     &     svertex(nvtx), isvtx1, isvtx2
      integer(8), intent(inout) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj)
      integer, intent(out) ::
     &     ireo(nvtx), nvtx_new, nvtx_bcres, merge_sign
      
      integer ::
     &     iord(nvtx)
      integer(8) ::
     &     topo0(nvtx,nvtx), xlines0(nvtx,nj)
      integer ::
     &     idx, jdx, kdx, jdxnd, ii, ivtx, ij,
     &     n_zero_vtx_res, n_zero_vtx, kdx_v, jdx_v,
     &     n_enclosed, isvtx_j, isvtx_k, ipass,
     &     iv_to, iv_fr
      logical ::
     &     merged(nvtx), reversed, right_sh

      integer ::
     &     occ_j(ngastp,2), occ_k(ngastp,2), occ_dum(ngastp,2),
     &     occ_fr(ngastp,2), occ_to(ngastp,2)

      logical, external ::
     &     zero_i8vec
      integer, external ::
     &     sign_hpvx, sign_merge

      ! check list
      do idx = 1, nlist
        do jdx = 1, idx-1
          if (vtx_list(jdx).ge.vtx_list(idx)) then
            write(luout,*) '>',vtx_list(1:nlist)
            call quit(1,'topo_merge_vtxs',
     &           'a unique, ascending list was expected')

          end if
        end do
      end do

      merge_sign = 1

      ! count number of external vertices that contain
      ! no lines (dummy vertices to keep same number of 
      ! joined vertices for all occupations of an operator,
      ! e.g. for densities ...)
      n_zero_vtx_res = 0
      do ij = 1, nj
        if (zero_i8vec(xlines(1,ij),nvtx,1))
     &       n_zero_vtx_res = n_zero_vtx_res+1
      end do
      n_zero_vtx = 0
      do ivtx = 1, nvtx
        if (zero_i8vec(xlines(ivtx,1),nj,nvtx))
     &       n_zero_vtx = n_zero_vtx+1
      end do
      
      do idx = 1, nvtx
        ireo(idx) = idx
        iord(idx) = idx
      end do

      merged(1:nvtx) = .false.
      nvtx_new = nvtx
      nvtx_bcres = nlist
      topo0 = topo
      xlines0 = xlines

      ! loop over list
      do idx = 1, nlist
        if (merged(idx)) cycle
c        ivtx1 = vtx_list(idx)
        ! find contiguos sequence of merge candidates:
        jdxnd = idx
        collect_list: do jdx = idx+1, nlist
          if (vtx_list(jdx)-vtx_list(jdx-1).gt.1) exit
          ! FIX -- check for self-contractions:
          do kdx = idx+1, jdx-1
            if (topo(kdx,jdx).gt.0) exit collect_list 
          end do
          jdxnd = jdx
        end do collect_list
c dbg
c        print *,'present list: ',idx,jdxnd
c dbg
        ! only one candidate? not much to do then
        if (jdxnd.eq.idx) cycle

        ! loop over pairs of vertices in contiguous list and try to merge
        jloop: do jdx = idx, jdxnd           
c dbg
c          print *,'jdx =',jdx
c          print *,'merged',merged
c          print *,'iord  ',iord
c          print *,'ireo  ',ireo
c dbg
          if (merged(jdx)) cycle jloop

          ! get supervertex number of vertex on which we merge
          ! we need to first carry out all merges with vertices
          ! having the same supervertex number (for sign-consistency)
          jdx_v = vtx_list(jdx)
          isvtx_j = svertex(jdx_v)

          do ipass = 1, 2

           do kdx = jdx+1, jdxnd
            if (merged(jdx)) cycle jloop
            if (merged(kdx)) cycle
            kdx_v = vtx_list(kdx)
            isvtx_k = svertex(kdx_v)

            if (ipass.eq.1.and.isvtx_k.ne.isvtx_j) cycle
            if (ipass.eq.2.and.isvtx_k.eq.isvtx_j) cycle

            if (may_merge(jdx_v,kdx_v,right_sh)) then
              if (.not.right_sh) then
                iv_fr = kdx_v ! shift from right
                iv_to = jdx_v ! to left
                merged(kdx) = .true.
              else
                iv_fr = jdx_v ! shift from left
                iv_to = kdx_v ! to right
                merged(jdx) = .true.
              end if
c              jdx_v = vtx_list(jdx)
c              kdx_v = vtx_list(kdx)

              ! reversed sequence of super-vertices?
              reversed = isvtx1.ne.isvtx2.and.
     &                   isvtx_j.eq.isvtx2 .and.
     &                   isvtx_k.eq.isvtx1
c dbg
c              print *,'isvtx1,isvtx2: ',isvtx1,isvtx2
c              print *,'svertex:       ',svertex(jdx_v),svertex(kdx_v)
c              print *,'  -> reversed = ',reversed
c              if (right_sh) print *,'Performing left -> right merge!'
c dbg

              ! hande sign change upon merge
              call topo_get_cnt_i0_j0(occ_dum,occ_j,occ_k,n_enclosed,
     &             jdx_v,kdx_v,
     &             topo,xlines,nvtx,nj)
c dbg
c        print *,'idx = ',jdx_v,kdx_v
c        print *,'0,I0,J0:'
c        call wrt_occ(6,occ_dum)
c        call wrt_occ(6,occ_j)
c        call wrt_occ(6,occ_k)
c        print *,'nenclosed = ',n_enclosed
c dbg

              if (.not.right_sh) then
                occ_fr = occ_k
                occ_to = occ_j
              else
                ! switch and transpose for correct approach
                occ_fr(1:ngastp,1) = occ_j(1:ngastp,2)
                occ_fr(1:ngastp,2) = occ_j(1:ngastp,1)
                occ_to(1:ngastp,1) = occ_k(1:ngastp,2)
                occ_to(1:ngastp,2) = occ_k(1:ngastp,1)
              end if

              ! get the sign for approaching 
              ! {I0C I0A}{}{J0C J0A} -> {I0C J0C I0A J0A}{}{}
              merge_sign = merge_sign *
     &             sign_merge(occ_to,occ_fr,n_enclosed,reversed) 
c dbg
c              print *,'after sign_merge: ',merge_sign
c dbg
              ! get the sign for HPVX reordering
              ! {I0C J0C I0A J0A} -> {IJC IJA}
              if (.not.reversed) then
                merge_sign = merge_sign *
     &             sign_hpvx(2,occ_j,.false.,occ_k,.false.)
              else
                merge_sign = merge_sign *
     &             sign_hpvx(2,occ_k,.false.,occ_j,.false.)
              end if
c dbg
c              print *,'after hpvx: ',merge_sign
c dbg

              ! also change position of already merged vtxs
              do ii = 1, nvtx
                if (ii.eq.iv_fr) cycle
                if (iord(ii).eq.iord(iv_fr)) then
                  iord(ii) = iord(iv_to)
                  ireo(ii) = ireo(iv_to)
                end if
              end do
              iord(iv_fr) = iord(iv_to)
              ireo(iv_fr) = ireo(iv_to)

c dbg
c              print *,'modified in iord: ',iv_fr
c              print *,'iord  ',iord
c              print *,'ireo  ',ireo
c dbg
              do ii = 1, nvtx
                if (ireo(ii).gt.iv_fr)
     &               iord(ii) = iord(ii)-1
              end do
c dbg
c              print *,'iord(final) ',iord
c dbg
              nvtx_new = nvtx_new-1
              nvtx_bcres = nvtx_bcres-1
              if (zero_i8vec(xlines(jdx_v,1),nj,nvtx) .or.
     &            zero_i8vec(xlines(kdx_v,1),nj,nvtx) )
     &             n_zero_vtx = n_zero_vtx-1

              ! add column topo(:,kdx) to topo(:,jdx)
              topo(1:nvtx,iv_to) = topo(1:nvtx,iv_to) 
     &                           + topo(1:nvtx,iv_fr)
              ! zero column topo(:,kdx)
              topo(1:nvtx,iv_fr) = 0
              ! add row topo(kdx,:) to topo(jdx,:)
              topo(iv_to,1:nvtx) = topo(iv_to,1:nvtx) 
     &                           + topo(iv_fr,1:nvtx)
              ! zero row topo(kdx,:)
              topo(iv_fr,1:nvtx) = 0

              ! add xlines(kdx,:) to xlines(jdx,:)
              xlines(iv_to,1:nj) = xlines(iv_to,1:nj) 
     &                           + xlines(iv_fr,1:nj)
              ! zero xlines(kdx,:)
              xlines(iv_fr,1:nj) = 0

            end if
           end do ! kdx
          end do ! ipass
        end do jloop ! jdx

      end do

      ireo = iord

      return

      contains

      logical function may_merge(ivtx1,ivtx2,right_shift)

      implicit none

      integer, intent(in) ::
     &     ivtx1, ivtx2
      logical, intent(out) ::
     &     right_shift

      integer ::
     &     icnt, ij, idx1, idx2, iel, jvtx, num1, num2,
     &     ivtx, iocc(ngastp,2), jocc(ngastp,2)
      logical ::
     &     horiz, verti

      integer ::
     &     cnt_line1(nvtx+nj), cnt_line2(nvtx+nj),
     &     ovl_line(nvtx+nj)
      integer(8) ::
     &     base
      
      integer, external ::
     &     imltlist, occ_p_el, idxmax, int8_expand
      logical, external ::
     &     zero_i8vec, zero_ivec

      right_shift = .false.
      base = pack_base

c dbg
c      print *,'may_merge>',ivtx1,ivtx2
c dbg
      ! very inital test: no self contraction between the two
      ! nodes to be merged:
      may_merge = topo(ivtx1,ivtx2).eq.0.and.topo(ivtx2,ivtx1).eq.0
      if (.not.may_merge) return
      ! check external lines
      ! preliminarly allow merge only if ...
      ! a) the number of vertices is still larger than the
      !    number of final result vertices
      may_merge = nvtx_new.gt.nj
c dbg
c      print *,'a) ',may_merge
c dbg
      if (.not.may_merge) return

      ! b) either vertex has no external lines or ...
      may_merge = zero_i8vec(xlines(ivtx1,1),nj,nvtx).or.
     &            zero_i8vec(xlines(ivtx2,1),nj,nvtx)

      ! (unless we need a certain number of zero vertices)
      if (may_merge.and.n_zero_vtx_res.gt.0.and.
     &                  n_zero_vtx.le.n_zero_vtx_res) then
        may_merge = .false.
        return
      end if

      ! do not shift non-zero vertex to zero vertex
      right_shift = may_merge.and.
     &              .not.zero_i8vec(xlines(ivtx2,1),nj,nvtx)

c dbg
c      print *,'b) ',may_merge
c dbg
      ! b) both vertices contribute to the same *single* vertex
      if (.not.may_merge) then
        icnt = 0
        do ij = 1, nj
          if (xlines(ivtx1,ij).gt.0) then
            icnt = icnt+1
            idx1 = ij
          end if
        end do
        if (icnt.eq.1) then
          icnt = 0
          do ij = 1, nj
            if (xlines(ivtx2,ij).gt.0) then
              icnt = icnt+1
              idx2 = ij
            end if
          end do
          may_merge = icnt.eq.1.and.idx1.eq.idx2
        end if
c dbg
c        print *,'c) ',may_merge
c dbg

        ! do right-shift if advantageous vertex order is obtained
        ! --> less reordering steps later on
        if (may_merge.and.ivtx2-ivtx1.gt.1) then
          ! right-shift if enclosed vtxs belong to higher result vtx ...
          do ij = idx1-1, 1, -1
            if (sum(xlines(ivtx1+1:ivtx2-1,ij)).gt.0) then
              right_shift = .true.
              cycle
            end if
          end do
          ! but not if there are some belonging to even lower vertices
          ! --> not trivially clear if a right-shift would be better
          if (right_shift) then
            do ij = idx1+1, nj
              if (sum(xlines(ivtx1+1:ivtx2-1,ij)).gt.0) then
                right_shift = .false.
                cycle
              end if
            end do
          end if
        end if

      end if

      ! exit, if these tests failed
      if (.not.may_merge) return

c dbg
c      print *,'final test'
c dbg
      ! loop over elements of occupation vector
      do iel = 1, 2*ngastp
c dbg
c        print *,'iel = ',iel
c dbg
        ! get contraction occupations for present element ...
        do jvtx = 1, nvtx
          cnt_line1(jvtx) = occ_p_el(topo(ivtx1,jvtx),iel)
          cnt_line2(jvtx) = occ_p_el(topo(ivtx2,jvtx),iel)
        end do
        ! ... including the external line info
        ij = 0
        do jvtx = nvtx+1, nvtx+nj
          ij = ij+1
          cnt_line1(jvtx) = occ_p_el(xlines(ivtx1,ij),iel)
          cnt_line2(jvtx) = occ_p_el(xlines(ivtx2,ij),iel)
        end do
c dbg
c        print *,'l1:',cnt_line1(1:nvtx+nj)
c        print *,'l2:',cnt_line2(1:nvtx+nj)
c dbg
        ! we may merge, if only either vertex contributes
        ! internal or external lines
        may_merge = zero_ivec(cnt_line1,nvtx+nj).or.
     &          zero_ivec(cnt_line2,nvtx+nj)

        ! ... else ...
        if (.not.may_merge) then
          ! check for number of entries
          num1 = nvtx+nj-imltlist(0,cnt_line1,nvtx+nj,1)
          num2 = nvtx+nj-imltlist(0,cnt_line2,nvtx+nj,1)
          ! at most 1 entry per node ...
          if (num1.eq.1.and.num2.eq.1) then
            ! ... in the same column!
            num1 = idxmax(cnt_line1,nvtx+nj,1)
            num2 = idxmax(cnt_line2,nvtx+nj,1)
            may_merge = num1.eq.num2
cmh
            ! exception only if merging vtxs from different supervtxs
            if (may_merge) then
              ivtx_loop: do ivtx = 1, nvtx
                if (ireo(ivtx).ne.ivtx1) cycle ivtx_loop
                do jvtx = 1, nvtx
                  if (ireo(jvtx).ne.ivtx2) cycle
                  if (svertex(ivtx).ne.svertex(jvtx)) cycle
                  iocc = 0
                  jocc = 0
                  if (num1.le.nvtx) then
                   icnt = int8_expand(topo0(ivtx,num1),base,iocc)
                   icnt = int8_expand(topo0(jvtx,num1),base,jocc)
                  else
                   icnt = int8_expand(xlines0(ivtx,num1-nvtx),base,iocc)
                   icnt = int8_expand(xlines0(jvtx,num1-nvtx),base,jocc)
                  end if
                  if (iocc_nonzero(iocc_overlap(iocc,.false.,
     &                jocc,.false.))) then
                    may_merge = .false.
                    exit ivtx_loop
                  end if
                end do
              end do ivtx_loop
            end if
cmhend
          else
            may_merge = .false.
          end if
        end if
        if (.not.may_merge) exit

          ! this seems not correct:
c        ! the occupations may either be horizontally only
c        ! or vertically only distributed
c        horiz = imltlist(0,cnt_line1,nvtx+nj,1).lt.nvtx+nj-1
c     &      .or.imltlist(0,cnt_line2,nvtx+nj,1).lt.nvtx+nj-1
c        call list_ovl(ovl_line,cnt_line1,cnt_line2,nvtx+nj,1)
c        verti = imltlist(0,ovl_line,nvtx+nj,1).lt.nvtx+nj
c        ! if both is the case, we cannot symmetrize (yet)
cc dbg
c        print *,'horiz, verti: ',horiz,verti
cc dbg
c        if (horiz.and.verti) then
c          may_merge = .false.
c          exit
c        end if
c dbg
c        print *,'OK'
c dbg

      end do
c dbg
c      if (.not.may_merge) print *,'not OK'
c dbg

      return
      end function

      end
