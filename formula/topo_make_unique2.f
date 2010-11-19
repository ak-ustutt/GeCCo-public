*----------------------------------------------------------------------*
      subroutine topo_make_unique2(ireo,vtx,svtx,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     (unrestricted) sort of vertices, topo matrix, and xline matrix
*     to ascending sequence as needed for comparing terms
*     not intended for generating the actual unique sequence for
*     the representation (which needs to consider restricted sort and
*     thus is much more expensive)
*     
*     the input vertices are reordered
*
*     ireo(1:nvtx) contains the reordering info: ireo(ivtx) = ivtx_old
*
*     version that can also deal with equivalent supervertices
*     topo should be symmetric on entry (not checked any more)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00,
     &     max_sweep = 4

      integer, intent(in) ::
     &     nvtx, nj
      integer, intent(out) ::
     &     ireo(nvtx)
      integer, intent(inout) ::
     &     svtx(nvtx)
      integer(8), intent(inout) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)

      integer(8) ::
     &     vtx0(nvtx)
      integer ::
     &     idx, jdx, sweep, neqv_blocks, svtx0(nvtx), nj_sub, tmp, tmp2,
     &     min_idx, final_sweep, final_block, kdx
      logical ::
     &     changed

      integer, external ::
     &     i8list_cmp, idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'topo_make_unique2')
        write(luout,*) 'topo on entry'
        call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)
      end if

      do idx = 1, nvtx
        ireo(idx) = idx
      end do
      vtx0 = vtx
      call idxsort8(vtx0,ireo,nvtx,+1)
      svtx0 = svtx
      call reoivec(svtx0,ireo,nvtx)
      do idx = 1, nvtx
        svtx(idx) = idxlist(svtx0(idx),svtx0,nvtx,1)
      end do
c      call idxsort(svtx,ireo,nvtx,+1)

      ! resorting svtx in asc. order such that vtx order is not disturbed
      do idx = 2, nvtx
        tmp = svtx(idx)
        tmp2 = ireo(idx)
        jdx = idx-1
        do while (jdx.ge.1)
          if (tmp.ge.svtx(jdx)) exit
          svtx(jdx+1) = svtx(jdx)
          ireo(jdx+1) = ireo(jdx)
          jdx = jdx-1
        end do
        svtx(jdx+1) = tmp
        ireo(jdx+1) = tmp2
      end do

      call reoi8mat(vtx,ireo,nvtx,1,1)
      call reoi8mat(xlines,ireo,nvtx,nj,1)

      ! look for vertices with equal entry on vtx:
      if (nj.gt.0) then
        min_idx = 1
        do idx = 1, nvtx
          if (idx.lt.min_idx) cycle
          jdx = idx
          nj_sub = 1
          vtx_loop: do while(jdx.lt.nvtx.and.(vtx(jdx+1).eq.vtx(idx).or.
     &             svtx(jdx+1).eq.svtx(jdx)))
            if (svtx(jdx+1).eq.svtx(jdx)) then
              nj_sub = nj_sub + 1
              jdx = jdx+1
            else
              ! check if all other vertices match
              do kdx = 2, nj_sub
                if (vtx(jdx+kdx).ne.vtx(idx-1+kdx)) exit vtx_loop
              end do
              jdx = jdx+nj_sub
            end if
          end do vtx_loop
          if (idx+nj_sub-1.ne.jdx) then
            ! sort such that xlines(,nj:1) give ascending seq
            call topo_sort_xlines2(xlines,ireo,idx,jdx,nvtx,nj,nj_sub)
c            min_idx = jdx+1
          end if
          min_idx = jdx+1
        end do
      end if

      ! apply reordering to topo
      call reoi8mat(topo,ireo,nvtx,nvtx,3)
c dbg
c      print *,'before first sweep:'
c      call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)
c dbg

      ! look for vertices with both equal entry on vtx and xlines
      final_sweep = 0
      final_block = 0
      do sweep = 1, max_sweep
        neqv_blocks = 0
        min_idx = 1
        do idx = 1, nvtx
          if (idx.lt.min_idx) cycle
          jdx = idx
          nj_sub = 1
          do while(jdx.lt.nvtx.and.svtx(jdx+1).eq.svtx(jdx))
            nj_sub = nj_sub + 1
            jdx = jdx+1
          end do
          jdx = idx+nj_sub-1
          do while(jdx.lt.nvtx.and.vtx(jdx+1).eq.vtx(idx).and.
     &         i8list_cmp(xlines(jdx+1:jdx+nj_sub,1:nj),
     &                   xlines(jdx-nj_sub+1:jdx,1:nj),nj_sub*nj).eq.0)
            jdx = jdx+nj_sub
          end do
          if (idx+nj_sub-1.ne.jdx) then
            neqv_blocks = neqv_blocks+1
            ! sort such that topo gives ascending seq
            call topo_sort_topo2(topo,changed,ireo,idx,jdx,nvtx,
     &           vtx,svtx,nj_sub)
            min_idx = jdx+1
c dbg
c            print *,'sorting: sweep = ',sweep,' block = ',neqv_blocks,
c     &           ' changed = ',changed
c            call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)
c dbg
            if (changed) then
              final_sweep = sweep + 1
              final_block = neqv_blocks - 1
            else if (neqv_blocks.eq.final_block
     &               .and.sweep.eq.final_sweep) then
              final_block = 0
              exit
            end if
          end if
        end do
        if (final_block.eq.0) exit
c        if (neqv_blocks.le.1.or..not.changed) exit
        if (sweep.eq.max_sweep) then
          write(luout,*) 'max_sweep = ',max_sweep
          write(luout,*) 'neqv_blocks: ',neqv_blocks
          call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)
          call warn('topo_make_unique2',
     &         'sort of topo matrix does not converge')
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'topo on exit'
        call prt_contr_p(luout,svtx,vtx,topo,xlines,nvtx,nj)
      end if

      return
      end

      subroutine topo_sort_xlines2(xlines,ireo,ist,ind,nvtx,nj,nj_sub)

      implicit none

      integer, intent(in) ::
     &     nvtx, nj, ist, ind, nj_sub
      integer, intent(inout) ::
     &     ireo(nvtx)
      integer(8), intent(inout) ::
     &     xlines(nvtx,nj)

      integer ::
     &     ij, ivtx, jvtx, ihlp(nj_sub)
      integer(8) ::
     &     iscr(nj,nj_sub), xl_r(nj,nvtx)

      integer, external ::
     &     i8mat_cmp

      if (nj.eq.0) return

      do ivtx = 1, nvtx
        do ij = 1, nj
          xl_r(ij,ivtx) = xlines(ivtx,nj+1-ij)
        end do
      end do

      do ivtx = ist+nj_sub, ind-nj_sub+1, nj_sub
        iscr(1:nj,1:nj_sub) = xl_r(1:nj,ivtx:ivtx+nj_sub-1)
        ihlp(1:nj_sub) = ireo(ivtx:ivtx+nj_sub-1)
        jvtx = ivtx-nj_sub
        do while (jvtx.ge.ist.and.
     &       i8mat_cmp(iscr(1:nj,1:nj_sub),
     &         xl_r(1:nj,jvtx:jvtx+nj_sub-1),nj,nj_sub).gt.0)
          xl_r(1:nj,jvtx+nj_sub:jvtx+2*nj_sub-1)
     &            = xl_r(1:nj,jvtx:jvtx+nj_sub-1)
          ireo(jvtx+nj_sub:jvtx+2*nj_sub-1)
     &            = ireo(jvtx:jvtx+nj_sub-1)
          jvtx = jvtx-nj_sub
        end do
        xl_r(1:nj,jvtx+nj_sub:jvtx+2*nj_sub-1)
     &            = iscr(1:nj,1:nj_sub)
        ireo(jvtx+nj_sub:jvtx+2*nj_sub-1) = ihlp(1:nj_sub)
      end do

      do ivtx = 1, nvtx
        do ij = 1, nj
          xlines(ivtx,nj+1-ij) = xl_r(ij,ivtx)
        end do
      end do

      return
      end

      subroutine topo_sort_topo2(topo,changed,ireo,ist,ind,nvtx,
     &     vtx,svtx,nj_sub)

      implicit none

      integer, intent(in) ::
     &     nvtx, ist, ind, nj_sub
      integer, intent(inout) ::
     &     ireo(nvtx), svtx(nvtx)
      integer(8), intent(inout) ::
     &     topo(nvtx,nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx)
      logical, intent(out) ::
     &     changed

      integer ::
     &     ivtx, jvtx, ihlp(nj_sub), jhlp(nj_sub), ireo_loc(nvtx)
      integer(8) ::
     &     lscr(nvtx,nj_sub), scr1(nvtx,nj_sub), scr2(nvtx,nj_sub)

      integer, external ::
     &     i8mat_cmp
      
      if (nvtx.eq.0) return

      changed = .false.
      do ivtx = ist+nj_sub, ind-nj_sub+1, nj_sub

        do jvtx = 1, nvtx
          ireo_loc(jvtx) = jvtx
        end do

        lscr(1:nvtx,1:nj_sub) = topo(1:nvtx,ivtx:ivtx+nj_sub-1)
        ihlp(1:nj_sub) = ireo_loc(ivtx:ivtx+nj_sub-1)
        jhlp(1:nj_sub) = ireo(ivtx:ivtx+nj_sub-1)
        jvtx = ivtx-nj_sub
        do while (jvtx.ge.ist)
          scr1 = lscr
          scr2 = topo(1:nvtx,jvtx:jvtx+nj_sub-1)
          scr1(jvtx:jvtx+nj_sub-1,1:nj_sub) = 0
          scr2(ivtx:ivtx+nj_sub-1,1:nj_sub) = 0
          if (i8mat_cmp(scr1,scr2,nvtx,nj_sub).ge.0) exit
          topo(1:nvtx,jvtx+nj_sub:jvtx+2*nj_sub-1)
     &          = topo(1:nvtx,jvtx:jvtx+nj_sub-1)
          ireo_loc(jvtx+nj_sub:jvtx+2*nj_sub-1)
     &          = ireo_loc(jvtx:jvtx+nj_sub-1)
          ireo(jvtx+nj_sub:jvtx+2*nj_sub-1)
     &          = ireo(jvtx:jvtx+nj_sub-1)
          jvtx = jvtx-nj_sub
        end do
        if (jvtx.ne.ivtx-nj_sub) then
          changed = .true.
          topo(1:nvtx,jvtx+nj_sub:jvtx+2*nj_sub-1)
     &          = lscr(1:nvtx,1:nj_sub)
          ireo_loc(jvtx+nj_sub:jvtx+2*nj_sub-1) = ihlp(1:nj_sub)
          ireo(jvtx+nj_sub:jvtx+2*nj_sub-1) = jhlp(1:nj_sub)
c dbg
c          write(*,'(a,12i4)') 'ireo_loc: ',ireo_loc
c dbgend
          ! reorder rows as well
          call reoi8mat(topo,ireo_loc,nvtx,nvtx,1)
        end if
      end do

      return
      end

