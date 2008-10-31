*----------------------------------------------------------------------*
      subroutine topo_make_unique(ireo,vtx,topo,xlines,nvtx,nj)
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
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00,
     &     max_sweep = 3

      integer, intent(in) ::
     &     nvtx, nj
      integer, intent(out) ::
     &     ireo(nvtx)
      integer(8), intent(inout) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)

      integer ::
     &     idx, jdx, sweep, neqv_blocks
      logical ::
     &     changed

      integer, external ::
     &     i8list_cmp

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'topo_make_unique')
        write(luout,*) 'topo on entry'
        call prt_contr_p(luout,-1,vtx,topo,xlines,nvtx,nj)
      end if
c test - symmetrize topo
      do idx = 1, nvtx
        do jdx = 1, idx-1
          topo(jdx,idx) = topo(idx,jdx)
        end do
      end do
c test

      do idx = 1, nvtx
        ireo(idx) = idx
      end do
      call idxsort8(vtx,ireo,nvtx,+1)
      call reoi8mat(xlines,ireo,nvtx,nj,1)

      ! look for vertices with equal entry on vtx:
      if (nj.gt.0) then
        do idx = 1, nvtx
          jdx = idx
c dbg fix by mh
          if (jdx.lt.nvtx) then
c dbg original
          do while(jdx.lt.nvtx.and.vtx(jdx+1).eq.vtx(idx))
            jdx = jdx+1
          end do
c dbg resume fix
          end if
c dbg end fix
          if (idx.ne.jdx) then
            ! sort such that xlines(,nj:1) give ascending seq
            call topo_sort_xlines(xlines,ireo,idx,jdx,nvtx,nj)
          end if
        end do
      end if

      ! apply reordering to topo
      call reoi8mat(topo,ireo,nvtx,nvtx,3)

      ! look for vertices with both equal entry on vtx and xlines
      do sweep = 1, max_sweep
        neqv_blocks = 0
        do idx = 1, nvtx
          jdx = idx
c dbg fix by mh
          if (jdx.lt.nvtx) then
c dbg original
          do while(jdx.lt.nvtx.and.vtx(jdx+1).eq.vtx(idx).and.
     &         i8list_cmp(xlines(jdx+1,1:nj),xlines(jdx,1:nj),nj).eq.0)
            jdx = jdx+1
          end do
c dbg resume fix
          end if
c dbg end fix
          if (idx.ne.jdx) then
            neqv_blocks = neqv_blocks+1
            ! sort such that topo gives ascending seq
            call topo_sort_topo(topo,changed,ireo,idx,jdx,nvtx,
     &           vtx,xlines,nj)
c dbg
c            print *,'sorting: sweep = ',sweep,' block = ',neqv_blocks,
c     &           ' changed = ',changed
c            call prt_contr_p(luout,-1,vtx,topo,xlines,nvtx,nj)
c dbg
          end if
        end do
        if (neqv_blocks.le.1.or..not.changed) exit
        if (sweep.eq.max_sweep) then
          write(luout,*) 'max_sweep = ',max_sweep
          write(luout,*) 'neqv_blocks: ',neqv_blocks
          call prt_contr_p(luout,-1,vtx,topo,xlines,nvtx,nj)
          call quit(1,'topo_make_unique',
     &         'sort of topo matrix does not converge')
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'topo on exit'
        call prt_contr_p(luout,-1,vtx,topo,xlines,nvtx,nj)
      end if

      return
      end

      subroutine topo_sort_xlines(xlines,ireo,ist,ind,nvtx,nj)

      implicit none

      integer, intent(in) ::
     &     nvtx, nj, ist, ind
      integer, intent(inout) ::
     &     ireo(nvtx)
      integer(8), intent(inout) ::
     &     xlines(nvtx,nj)

      integer ::
     &     ij, ivtx, jvtx, ihlp
      integer(8) ::
     &     iscr(nj), xl_r(nj,nvtx)

      integer, external ::
     &     i8list_cmp
      
      if (nj.eq.0) return

      do ivtx = 1, nvtx
        do ij = 1, nj
          xl_r(ij,ivtx) = xlines(ivtx,nj+1-ij)
        end do
      end do

      do ivtx = ist+1, ind
        iscr(1:nj) = xl_r(1:nj,ivtx)
        ihlp = ireo(ivtx)
        jvtx = ivtx-1
        do while (jvtx.ge.ist.and.
     &       i8list_cmp(iscr(1:nj),xl_r(1:nj,jvtx),nj).gt.0)
          xl_r(1:nj,jvtx+1) = xl_r(1:nj,jvtx)
          ireo(jvtx+1) = ireo(jvtx)
          jvtx = jvtx-1
        end do
        xl_r(1:nj,jvtx+1) = iscr(1:nj)
        ireo(jvtx+1) = ihlp
      end do

      do ivtx = 1, nvtx
        do ij = 1, nj
          xlines(ivtx,nj+1-ij) = xl_r(ij,ivtx)
        end do
      end do

      return
      end

      subroutine topo_sort_topo(topo,changed,ireo,ist,ind,nvtx,
     &     vtx,xlines,nj)

      implicit none

      integer, intent(in) ::
     &     nvtx, ist, ind, nj
      integer, intent(inout) ::
     &     ireo(nvtx)
      integer(8), intent(inout) ::
     &     topo(nvtx,nvtx)
      integer(8), intent(in) ::
     &     vtx(nvtx), xlines(nvtx,nj)
      logical, intent(out) ::
     &     changed

      integer ::
     &     ivtx, jvtx, ihlp, jhlp, ireo_loc(nvtx)
      integer(8) ::
     &     lscr(nvtx), iscr(nvtx), jscr(nvtx)

      integer, external ::
     &     i8list_cmp
      
      if (nvtx.eq.0) return

      do ivtx = 1, nvtx
        ireo_loc(ivtx) = ivtx
      end do

      changed = .false.
      do ivtx = ist+1, ind
        lscr(1:nvtx) = topo(1:nvtx,ivtx)
c        iscr(1:ivtx) = topo(1:ivtx,ivtx)
c        iscr(ivtx+1:nvtx) = topo(ivtx,ivtx+1:nvtx)
        ihlp = ireo_loc(ivtx)
        jhlp = ireo(ivtx)
        jvtx = ivtx-1
        do while (jvtx.ge.ist)
          if (i8list_cmp(lscr,topo(1,jvtx),ind).le.0) exit
c          jscr(1:jvtx) = topo(1:jvtx,jvtx)
c          jscr(jvtx+1:nvtx) = topo(jvtx,jvtx+1:nvtx)
c          if (i8list_cmp(iscr,topo(1,jvtx),nvtx).le.0) exit
c          if (i8list_cmp(iscr,jscr,nvtx).le.0) exit
          topo(1:nvtx,jvtx+1) = topo(1:nvtx,jvtx)
          ireo_loc(jvtx+1) = ireo_loc(jvtx)
          ireo(jvtx+1) = ireo(jvtx)
          jvtx = jvtx-1
        end do
        changed = changed.or.jvtx.ne.ivtx-1
        topo(1:nvtx,jvtx+1) = lscr(1:nvtx)
        ireo_loc(jvtx+1) = ihlp
        ireo(jvtx+1) = jhlp
      end do

      ! reorder rows as well
      call reoi8mat(topo,ireo_loc,nvtx,nvtx,1)      

c      ! re-check that sequence is OK
c      do ivtx = ist, ind-1
c        if (i8list_cmp(topo(1,ivtx),topo(1,ivtx+1),nvtx).lt.0)
cc        if (i8list_cmp(topo(1,ivtx),topo(1,ivtx+1),ind).lt.0)
c     &       call quit(1,'topo_make_unique',
c     &       'topo_reo in trouble')
c      end do

      return
      end

