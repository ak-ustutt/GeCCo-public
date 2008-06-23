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
     &     ntest = 00

      integer, intent(in) ::
     &     nvtx, nj
      integer, intent(out) ::
     &     ireo(nvtx)
      integer(8), intent(inout) ::
     &     vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,nj)

      integer ::
     &     idx, jdx

      integer, external ::
     &     i8list_cmp

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'topo_make_unique')
        write(luout,*) 'topo on entry'
        call prt_contr_p(luout,-1,vtx,topo,xlines,nvtx,nj)
      end if

      do idx = 1, nvtx
        ireo(idx) = idx
      end do
      call idxsort8(vtx,ireo,nvtx,+1)
      call reoi8mat(xlines,ireo,nvtx,nj,1)

      ! look for vertices with equal entry on vtx:
      if (nj.gt.0) then
        do idx = 1, nvtx
          jdx = idx
          do while(jdx.lt.nvtx.and.vtx(jdx+1).eq.vtx(idx))
            jdx = jdx+1
          end do
          if (idx.ne.jdx) then
            ! sort such that xlines(,nj:1) give ascending seq
            call topo_sort_xlines(xlines,ireo,idx,jdx,nvtx,nj)
          end if
        end do
      end if

      ! apply reordering to topo
      call reoi8mat(topo,ireo,nvtx,nvtx,3)

      ! look for vertices with both equal entry on vtx and xlines
      do idx = 1, nvtx
        jdx = idx
        do while(jdx.lt.nvtx.and.vtx(jdx+1).eq.vtx(idx).and.
     &       i8list_cmp(xlines(jdx+1,1:nj),xlines(jdx,1:nj),nj).eq.0)
          jdx = jdx+1
        end do
        if (idx.ne.jdx) then
          ! sort such that topo gives ascending seq
          call topo_sort_topo(topo,ireo,idx,jdx,nvtx,
     &         vtx,xlines,nj)
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

      subroutine topo_sort_topo(topo,ireo,ist,ind,nvtx,
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


      integer ::
     &     ivtx, jvtx, ihlp, jhlp, ireo_loc(nvtx)
      integer(8) ::
     &     iscr(nvtx)

      integer, external ::
     &     i8list_cmp
      
      if (nvtx.eq.0) return

      do ivtx = 1, nvtx
        ireo_loc(ivtx) = ivtx
      end do

      do ivtx = ist+1, ind
        iscr(1:nvtx) = topo(1:nvtx,ivtx)
        ihlp = ireo_loc(ivtx)
        jhlp = ireo(ivtx)
        jvtx = ivtx-1
        do while (jvtx.ge.ist)
          if (i8list_cmp(iscr,topo(1,jvtx),nvtx).le.0) exit
          topo(1:nvtx,jvtx+1) = topo(1:nvtx,jvtx)
          ireo_loc(jvtx+1) = ireo_loc(jvtx)
          ireo(jvtx+1) = ireo(jvtx)
          jvtx = jvtx-1
        end do
        topo(1:nvtx,jvtx+1) = iscr(1:nvtx)
        ireo_loc(jvtx+1) = ihlp
        ireo(jvtx+1) = jhlp
      end do

      ! reorder rows as well
      call reoi8mat(topo,ireo_loc,nvtx,nvtx,1)      

      ! re-check that sequnce is OK
      do ivtx = ist, ind-1
        if (i8list_cmp(topo(1,ivtx),topo(1,ivtx+1),nvtx).lt.0)
     &       call quit(1,'topo_make_unique',
     &       'topo_reo in trouble')
      end do

      return
      end

