*----------------------------------------------------------------------*
      subroutine identify_vertices3(vtxmap,success,
     &     vtx_sub,topo_sub,nvtx_sub,
     &     vtx,topo,nvtx)
*----------------------------------------------------------------------*
*     initial (and astonishingly simple) routine to identify a 
*     sub-contraction in a contraction
*     topomaps are fantastic ...
*     we are not assuming canonically ordered contractions
*     NB We assume that the first intm vertex is unique in both
*     contractions.
*----------------------------------------------------------------------*

      implicit none

      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     nvtx, nvtx_sub, 
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     vtx_sub(nvtx_sub), topo_sub(nvtx_sub,nvtx_sub)
      integer, intent(out) ::
     &     vtxmap(nvtx)

      logical ::
     &     ok
      integer ::
     &     ivtx, jvtx, kvtx, idx


      success = .false.

      vtxmap = 0
      idx = 0

      do jvtx = 1, nvtx_sub
        contr_loop: do ivtx = 1, nvtx
          if (vtx(ivtx).ne.vtx_sub(jvtx))cycle
          if(jvtx.eq.1)then
            vtxmap(ivtx) = jvtx
            idx = idx + 1
            exit contr_loop
          endif
          if(vtxmap(ivtx).eq.0)then
            ok = .true.
            do kvtx = 1, nvtx
              if(vtxmap(kvtx).ne.0)then
                ok = ok.and.
     &               topo(ivtx,kvtx).eq.topo_sub(jvtx,vtxmap(kvtx))
              endif
            enddo
            if(ok)then
              vtxmap(ivtx) = jvtx
              idx = idx + 1
              success = idx.eq.nvtx_sub
              if(success) return
              exit contr_loop
            endif
          endif
        enddo contr_loop
        if(jvtx.eq.1.and.idx.eq.0) return
      enddo

      return
      end
