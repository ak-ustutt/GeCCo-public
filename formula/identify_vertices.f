*----------------------------------------------------------------------*
      subroutine identify_vertices(vtxmap,success,
     &     vtx_sub,topo_sub,nvtx_sub,
     &     vtx,topo,nvtx)
*----------------------------------------------------------------------*
*     initial (and astonishingly simple) routine to identify a 
*     sub-contraction in a contraction
*     topomaps are fantastic ...
*     we are assuming canonically ordered contractions, as we expect
*     that the vertices in both the contraction and sub-contraction
*     should come in the same order (so we need not consider
*     permutations)
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
     &     ivtx, jvtx, ivtx_sub, jvtx_sub


      success = .false.

      vtxmap = 0
      ivtx_sub = 1
      ! compare vtx-type
      do ivtx = 1, nvtx
        if (vtx(ivtx).ne.vtx_sub(ivtx_sub)) cycle
        ! compare topo of vtx
        jvtx_sub = 1
        do jvtx = 1, nvtx
          ! must connect to same vertex ...
          if (vtx(jvtx).ne.vtx_sub(jvtx_sub)) cycle
          ! ... with same connection
          if (topo(jvtx,ivtx).ne.topo_sub(jvtx_sub,ivtx_sub)) cycle
          if (jvtx_sub.eq.nvtx_sub) exit
          jvtx_sub = jvtx_sub+1
        end do
        ! all topos found:
        if (jvtx_sub.eq.nvtx_sub) then
          ! set the entry
          vtxmap(ivtx) = ivtx_sub 
          success = ivtx_sub.eq.nvtx_sub
          if (success) exit
          ivtx_sub = ivtx_sub+1
        end if
      end do      

      return
      end
