*----------------------------------------------------------------------*
      subroutine identify_vertices_i8(vtxmap,success,
     &     vtx_sub,topo_sub,nvtx_sub,
     &     vtx,topo,nvtx)
*----------------------------------------------------------------------*
*     initial (and astonishingly simple) routine to identify a 
*     sub-contraction in a contraction
*     topomaps are fantastic ...
*c     we are assuming canonically ordered contractions, as we expect
*c     that the vertices in both the contraction and sub-contraction
*c     should come in the same order (so we need not consider
*c     permutations)
*     extended to arbitrary order: any permutation is allowed
*----------------------------------------------------------------------*

      implicit none

      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     nvtx, nvtx_sub
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     vtx_sub(nvtx_sub), topo_sub(nvtx_sub,nvtx_sub)
      integer, intent(out) ::
     &     vtxmap(nvtx)

      logical ::
     &     ok, matched(nvtx), found1, found2
      integer ::
     &     ivtx, jvtx, ivtx_sub, jvtx_sub, ivtx2

      success = .false.

      vtxmap = 0
      ! compare vtx-type
      isub_loop: do ivtx_sub = 1, nvtx_sub
        found1 = .false.
        i_loop: do ivtx = 1, nvtx
          if (vtxmap(ivtx).gt.0) cycle i_loop
          if (vtx(ivtx).ne.vtx_sub(ivtx_sub)) cycle i_loop
          ! compare topo of vtx
          matched = .false.
          do jvtx_sub = 1, nvtx_sub
            found2 = .false.
            do jvtx = 1, nvtx
              if (matched(jvtx)) cycle
              ! must connect to same vertex ...
              if (vtx(jvtx).ne.vtx_sub(jvtx_sub)) cycle
              ! ... with same connection
              if (topo(jvtx,ivtx).eq.topo_sub(jvtx_sub,ivtx_sub)) then
                matched(jvtx) = .true.
                found2 = .true.
                exit ! increase jvtx_sub
              end if
            end do
            if (.not.found2) cycle i_loop ! try next vertex: increase ivtx
          end do
          ! all topos found: set the entry
          vtxmap(ivtx) = ivtx_sub 
          found1 = .true.
          success = ivtx_sub.eq.nvtx_sub
          exit i_loop ! increase ivtx_sub
        end do i_loop 
        if (.not.found1) exit isub_loop ! no identical vertex
      end do isub_loop

      return
      end
