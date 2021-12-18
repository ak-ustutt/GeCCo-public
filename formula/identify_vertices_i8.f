*----------------------------------------------------------------------*
      subroutine identify_vertices_i8(vtxmap,success,nskip,
     &     svtx_sub,vtx_sub,topo_sub,nvtx_sub,
     &     svtx,vtx,topo,nvtx)
*----------------------------------------------------------------------*
*     initial (and astonishingly simple) routine to identify a 
*     sub-contraction in a contraction
*     topomaps are fantastic ...
*
*     rewritten for arbitrarily ordered contractions,
*     and with the possibility to skip nskip matching possibilities
*     (this way, a loop over nskip wrapped around this routine
*      may be used to systematically get all successful possibilities)
*----------------------------------------------------------------------*

      implicit none

      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     nvtx, nvtx_sub, svtx(nvtx), svtx_sub(nvtx_sub), nskip
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     vtx_sub(nvtx_sub), topo_sub(nvtx_sub,nvtx_sub)
      integer, intent(out) ::
     &     vtxmap(nvtx)

      integer ::
     &     ivtx, jvtx, ivtx_sub, jvtx_sub,
     &     lenmap_sub(nvtx_sub), vtxmap_sub(nvtx_sub,nvtx),
     &     idist(nvtx_sub), one(nvtx_sub), iskip
      logical, external ::
     &     next_dist2

      success = .false.

      vtxmap = 0
      one = 1
      lenmap_sub = 0
      vtxmap_sub = 0
      iskip = 0

      ! set up vtxmap_sub: to which vtxs could a sub-vtx be assigned?
      do ivtx_sub = 1, nvtx_sub
        do ivtx = 1, nvtx
          if (vtx(ivtx).eq.vtx_sub(ivtx_sub)) then
            lenmap_sub(ivtx_sub) = lenmap_sub(ivtx_sub)+1
            vtxmap_sub(ivtx_sub,lenmap_sub(ivtx_sub)) = ivtx
          end if
        end do
      end do

      ! for all sub-vtxs there must be at least one possibility
      if (any(lenmap_sub(1:nvtx_sub).eq.0)) return

      idist = 1
      idist(1) = 0

      dist_loop: do while(next_dist2(idist,nvtx_sub,one,lenmap_sub,1))
        vtxmap = 0
        do ivtx_sub = 1, nvtx_sub
          if (vtxmap(vtxmap_sub(ivtx_sub,idist(ivtx_sub))).ne.0)
     &          cycle dist_loop
          vtxmap(vtxmap_sub(ivtx_sub,idist(ivtx_sub))) = ivtx_sub
        end do

        do ivtx = 1, nvtx-1
          ivtx_sub = vtxmap(ivtx)
          if (ivtx_sub.eq.0) cycle
          do jvtx = ivtx, nvtx
            jvtx_sub = vtxmap(jvtx)
            if (jvtx_sub.eq.0) cycle

              ! either same supervtx or different supervtx
              if (svtx(ivtx).eq.svtx(jvtx).neqv.
     &            svtx_sub(ivtx_sub).eq.svtx_sub(jvtx_sub))
     &          cycle dist_loop

              ! same connectivity
              if (topo(ivtx,jvtx).ne.topo_sub(ivtx_sub,jvtx_sub))
     &          cycle dist_loop
          end do
        end do

        ! skip if requested
        if (iskip.lt.nskip) then
          iskip = iskip + 1
          cycle dist_loop
        end if

        ! this is the vtxmap we were looking for
        success = .true.
        exit dist_loop
      end do dist_loop

      return
      end
