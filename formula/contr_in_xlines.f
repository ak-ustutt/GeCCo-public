*----------------------------------------------------------------------*
      subroutine contr_in_xlines(vtxmap,success,
     &     xlines_sub,nvtx_sub,nj_sub,
     &     topo,nvtx)
*----------------------------------------------------------------------*
*     checks whether the contractions on topo between the sub-vertices
*     and the surrounding vertices are contained in the open lines
*     of the sub-vertices and give a valid line structure
*     (no "ring" contractions)
*
*     matthias, dez 2009 (adopted from identify_vertices_i8)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ifc_operators.h'

      logical, intent(out) ::
     &     success
      integer, intent(in) ::
     &     nvtx, nvtx_sub, nj_sub, vtxmap(nvtx)
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx),
     &     xlines_sub(nvtx_sub,nj_sub)

      logical ::
     &     found, down(nj_sub)
      integer ::
     &     ivtx, ivtx_sub, submap(nvtx_sub), ij,
     &     occ_topo(ngastp,2), occ_x(ngastp,2), icnt,
     &     jvtx, connected(nvtx,nvtx), occ_over(ngastp,2), kvtx
      integer(8) ::
     &     base, xlines(nvtx_sub,nj_sub), topo_cur, overlap

      integer, external ::
     &     int8_expand
      integer(8), external ::
     &     int8_pack

      base = pack_base
      xlines = xlines_sub

      ! invert vtxmap
      submap = 0
      do ivtx = 1, nvtx
        if (vtxmap(ivtx).gt.0) then
          if (submap(vtxmap(ivtx)).ne.0) call quit(1,'contr_in_xlines',
     &            'vtxmap has no unique inverse!')
          submap(vtxmap(ivtx)) = ivtx
        end if
      end do
      if (any(submap(1:nvtx_sub).eq.0)) call quit(1,'contr_in_xlines',
     &            'incomplete vtxmap?')

      connected = 0
      success = .true.

      ! find out to which result-sub-vtxs the vtxs are connected
      do ivtx_sub = 1, nvtx_sub
        do ivtx = 1, nvtx
          if (vtxmap(ivtx).ne.0) cycle
          topo_cur = topo(submap(ivtx_sub),ivtx)
          if (topo_cur.eq.0) cycle
          occ_topo = 0
          icnt = int8_expand(topo_cur,base,occ_topo)
          found = .false.
          do ij = 1, nj_sub
            occ_x = 0
            icnt = int8_expand(xlines(ivtx_sub,ij),base,occ_x)
            occ_over = iocc_overlap(occ_x,.false.,occ_topo,.false.)
            overlap = int8_pack(occ_over,ngastp*2,base)
            if (.not.iocc_zero(occ_over)) then
              if (connected(submap(ivtx_sub),ivtx).eq.0.or.
     &            submap(ivtx_sub).lt.ivtx)
     &               connected(submap(ivtx_sub),ivtx) = ij
              occ_topo = occ_topo - occ_over
              xlines(ivtx_sub,ij) = xlines(ivtx_sub,ij) - overlap
              if (iocc_zero(occ_topo)) then
                found = .true.
                exit
              end if
            end if
          end do
          ! test already fails if any contraction is not found in xlines
          if (.not.found) then
            success = .false.
            return
          end if
        end do
      end do

c dbg
c      print *,'contr_in_xlines at work'
c      print *,'submap: ',submap
c      print *,'initial connectivity chart:'
c      call wrtimat2(connected,nvtx,nvtx,nvtx,nvtx)
c dbgend

      ! now spread the connectivity information downwards ...
      do ivtx = 1, nvtx
        do jvtx = ivtx+1, nvtx !upper right triangle
          if (connected(ivtx,jvtx).ne.0) then
            do kvtx = jvtx+1,nvtx
              if (topo(jvtx,kvtx).ne.0.and.
     &            connected(jvtx,kvtx).lt.connected(ivtx,jvtx))
     &             connected(jvtx,kvtx) = connected(ivtx,jvtx)
            end do
          end if
        end do
      end do
      ! and upwards!
      do ivtx = nvtx, 1, -1
        do jvtx = ivtx-1, 1, -1 !lower left triangle
          if (connected(ivtx,jvtx).ne.0) then
            do kvtx = jvtx-1, 1, -1
              if (topo(jvtx,kvtx).ne.0.and.
     &            (connected(jvtx,kvtx).eq.0.or.
     &             connected(jvtx,kvtx).gt.connected(ivtx,jvtx)))
     &             connected(jvtx,kvtx) = connected(ivtx,jvtx)
            end do
          end if
        end do
      end do

c dbg
c      print *,'updated connectivity chart:'
c      call wrtimat2(connected,nvtx,nvtx,nvtx,nvtx)
c dbgend

      ! check that no vtx is connected to the same result-vtx both
      ! down- and upwards!
      do ivtx = 2, nvtx-1
        down = .false.
        do jvtx = 1, ivtx-1
          if (connected(jvtx,ivtx).ne.0)
     &            down(connected(jvtx,ivtx)) = .true.
        end do
        do jvtx = ivtx+1, nvtx
          if (connected(jvtx,ivtx).ne.0.and.down(connected(jvtx,ivtx)))
     &    then
            success = .false.
            return
          end if
        end do
      end do

      return
      end
