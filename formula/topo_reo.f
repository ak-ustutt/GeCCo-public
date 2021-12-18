      subroutine topo_reo(
     &              svertex_new,vtx_new,topo_new,xlines_new,nvtx_new,
     &              svertex,    vtx,    topo,    xlines,    nvtx,
     &              ireo, nj)

      implicit none

      include 'stdunit.h'
      
      integer, intent(in) ::
     &     nvtx, nvtx_new, nj, ireo(nvtx)
      integer(8), intent(out) ::
     &     vtx_new(nvtx_new), topo_new(nvtx_new,nvtx_new),
     &     xlines_new(nvtx_new,nj)
      integer, intent(out) ::
     &     svertex_new(nvtx_new)
      integer(8), intent(in) ::
     &     vtx(nvtx), topo(nvtx,nvtx),
     &     xlines(nvtx,nj)
      integer, intent(in) ::
     &     svertex(nvtx)

      integer ::
     &     idx, jdx

      do idx = 1, nvtx
        if (ireo(idx).lt.1.or.ireo(idx).gt.nvtx_new) then
          write(lulog,*) 'reo: ',ireo(1:nvtx)
          call quit(1,'topo_reo','invalid reordering array')
        end if
      end do

      do idx = 1, nvtx
        svertex_new(ireo(idx)) = svertex(idx)
        vtx_new(ireo(idx)) = vtx(idx)
      end do
      
      topo_new = 0
      do jdx = 1, nvtx
        do idx = 1, nvtx
          topo_new(ireo(idx),ireo(jdx)) =
     &         topo_new(ireo(idx),ireo(jdx)) + topo(idx,jdx)
        end do
      end do

      xlines_new = 0
      do jdx = 1, nj
        do idx = 1, nvtx
          xlines_new(ireo(idx),jdx) =
     &         xlines_new(ireo(idx),jdx) + xlines(idx,jdx)
        end do
      end do

      return
      end
