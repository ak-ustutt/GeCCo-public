      subroutine topo_extract_bcres(bcres,
     &     svertex,topo,xlines,vtx_list,
     &     nj,nvtx,nvtx_bcr)

      implicit none

      integer, intent(in) ::
     &     nj, nvtx, nvtx_bcr
      integer(8), intent(in) ::
     &     topo(nvtx,nvtx), xlines(nvtx,nj)
      integer, intent(in) ::
     &     svertex(nvtx), vtx_list(nvtx_bcr)
      integer(8), intent(out) ::
     &     bcres(nvtx_bcr)

      integer ::
     &     idx, jdx, kdx
      integer, external ::
     &     imltlist

      bcres = 0

      do jdx = 1, nvtx_bcr
        idx = vtx_list(jdx)
        do kdx = 1, nj
          bcres(jdx) = bcres(jdx) + xlines(idx,kdx)
        end do
        do kdx = 1, nvtx
          if (imltlist(kdx,vtx_list,nvtx_bcr,1).gt.0) cycle
          bcres(jdx) = bcres(jdx) + topo(idx,kdx)
        end do
      end do

      return
      end
