*----------------------------------------------------------------------*
      integer(8) function topo_hash(vtx,topo,xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     get hash function value for arrays defining a contraction
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nvtx,nj
      integer(8), intent(in) ::
     &     vtx(nvtx),topo(nvtx,nvtx),xlines(nvtx,nj)

      integer ::
     &     idx, jdx
      integer(8) ::
     &     ifac

      topo_hash = 0
      ifac = 0
      do idx = 1, nvtx
        ifac = ifac + 1
        topo_hash = topo_hash + vtx(idx)/ifac + mod(vtx(idx),ifac)
      end do
      do jdx = 1, nvtx
        do idx = 1, nvtx
          ifac = ifac + 1
          topo_hash = topo_hash + topo(idx,jdx)/ifac
     &                + mod(topo(idx,jdx),ifac)
        end do
      end do
      do jdx = 1, nj
        do idx = 1, nvtx
          ifac = ifac + 1
          topo_hash = topo_hash + xlines(idx,jdx)/ifac
     &                + mod(xlines(idx,jdx),ifac)
        end do
      end do

      end
