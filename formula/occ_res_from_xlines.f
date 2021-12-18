*----------------------------------------------------------------------*
      subroutine occ_res_from_xlines(occ_res,
     &     xlines,nvtx,nj)
*----------------------------------------------------------------------*
*     generate result of the full contraction
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nj, nvtx
      integer(8), intent(in) ::
     &     xlines(nvtx,nj)
      integer(8), intent(out) ::
     &     occ_res(nj)

      integer ::
     &     idx, jdx, kdx, jdxnd
      integer, external ::
     &     imltlist

      occ_res = 0

      do kdx = 1, nj
        do idx = 1, nvtx
          occ_res(kdx) = occ_res(kdx) + xlines(idx,kdx)
        end do
      end do

      return
      end
