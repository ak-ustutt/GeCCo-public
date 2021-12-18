*----------------------------------------------------------------------*
      subroutine set_nca_vtx(nca_vtx,occ_vtx,nvtx)
*----------------------------------------------------------------------*
*     just report the number of C/A on each vertex of occ_vtx on
*     array occ_vtx
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     nvtx, occ_vtx(ngastp,2,nvtx)

      integer, intent(out) ::
     &     nca_vtx(nvtx)

      integer ::
     &     ivtx, ica, igastp

      do ivtx = 1, nvtx
        nca_vtx(ivtx) = 0
        do ica = 1, 2
          do igastp = 1, ngastp
            nca_vtx(ivtx) = nca_vtx(ivtx) + occ_vtx(igastp,ica,ivtx)
          end do
        end do
      end do

      return
      end
