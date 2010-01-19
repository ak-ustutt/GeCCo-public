      subroutine pack_occ(occ_pack,occ_full,nvtx)

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     nvtx
      integer(8), intent(out) ::
     &     occ_pack(nvtx)
      integer, intent(in) ::
     &     occ_full(ngastp*2,nvtx)

      integer ::
     &     ivtx

      integer, external ::
     &     int8_pack

      do ivtx = 1, nvtx
        occ_pack(ivtx) = int8_pack(occ_full(1:,ivtx),ngastp*2,pack_base)
      end do

      end
