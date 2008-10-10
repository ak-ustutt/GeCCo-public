      subroutine unpack_occ(occ_full,occ_pack,nvtx)

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     nvtx
      integer(8), intent(in) ::
     &     occ_pack(nvtx)
      integer, intent(out) ::
     &     occ_full(ngastp*2,nvtx)

      integer(8) ::
     &     scr1, scr2
      integer ::
     &     ivtx, iel

      do ivtx = 1, nvtx
        scr1 = occ_pack(ivtx)
        do iel = 1, ngastp*2
          scr2 = mod(scr1,pack_base)
          occ_full(iel,ivtx) = scr2
          scr1 = scr1 / pack_base
        end do        
      end do

      end
