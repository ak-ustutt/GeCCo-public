      integer(8) function occ_overlap_p(occ1,occ2)

      implicit none

      include 'opdim.h'

      integer(8), intent(in) ::
     &     occ1, occ2

      integer ::
     &     occ_scr(ngastp*2)
      integer(8) ::
     &     scr1, scr2
      integer ::
     &     ivtx, iel

      integer(8), external ::
     &     int8_pack

      scr1 = occ1
      scr2 = occ2
      do iel = 1, ngastp*2
        occ_scr(iel) = min(mod(scr1,pack_base),
     &                     mod(scr2,pack_base))
        scr1 = scr1 / pack_base
        scr2 = scr2 / pack_base
      end do        

      occ_overlap_p = int8_pack(occ_scr,ngastp*2,pack_base)

      end
