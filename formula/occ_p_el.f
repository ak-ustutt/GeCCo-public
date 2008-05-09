      integer function occ_p_el(occ_p,iel)

      implicit none

      include 'opdim.h'

      integer(8), intent(in) ::
     &     occ_p
      integer, intent(in) ::
     &     iel

      integer(8) ::
     &     scr
      integer ::
     &     i

      scr = occ_p
      do i = 1, iel-1
        scr = scr / pack_base
      end do

      occ_p_el = mod(scr,pack_base)

      return
      end
