*----------------------------------------------------------------------*
      integer function igamstr(nel,idorb,idgam,igamorb)
*----------------------------------------------------------------------*
*     set idgam array (from info on igamorb) and get total IRREP
*     of string (returned as function value)
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     nel, idorb(nel), igamorb(*)
      integer, intent(out) ::
     &     idgam(nel)

      integer ::
     &     iel, itgam

      itgam = 1
      do iel = 1, nel
        idgam(iel) = igamorb(idorb(iel))
        itgam = multd2h(itgam,idgam(iel))
      end do

      igamstr = itgam

      return
      end
