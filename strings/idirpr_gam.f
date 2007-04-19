*----------------------------------------------------------------------*
      integer function idirpr_gam(n,idgam)
*----------------------------------------------------------------------*
*     get total IRREP from direct product of n IRREPs on idgam(1:n)
*     only D2h subgroups, of course ...
*----------------------------------------------------------------------*
      implicit none

      include "multd2h.h"

      integer, intent(in) ::
     &     n, idgam(n)

      integer ::
     &     i

      idirpr_gam = 1
      do i = 1, n
        idirpr_gam = multd2h(idirpr_gam,idgam(i))
      end do

      end 
