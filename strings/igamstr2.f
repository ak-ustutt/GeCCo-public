*----------------------------------------------------------------------*
      integer function igamstr2(nel,idorb,igamorb)
*----------------------------------------------------------------------*
*     as igamstr, but the IRREP distribution is not set
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     nel, idorb(nel), igamorb(*)

      integer ::
     &     iel, itgam

      itgam = 1
      do iel = 1, nel
c dbg
c        print *,'.  ',iel
c        print *,'.. ',idorb(iel)
c        print *,'...',igamorb(idorb(iel))
c dbg
        itgam = multd2h(itgam,igamorb(idorb(iel)))
      end do

      igamstr2 = itgam

      return
      end
