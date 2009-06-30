*----------------------------------------------------------------------*
      integer function idx4sg_sb(nel,idx_in,nm,ig,
     &     idorb,idspn,idgam,iy4sg,
     &     iorboff,nelmax,msmax,ngam)
*----------------------------------------------------------------------*
*     get index of string of length nel, with orbitals idorb(1:n),
*     spins idspn(1:n), and orbital IRREPS idgam(1:n) within current
*     subspace
*
*     initial idx, ms, and GAMMA are provided in idx_in, nm, and ig 
*
*     doubly occupied orbitals occur twice, e.g.
*
*      idorb() = (/11,14,14,17/)  idspn() = (/-1,2,2,+1/)
*
*     iy4sg contains arc weights of graph
*     no checks for consistency !
*
*     andreas, june 2006
*
*----------------------------------------------------------------------*
      implicit none

      include "multd2h.h"

      integer, parameter ::
     &     ianum(-1:2) = (/2,0,1,3/)

      integer, intent(in) ::
     &     nel, idorb(nel), idspn(nel), idgam(nel),
     &     nelmax, ngam, idx_in, msmax, iorboff,
     &     iy4sg(1:3,0:nelmax,-msmax:msmax,1:ngam,*)

      integer, intent(inout) ::
     &     nm, ig

      integer ::
     &     idx, nn, i, ispni, ia

      idx = idx_in
      nn = 0
c dbg
c      print *,'top: idx = ',idx
c dbg

      i = 0
      do while(i.lt.nel)
        i = i+1
        ispni = idspn(i)
        ia = ianum(ispni)
        idx = idx + iy4sg(ia,nn,nm,ig,idorb(i)-iorboff)
c dbg
c        print *,'ia,nn,nm,ig,idorb(i),iorboff: ',
c     &           ia,nn,nm,ig,idorb(i),iorboff
c        print *,'idx = ',idx
c dbg
        nn = nn + (iabs(ispni))
        if ((iabs(ispni)).eq.1) nm = nm + ispni
        if ((iabs(ispni)).eq.1) ig = multd2h(ig,idgam(i))
        if (ispni.eq.2) i = i+1
      end do

      idx4sg_sb = idx

      return

      end
