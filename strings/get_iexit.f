*----------------------------------------------------------------------*
      subroutine get_iexit(iwexit,iw4sg,
     &     nel,maxms,ngam,norb)
*----------------------------------------------------------------------*
*     little wrapper routine to make addressing more legible
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel,maxms,norb,ngam,
     &     iw4sg(0:nel,-maxms:maxms,ngam,norb)

      integer, intent(out) ::
     &     iwexit(-maxms:maxms,ngam,0:nel)

      integer ::
     &     iel

      ! simply copy the last level (and reorder nel --> last index)
      do iel = 0, nel
        iwexit(-maxms:maxms,1:ngam,iel) =
     &     iw4sg(iel,-maxms:maxms,1:ngam,norb)
      end do

      return
      end
