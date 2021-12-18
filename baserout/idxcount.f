*----------------------------------------------------------------------*
      integer function idxcount(inum,ilist,nel,inc)
*----------------------------------------------------------------------*
*     return number of occurence of integer inum in list ilist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inum, nel, inc, ilist(nel)

      integer ::
     &     idx, ipos

      idxcount = 0

      ipos = 1
      do idx = 1, nel
        if (inum.eq.ilist(ipos)) idxcount = idxcount + 1
        ipos = ipos + inc
      end do

      return
      end 
