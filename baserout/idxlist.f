*----------------------------------------------------------------------*
      integer function idxlist(inum,ilist,nel,inc)
*----------------------------------------------------------------------*
*     return index of first occurence of integer inum in list ilist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     inum, nel, inc, ilist(nel)

      integer ::
     &     idx, ipos

      ipos = 1
      do idx = 1, nel
        if (inum.eq.ilist(ipos)) then
          idxlist = idx
          return
        end if
        ipos = ipos + inc
      end do

      idxlist = -1

      return
      end 
