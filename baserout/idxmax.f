*----------------------------------------------------------------------*
      integer function idxmax(ilist,nel,inc)
*----------------------------------------------------------------------*
*     return index of maxinum number in list ilist
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel, inc, ilist(nel)

      integer ::
     &     idx, ipos, imax

      imax = ilist(1)
      idxmax = 1
      ipos = 1 + inc
      do idx = 2, nel
        if (ilist(ipos).gt.imax) then
          idxmax = idx
          imax = ilist(ipos)
        end if
        ipos = ipos + inc
      end do

      return
      end 
