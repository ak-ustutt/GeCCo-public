*----------------------------------------------------------------------*
      integer function lensubspc(ispc,idspc,n)
*----------------------------------------------------------------------*
*     return number of integers of value<=ispc that occur leftmost
*     in integer string idspc(1:n)
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ispc, n, idspc(n)

      integer ::
     &     idx

      idx = 0
      do while(idx.lt.n)
        if (idspc(idx+1).gt.ispc) exit
        idx = idx+1
      end do
      lensubspc=idx
      
      return
      end
