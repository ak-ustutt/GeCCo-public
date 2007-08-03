*----------------------------------------------------------------------*
      subroutine ms2idxms(idxms,ms,occ,n)
*----------------------------------------------------------------------*
*     convert Ms-value (actually 2 times Ms) to index value
*     e.g.  2, 0, -2  to 1, 2, 3 ....
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     n, ms(n), occ(n)
      integer, intent(out) ::
     &     idxms(n)

      integer ::
     &     i

      do i = 1, n
        idxms(i) = (occ(i) - ms(i))/2 + 1
      end do

      return
      end
