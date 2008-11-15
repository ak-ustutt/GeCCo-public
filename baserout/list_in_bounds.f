      logical function list_in_bounds(list,len,low,high)
*
*     return TRUE if any element of list is within bounds
*

      implicit none

      integer, intent(in) ::
     &     len, list(len), low, high

      integer ::
     &     iel

      list_in_bounds = .false.
      do iel = 1, len
        list_in_bounds = list_in_bounds.or.
     &       (list(iel).ge.low.and.list(iel).le.high)
      end do

      end
