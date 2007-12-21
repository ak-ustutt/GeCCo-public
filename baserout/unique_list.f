*----------------------------------------------------------------------*
      subroutine unique_list(list,n)
*----------------------------------------------------------------------*
*     input: list of n(input) integer values
*     output: sorted list (incr. values) of n(output) unique integers
*----------------------------------------------------------------------*

      implicit none

      integer, intent(inout) ::
     &     n, list(n)

      integer ::
     &     idx, jdx

      ! sort to ascending order
      call isort(list,n,+1)

      jdx = 1
      do idx = 2, n
        if (list(idx).eq.list(jdx)) cycle
        jdx = jdx+1        
        list(jdx) = list(idx)
      end do

      n = jdx

      return
      end
