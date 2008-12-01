*----------------------------------------------------------------------*
      integer function factorial(arg)
*----------------------------------------------------------------------*
*     get factorial of argument arg
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     arg

      integer ::
     &     ii

      factorial = 1

      do ii = 2,arg
        factorial = factorial*ii
      end do

      return
      end
