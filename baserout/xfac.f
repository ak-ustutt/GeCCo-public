      real(8) function xfac(n)
*
* n !  as double precision real
*
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     n

      integer ::
     &     k

      if( n .lt. 0 ) then
        xfac = 0d0
        write(lulog,*)
     &       ' warning faculty of negative number set to zero '
      else
c
        xfac = 1.0d0
        do k = 2,n
          xfac = xfac * dble(k)
        end do
      end if
c
      return
      end
