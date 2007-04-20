      integer function ifac(n)
c
c n !
c
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     n
      integer ::
     &     k

      if( n .lt. 0 ) then
        ifac = 0
        write(luout,*)
     &       ' warning: faculty of negative number set to zero'
      else
c
        ifac = 1
        do k = 2,n
          ifac = ifac * k
        end do
      end if
c
      return
      end
