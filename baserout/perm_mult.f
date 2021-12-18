      subroutine perm_mult(ip_res,ip1,ip2,nel)
      implicit none

      integer, intent(in) ::
     &     nel, ip1(nel), ip2(nel)
      integer, intent(out) ::
     &     ip_res(nel)
      integer ::
     &     ip_scr(nel)
      integer ::
     &     idx

      ip_scr(1:nel) = 0
      do idx = 1, nel
        if (ip2(idx).eq.0) cycle
        ip_scr(idx) = ip1(ip2(idx))
      end do

      ip_res = ip_scr

      return
      end

      
