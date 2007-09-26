      subroutine perm_inv(ip_res,ip1,nel)
      implicit none

      integer, intent(in) ::
     &     nel, ip1(nel)
      integer, intent(out) ::
     &     ip_res(nel)
      integer ::
     &     ip_scr(nel),
     &     idx, jdx, kdx

      ip_scr(1:nel) = 0
      do idx = 1, nel
        if (ip1(idx).eq.0) cycle        
        do jdx = 1, nel
          if (ip1(jdx).ne.idx) cycle
          kdx = jdx
          exit
        end do
        ip_scr(idx) = kdx
      end do
      ip_res(1:nel) = ip_scr(1:nel)

      return
      end

      
