      integer function ibico(n,k)
*     binomial coefficient n over k

      implicit none

      integer, intent(in) ::
     &     n, k
      integer ::
     &     l, nmkp1
      integer, external ::
     &     ifac

      if (k.lt.0.or.n.lt.k)
     &   call quit(1,'ibico','k<0 or n<k')

      if (k.gt.n/2) then
        nmkp1 = k+1 !use (n over k) = (n over (n-k))
      else
        nmkp1 = n-k+1
      end if

      ibico = 1
      do l = nmkp1, n
        ibico = ibico*l
      end do
      ibico = ibico/ifac(n-nmkp1+1)

      return
      end
