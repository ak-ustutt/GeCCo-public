      subroutine set_op_ldim_c(ldimc,ldima,hpvxc,hpvxa,nstr,nc,na)
      implicit none
      
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     ldimc(nc), ldima(na)
      integer, intent(in) ::
     &     hpvxc(nc), hpvxa(na),
     &     nc, na, nstr(nc+na)

      integer ::
     &     idx, ldim, ioff, idx_hpvx, hpvx

      ldim = 1
      do idx_hpvx = 1, ngastp
        hpvx = hpvxseq(idx_hpvx)

        do idx = 1, nc
          if (hpvxc(idx).ne.hpvx) cycle
          ldimc(idx) = ldim
          ldim = ldim * nstr(idx)
        end do

        ioff = nc
        do idx = 1, na
          if (hpvxa(idx).ne.hpvx) cycle
          ldima(idx) = ldim
          ldim = ldim * nstr(ioff+idx)
        end do
      end do

      return
      end
      
