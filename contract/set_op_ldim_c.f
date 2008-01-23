*----------------------------------------------------------------------*
      subroutine set_op_ldim_c(ldimc,ldima,
     &     hpvxc,hpvxa,nstr,nc,na,transpose)
*----------------------------------------------------------------------*
*     set the sequence of dimensions for addressing the operator list
*----------------------------------------------------------------------*

      implicit none
      
      include 'opdim.h'
      include 'hpvxseq.h' ! <- determines the sequence of HPVX

      logical, intent(in) ::
     &     transpose
      integer, intent(out) ::
     &     ldimc(nc), ldima(na)
      integer, intent(in) ::
     &     hpvxc(nc), hpvxa(na),
     &     nc, na, nstr(nc+na)

      integer ::
     &     idx, ldim, ioff, idx_hpvx, hpvx

      ! C A sequence
      ! NB: for pure (de)excitation operators, transpose has no effect
      !     as hpvxseq enforces P major storage
      if (.not.transpose) then

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
        
      ! A C sequence
      else

        ldim = 1
        do idx_hpvx = 1, ngastp
          hpvx = hpvxseq(idx_hpvx)        
          
          ioff = nc
          do idx = 1, na
            if (hpvxa(idx).ne.hpvx) cycle
            ldima(idx) = ldim
            ldim = ldim * nstr(ioff+idx)
          end do

          do idx = 1, nc
            if (hpvxc(idx).ne.hpvx) cycle
            ldimc(idx) = ldim
            ldim = ldim * nstr(idx)
          end do          
        end do
        
      end if

      return
      end
      
