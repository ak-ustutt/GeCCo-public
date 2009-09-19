      subroutine reo_full2sym(xblk,xfull,nfull,lblk,nblk,gamma,psign)
*
*     reorder triangular full matrix to symmetry blocked
*     (column symmetry ordered)
*
*     adapted from Jeppe
*      

      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     nblk, gamma, lblk(nblk), nfull
      real(8), intent(in) ::
     &     xfull(nfull*(nfull+1)/2), psign
      real(8), intent(out) ::
     &     xblk(*)

      integer ::
     &     ioff_blk, ioff_row, ioff_col,
     &     irsym, icsym, nrow, ncol,
     &     icdx, irdx, icdx_tot, irdx_tot, icr_min, icr_max
      real(8) ::
     &     fac

      if (abs(abs(psign)-1d0).gt.1d-12)
     &     call quit(1,'reo_full2sym','phase sign must be +1 or -1')

      ! loop over symmetry blocks
      ioff_blk = 0
      do icsym = 1, nblk
        irsym = multd2h(icsym,gamma)
        nrow = lblk(irsym)
        ncol = lblk(icsym)

        ! offsets
        ioff_row = 0
        if (irsym.gt.1) ioff_row = sum(lblk(1:irsym-1))
        ioff_col = 0
        if (icsym.gt.1) ioff_col = sum(lblk(1:icsym-1))

        do icdx = 1, ncol
          icdx_tot = icdx+ioff_col
          do irdx = 1, nrow
            irdx_tot = irdx+ioff_row
            icr_max = max(icdx_tot,irdx_tot)
            icr_min = min(icdx_tot,irdx_tot)

            fac = 1d0
            if (icdx_tot.gt.irdx_tot) fac = psign ! sign consistent to DALTON
                                      ! see DALTON/pdpack/linextra.F line 109

            xblk(ioff_blk + (icdx-1)*nrow+irdx) =
     &           fac*xfull(icr_max*(icr_max-1)/2+icr_min)
          end do
        end do

        ioff_blk = ioff_blk + nrow*ncol

      end do

      return
      end
