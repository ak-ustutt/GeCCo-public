*----------------------------------------------------------------------*
      integer function iblk_occ_inter(iocc,dag,oper)
*----------------------------------------------------------------------*
*     Returns the block of an intermediate-type operator (njoined>1)
*     that corresponds to the hpvx-occupation, iocc.
*     GWR December 2007
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      logical, intent(in) ::
     &     dag
      integer, intent(in) ::
     &     iocc(ngastp,2)
      type(operator), intent(in) ::
     &     oper

      logical ::
     &     dago
      integer ::
     &     idx, ijoin, njoined, iblkoff
      integer ::
     &     iocc_temp(ngastp,2)
      logical, external ::
     &     iocc_equal

      iblk_occ_inter = -1
      dago = oper%dagger
      njoined = oper%njoined

      idx_loop: do idx = 1, oper%n_occ_cls

        iblkoff = (idx-1)*njoined

        iocc_temp(1:ngastp,1:2) =
     &       oper%ihpvca_occ(1:ngastp,1:2,iblkoff+1)
        do ijoin = 2, njoined
          iocc_temp(1:ngastp,1:2) = iocc_temp(1:ngastp,1:2) +
     &         oper%ihpvca_occ(1:ngastp,1:2,iblkoff+ijoin)
        enddo

        if(.not.iocc_equal(iocc,dag,iocc_temp,dago)) cycle idx_loop

        iblk_occ_inter = idx
        exit idx_loop

      enddo idx_loop

      return
      end
