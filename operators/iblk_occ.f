*----------------------------------------------------------------------*
      integer function iblk_occ(iocc,dag,oper)
*----------------------------------------------------------------------*
*     return block of operator oper that correspondents to the
*     hpv-occupation iocc
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      logical, intent(in) ::
     &     dag
      integer, intent(in) ::
     &     iocc(ngastp,2,*)
      type(operator) ::
     &     oper

      logical ::
     &     dago
      integer ::
     &     idx, ijoin, njoined, iblkoff
      integer, external ::
     &     iocc_equal

      iblk_occ = -1
      dago = oper%dagger
      njoined = oper%njoined

      idx_loop: do idx = 1, oper%n_occ_cls

        iblkoff = (idx-1)*njoined
        if (iocc_equal(iocc,dag,oper%ihpvca_occ(1,1,iblkoff+1),dago))
     &                                                              then

          do ijoin = 2, njoined
            if (.not.iocc_equal(iocc(1,1,ijoin),dag,
     &                         oper%ihpvca_occ(1,1,iblkoff+ijoin),dago))
     &           cycle idx_loop
          end do
          iblk_occ = idx
          exit idx_loop
        end if

      end do idx_loop

      return

      end
