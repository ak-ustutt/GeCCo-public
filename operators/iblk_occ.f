*----------------------------------------------------------------------*
      integer function iblk_occ(iocc,dag,oper)
*----------------------------------------------------------------------*
*     return block of operator oper that corresponds to the
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
     &     iocc_equal_n

      iblk_occ = -1
      dago = oper%dagger
      njoined = oper%njoined

      idx_loop: do idx = 1, oper%n_occ_cls

        iblkoff = (idx-1)*njoined
        if (iocc_equal_n(iocc,dag,
     &                   oper%ihpvca_occ(1,1,iblkoff+1),dago,njoined))
     &                                                              then
          iblk_occ = idx
          exit idx_loop
        end if

      end do idx_loop

      return

      end
