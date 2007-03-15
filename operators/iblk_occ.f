
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
     &     iocc(ngastp,2)
      type(operator) ::
     &     oper

      logical ::
     &     dago
      integer ::
     &     idx
      integer, external ::
     &     iocc_equal

      iblk_occ = -1
      dago = oper%dagger

      do idx = 1, oper%n_occ_cls

        if (iocc_equal(iocc,dag,oper%ihpvca_occ(1,1,idx),dago)) then
          iblk_occ = idx
          exit
        end if

      end do

      return

      end
