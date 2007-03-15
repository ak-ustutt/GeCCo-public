*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
      integer function maxxlvl_op(oper)
*----------------------------------------------------------------------*
*     return maximum excitation level
*     NOTE: works currently only for pure p/h (active spaces: unclear)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      type(operator) ::
     &     oper

      logical ::
     &     dago
      integer ::
     &     idx, ixlvl

      maxxlvl_op = -100000
      dago = oper%dagger

      do idx = 1, oper%n_occ_cls

        ixlvl = (-oper%ihpvca_occ(1,1,idx)
     &           +oper%ihpvca_occ(2,1,idx)
     &           +oper%ihpvca_occ(1,2,idx)
     &           -oper%ihpvca_occ(2,2,idx))/2
        if (dago) ixlvl = -ixlvl
        maxxlvl_op = max(maxxlvl_op,ixlvl)

      end do

      return

      end
