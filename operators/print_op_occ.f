*----------------------------------------------------------------------*
      subroutine print_op_occ(luout,op)
*----------------------------------------------------------------------*
*     print occupations of the operator op
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      integer, intent(in) ::
     &     luout
      type(operator) ::
     &     op

      integer ::
     &     idx, iblk

      do idx = 1, op%n_occ_cls
        iblk = (idx-1)/op%njoined + 1
        call wrt_occ_rstr(luout,iblk,
     &       op%ihpvca_occ(1,1,idx),
     &       op%igasca_restr(1,1,1,1,1,idx),op%ngas,op%nspin)
      end do

      return

      end
