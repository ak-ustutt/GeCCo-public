*----------------------------------------------------------------------*
      subroutine print_op_occ(lulog,op)
*----------------------------------------------------------------------*
*     print occupations of the operator op
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'

      integer, intent(in) ::
     &     lulog
      type(operator) ::
     &     op

      integer ::
     &     idx, iblk

      do idx = 1, op%n_occ_cls*op%njoined
        iblk = (idx-1)/op%njoined + 1
        if (op%formal_blk(iblk).and.mod(idx-1,op%njoined).eq.0)
     &       write(lulog,*) '--- formal: ---'
        call wrt_occ_rstr(lulog,iblk,
     &       op%ihpvca_occ(1,1,idx),
     &       op%igasca_restr(1,1,1,1,1,idx),op%ngas,op%nspin)
      end do

      return

      end
