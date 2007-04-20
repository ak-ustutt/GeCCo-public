*----------------------------------------------------------------------*
      subroutine init_op_files(ffops,ops,nops)
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_operator.h'
      include 'def_filinf.h'

      integer, intent(in) ::
     &     nops

      type(operator), intent(in) ::
     &     ops(nops)
      type(filinf), intent(in) ::
     &     ffops(nops)

      integer ::
     &     iop

      do iop = 1, nops
        call file_init(ffops(iop),
     &       'op_'//trim(ops(iop)%name)//'_elements.da',ftyp_da_unf,
     &       lblk_da)
      end do

      return
      end
