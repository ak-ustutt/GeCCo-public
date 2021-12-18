*----------------------------------------------------------------------*
      subroutine set_hermitian(op,hermitian)
*----------------------------------------------------------------------*
*     set "hermitian" attribute
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'

      type(operator), intent(inout) ::
     &     op

      integer, intent(in) ::
     &     hermitian

      op%hermitian = hermitian

      return
      end
