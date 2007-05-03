*----------------------------------------------------------------------*
      subroutine wrt_contr(lu,contr)
*----------------------------------------------------------------------*
*     write formula on unit lu in condensed form (assuming that 
*     usually word numbers (<32768) are sufficient)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     lu
      type(contraction), intent(in) ::
     &     contr
      
      ! add command record
      write(lu) command_add_contribution, contr%idx_res

      call rw_contr_kernel(-1,lu,contr)

      end
