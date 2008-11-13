*----------------------------------------------------------------------*
      subroutine set_pert_order(op,iorder,spec,freq_idx)
*----------------------------------------------------------------------*
*     set perturbation order and species of operator
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'

      type(operator), intent(inout) ::
     &     op

      integer, intent(in) ::
     &     iorder, spec, freq_idx(*)

      op%order = iorder
      op%species = spec
      allocate(op%ifreq(iorder))
      op%ifreq=freq_idx(1:iorder)

      return
      end
