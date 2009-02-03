*----------------------------------------------------------------------*
      subroutine set_pert_order(op,iorder,spec,freq_idx)
*----------------------------------------------------------------------*
*     set perturbation order, species of operator, and frequency indices
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
      if (associated(op%ifreq)) deallocate(op%ifreq)
      allocate(op%ifreq(iorder))
      op%ifreq=freq_idx(1:iorder)

      return
      end
