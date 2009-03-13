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
c      print *,'order = ',iorder
      op%species = spec
c      print *,'species = ',spec
      if (associated(op%ifreq)) deallocate(op%ifreq)
      allocate(op%ifreq(iorder))
      op%ifreq=freq_idx(1:iorder)
c      print *,'ifreq = ',op%ifreq(1:iorder)

      return
      end
