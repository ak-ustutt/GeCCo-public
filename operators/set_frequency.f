*----------------------------------------------------------------------*
      subroutine set_frequency(mel,order,ifreq)
*----------------------------------------------------------------------*
*     set perturbation order and species of operator
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'
      include 'ifc_input.h'
      include 'def_me_list.h'
      include 'def_filinf.h'

      type(me_list), intent(inout) ::
     &     mel

      integer, intent(in) ::
     &     order, ifreq

      real(8) ::
     &     freq(order)

      freq = 0d0
      call get_argument_value('calculate.experimental','freq',
     &                        xarr=freq(1:order-1))
      freq(order) = -sum(freq)
      mel%frequency = freq(ifreq)

      return
      end
