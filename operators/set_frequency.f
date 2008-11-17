*----------------------------------------------------------------------*
      subroutine set_frequency(mel,order)
*----------------------------------------------------------------------*
*     set frequency assigned to ME-list
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
     &     order

      real(8) ::
     &     freq(order)

      integer ::
     &     ii

      if (order.gt.0) then

        ! get complete user defined frequency array, sum of frequencies is zero
        freq = 0d0
        call get_argument_value('calculate.experimental','freq',
     &                          xarr=freq(1:order-1))
        freq(order) = -sum(freq)

        ! frequency is sum of frequencies associated with frequency indices
        if (.not.associated(mel%op%ifreq)) call quit(1,'set_frequency',
     &       'no frequency index associated to operator')
        mel%frequency = 0d0
        do ii = 1,mel%op%order
          mel%frequency = mel%frequency + freq(mel%op%ifreq(ii))
        end do
        ! flip sign if operator species = 1 (T-amplitudes)
        if (mel%op%species.eq.1) then
          mel%frequency = -mel%frequency
        else if (mel%op%species.ne.2) then
          call quit(1,'set_frequency',
     &       'no sign associated with this operator species')
        end if
c some output would be nice here (if ntest.ge.100)
c        print *,'associated frequency: ',mel%frequency
      end if

      return
      end
