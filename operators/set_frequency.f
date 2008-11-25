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
      include 'stdunit.h'

      type(me_list), intent(inout) ::
     &     mel

      integer, intent(in) ::
     &     order

      integer, parameter ::
     &     maximum_order = 10, ntest = 00

      real(8) ::
     &     freq(maximum_order)

      integer ::
     &     ii, iprint

      iprint = max(iprlvl, ntest)

      if (order.gt.0) then

        ! get complete user defined frequency array, sum of frequencies is zero
        call get_argument_value('calculate.experimental','freq',
     &                          xarr=freq(1:maximum_order))
        freq(order:maximum_order) = 0d0
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

        if (iprint.ge.10) 
     &        write(luout,*)
     &             'Frequency associated with ',trim(mel%label),
     &             ': ',mel%frequency
      end if

      return
      end
