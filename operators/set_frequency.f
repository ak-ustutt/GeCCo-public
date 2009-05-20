*----------------------------------------------------------------------*
      subroutine set_frequency(mel,freq)
*----------------------------------------------------------------------*
*     set frequency assigned to ME-list
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_me_list.h'
      include 'stdunit.h'

      type(me_list), intent(inout) ::
     &     mel
      real(8), intent(in) ::
     &     freq

      integer, parameter ::
     &     ntest = 00

      integer ::
     &     iprint

      iprint = max(iprlvl, ntest)


      mel%frequency = freq

      ! flip sign if operator species = 1 (T-amplitudes)
      if (mel%op%species.eq.1) then
        mel%frequency = -mel%frequency
      else if (mel%op%species.ne.2) then
        call quit(1,'set_frequency',
     &     'no sign associated with this operator species')
      end if

      if (iprint.ge.10) 
     &      write(luout,*)
     &           'Frequency associated with ',trim(mel%label),
     &           ': ',mel%frequency

      return
      end
