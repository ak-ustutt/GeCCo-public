*----------------------------------------------------------------------*
      subroutine set_frequency(mel)
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

      integer, parameter ::
     &     maximum_order = 10, ntest = 00

      real(8), allocatable ::
     &     freq(:,:)

      integer ::
     &     ii, iprint, icnt, ncnt, pos

      integer, allocatable ::
     &     maxord(:)

      iprint = max(iprlvl, ntest)


      ! get complete user defined frequency array, sum of frequencies is zero
      ncnt = is_keyword_set('calculate.experimental')
      allocate(freq(ncnt,maximum_order),maxord(ncnt))
      freq = 0d0
      do icnt = 1,ncnt
        call get_argument_value('calculate.experimental','order',
     &       keycount=icnt,ival=maxord(icnt))
        if (maxord(icnt).gt.0) then
          call get_argument_value('calculate.experimental','freq',
     &         keycount=icnt,xarr=freq(icnt,1:maximum_order))
          freq(icnt,maxord(icnt):maximum_order) = 0d0
          freq(icnt,maxord(icnt)) = -sum(freq(icnt,:))
        end if
      end do

      ! frequency is sum of frequencies associated with frequency indices
      if (.not.associated(mel%op%ifreq)) call quit(1,'set_frequency',
     &     'no frequency index associated to operator')
      mel%frequency = 0d0
      do ii = 1,mel%op%order
        icnt = (mel%op%ifreq(ii)-1)/maxval(maxord)+1
        pos = mel%op%ifreq(ii)-(icnt-1)*maxval(maxord)
        mel%frequency = mel%frequency + freq(icnt,pos)
      end do

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
