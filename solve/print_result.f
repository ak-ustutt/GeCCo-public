*----------------------------------------------------------------------*
      subroutine print_result(maxord,order,ifreq,mel,zero)
*----------------------------------------------------------------------*
*     print calculated property
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'
      include 'stdunit.h'
      include 'def_me_list.h'
      include 'def_filinf.h'
      include 'ifc_input.h'

      type(me_list), intent(in), optional ::
     &     mel

      integer, intent(in) ::
     &     maxord, order, ifreq(order)

      logical, intent(in) ::
     &     zero

      integer, parameter ::
     &     len_short = 20, maximum_order = 10, len_long = 200

      character(len_short) ::
     &     pert

      character(len_long) ::
     &     part1, part2, part3, part4, line

      integer ::
     &     freq_idx, last1, last2, last3, last4

      real(8) ::
     &     freq(maximum_order), abssum, value

      logical ::
     &     closeit

      pert(1:len_short) = ' '
      call get_argument_value('calculate.experimental','pert',
     &                        str=pert(1:len_short))

      ! duplicate values for pert if necessary and not specified
      do freq_idx = 2,maxord
        if (pert(freq_idx:freq_idx).eq.' ')
     &      pert(freq_idx:freq_idx) = pert(freq_idx-1:freq_idx-1)
      end do

      if (order.gt.0) then
        ! get complete user defined frequency array, sum of frequencies is zero
        call get_argument_value('calculate.experimental','freq',
     &                          xarr=freq(1:maximum_order))
        freq(maxord:maximum_order) = 0d0
        freq(maxord) = -sum(freq)
      end if

      if (order.eq.0) then
        last1 = 6
        part1(1:last1) = 'Energy'
      else if (order.eq.1) then
        last1 = 41
        part1(1:last1) = 'Dipole moment (correlation contribution) '
      else if (order.eq.2) then
        last1 = 15
        part1(1:last1) = 'Polarizability '
      else if (order.eq.3) then
        last1 = 20
        part1(1:last1) = 'Hyperpolarizability '
      else if (order.eq.4) then
        last1 = 24
        part1(1:last1) = '2nd Hyperpolarizability '
      else if (order.eq.5) then
        last1 = 24
        part1(1:last1) = '3rd Hyperpolarizability '
      else
        last1 = 24
        part1(1:last1) = 'xth Hyperpolarizability '
        write(part1(1:1),'(i1)') order
      end if

      if (order.gt.0) then
        part2(1:2) = '<<'
        do freq_idx = 1,order
          part2(2*freq_idx+1:2*freq_idx+2) =
     &                    pert(ifreq(freq_idx):ifreq(freq_idx))//','
        end do
        part2(2*order+2:2*order+3) = '>>'
        last2 = 2*order+3

        part3(1:19) = ' with frequencies ('
        abssum = 0d0
        do freq_idx = 1,order
          abssum = abssum + abs(freq(ifreq(freq_idx)))
          write(part3(10*freq_idx+10:10*freq_idx+18),'(f9.6)')
     &            freq(ifreq(freq_idx))
          part3(10*freq_idx+19:10*freq_idx+19) = ','
        end do
        if (abssum.gt.0d0) then
          part3(10*order+19:10*order+19) = ')'
          last3 = 10*order+19
        else
          last3 = 9
          part3(1:last3) = ' (static)'
        end if
      end if

      ! read result from file if not zero by symmetry
      value = 0d0
      if (.not.zero) then
        closeit = .false.
        if (mel%fhand%unit.lt.0) then
          call file_open(mel%fhand)
          closeit = .true.
        end if
        call get_vec(mel%fhand,value,1,1)
        if (closeit)
     &       call file_close_keep(mel%fhand)
      end if

      last4 = 28
      write(part4(1:last4),'(": >>> ",f18.12," <<<")') value

      ! print result
      line(1:79) = '>============================================='//
     &             '===============================<'
      write(luout,*)
      write(luout,*) line(1:79)
      if (order.eq.0) then
        write(luout,*) part1(1:last1)//part4(1:last4)
      else if (order.eq.1) then
        write(luout,*) part1(1:last1)//part2(1:last2)//part4(1:last4)
      else
        write(luout,*) part1(1:last1)//part2(1:last2)//part3(1:last3)//
     &                 part4(1:last4)
      end if
      write(luout,*) line(1:79)
      write(luout,*)

      return
      end

