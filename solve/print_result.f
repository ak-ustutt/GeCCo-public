*----------------------------------------------------------------------*
      subroutine print_result(order,ifreq,mel,zero)
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
     &     order, ifreq(order)

      logical, intent(in) ::
     &     zero

      integer, parameter ::
     &     maximum_order = 10, len_long = 200

      character(len_long) ::
     &     pert, part1, part2, part3, part4, line, phrase

      integer ::
     &     freq_idx, last1, last2, last3, last4, iostatus,
     &     icnt, ncnt, pos

      real(8) ::
     &     abssum, value, nucmom

      integer, allocatable ::
     &     maxord(:)

      real(8), allocatable ::
     &     freq(:,:)

      logical ::
     &     closeit, file_exists, got_nucmom

      ncnt = is_keyword_set('calculate.experimental')
      allocate(maxord(ncnt),freq(ncnt,maximum_order))

      pert(1:len_long) = ' '
      do icnt = 1,ncnt
        call get_argument_value('calculate.experimental','order',
     &       keycount=icnt,ival=maxord(icnt))
      end do
      do icnt = 1,ncnt
        pos = (icnt-1)*maxval(maxord) + 1
        call get_argument_value('calculate.experimental','pert',
     &       keycount=icnt,str=pert(pos:len_long))

        ! duplicate values for pert if necessary and not specified
        do freq_idx = 2,maxord(icnt)
          pos = (icnt-1)*maxval(maxord) + freq_idx
          if (pert(pos:pos).eq.' ')
     &        pert(pos:pos) = pert(pos-1:pos-1)
        end do

        ! get complete user defined frequency array, sum of frequencies is zero
        call get_argument_value('calculate.experimental','freq',
     &       keycount=icnt,xarr=freq(icnt,1:maximum_order))
        freq(icnt,maxord(icnt):maximum_order) = 0d0
        freq(icnt,maxord(icnt)) = -sum(freq(icnt,:))
      end do

      if (order.eq.0) then
        last1 = 6
        part1(1:last1) = 'Energy'
      else if (order.eq.1) then
        got_nucmom = .false.
        nucmom = 0d0
        inquire(file='4traf.out',exist=file_exists)
        if (file_exists) then
          ! read nuclear dipole moment from file 4traf.out if existing
          open(1,file='4traf.out')
          do while (iostatus.ge.0)
            read(1,'(a25)',iostat=iostatus) phrase(1:25)
            if (phrase(1:25).eq.' dipole moment of nuclei:') then
              if (pert(ifreq(1):ifreq(1)).ne.'X') read(1,*) phrase(1:25)
              if (pert(ifreq(1):ifreq(1)).eq.'Z') read(1,*) phrase(1:25)
              read(1,*) phrase, nucmom
              got_nucmom = .true.
              write(luout,'("adding nuclear dipole moment <<",a1,
     &                       "nuc>> = ",f18.12)')
     &                       pert(ifreq(1):ifreq(1)), nucmom
              exit
            end if
          end do
          close(1)
        end if
        if (got_nucmom) then
          last1 = 22
          part1(1:last1) = 'Dipole moment (total) '
        else
          last1 = 27
          part1(1:last1) = 'Dipole moment (electronic) '
        end if
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
        write(part1(1:1),'(i1)') order-2
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
          icnt = (ifreq(freq_idx)-1)/maxval(maxord)+1
          pos = ifreq(freq_idx)-(icnt-1)*maxval(maxord)
          abssum = abssum + abs(freq(icnt,pos))
          write(part3(10*freq_idx+10:10*freq_idx+18),'(f9.6)')
     &            freq(icnt,pos)
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

      ! flip sign acc. to Taylor expansion convention: 
      ! E = E0 - mu*e - (1/2)alpha*e^2 - (1/6)beta*e^3 - ...
      if (order.gt.0) value = -value
      ! add nuclear dipole moment (if not found, add zero)
      if (order.eq.1) value = value + nucmom

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

