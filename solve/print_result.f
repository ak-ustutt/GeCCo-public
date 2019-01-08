*----------------------------------------------------------------------*
      subroutine print_result(order,ifreq,mel,zero,orb_info)
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
      include 'def_orbinf.h'
      include 'def_pert_info.h'
      include 'multsign.h'

      type(me_list), intent(in) ::
     &     mel
      type(orbinf), intent(in) ::
     &     orb_info

      integer, intent(in) ::
     &     order, ifreq(order)

      logical, intent(in) ::
     &     zero

      integer, parameter ::
     &     len_long = 200

      character(len_long) ::
     &     part1, part2, part3, part4, part5, line, phrase

      character(len=4) ::
     &     property

      integer ::
     &     idx, last1, last2, last3, last4, last5, iostatus,
     &     icnt, ncnt, npop, ncmp, glob_sign, temp_sign

      real(8) ::
     &     abssum, value, nucmom, fac

      type(pert_op_info) ::
     &     pop(3*maxpop)

      type(pert_component_info) ::
     &     cmp(maxcmp)

      integer, allocatable ::
     &     maxord(:)

      real(8), allocatable ::
     &     freq(:,:)

      logical ::
     &     closeit, file_exists, got_nucmom

      ncnt = is_keyword_set('method.response')
      allocate(maxord(ncnt),freq(ncnt,maximum_order))

      do icnt = 1,ncnt
        call get_argument_value('method.response','order',
     &       keycount=icnt,ival=maxord(icnt))
      end do

      call get_response_input(ncnt,maxord,ncmp,cmp,npop,pop,orb_info)

      ! find out global sign and factor
      glob_sign = 1
      fac = 1d0
      do idx = 1, order
        temp_sign = pop(cmp(ifreq(idx))%pop_idx)%sign
        fac = fac * 1d0/real((temp_sign+3)/4)
        temp_sign = mod(temp_sign-1,4)+1
        glob_sign = multsign(glob_sign,temp_sign)
      end do

      if (order.gt.0) then
        part2(1:2) = '<<'
        last2 = 0
        if (glob_sign.eq.3.or.glob_sign.eq.4) then
          part2(1:5) = 'Im(<<'
          last2 = 3
        end if
        do idx = 1,order
          part2(3*idx+last2:3*idx+2+last2) =
     &                    pop(cmp(ifreq(idx))%pop_idx)%name//
     &                    pop(cmp(ifreq(idx))%pop_idx)%comp//','
          part5(idx:idx) = pop(cmp(ifreq(idx))%pop_idx)%comp
        end do
        part2(5+last2:5+last2) = ';'
        part2(3*order+2+last2:3*order+3+last2) = '>>'
        last2 = 3*order+3+last2
        if (glob_sign.eq.3.or.glob_sign.eq.4) then
          part2(last2+1:last2+1) = ')'
          last2 = last2+1
        end if
        last5 = order
        ! only one ‹› for dipole moment
        if (order.eq.1) then
          if (glob_sign.eq.3.or.glob_sign.eq.4) then
            part2(last2-5:last2-3) = part2(last2-4:last2-2)
            part2(last2-2:last2-2) = ')'
          else
            part2(last2-4:last2-2) = part2(last2-3:last2-1)
          end if
          last2 = last2 - 2      
        end if

        part3(1:1) = '('
        abssum = 0d0
        do idx = 2,order
          abssum = abssum + abs(cmp(ifreq(idx))%freq)
          write(part3(10*idx-18:10*idx-10),'(f9.6)')
     &            cmp(ifreq(idx))%freq
          part3(10*idx-9:10*idx-9) = ','
        end do
        if (abssum.gt.0d0) then
          part3(10*order-9:10*order-9) = ')'
          last3 = 10*order-9
        else
          last3 = 3
          part3(1:last3) = '(0)'
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

      ! flip sign if - or -i
      if (glob_sign.eq.2.or.glob_sign.eq.4) value = -value
      value = fac * value

      last4 = 28
      if (abs(value).lt.100000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(" = >>",f19.12," << ")') value
      else if (abs(value).lt.10000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(" = >>",f19.10," << ")') value
      else if (abs(value).lt.1000000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(" = >>",f19.8," << ")') value
      else if (abs(value).lt.100000000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(" = >>",f19.6," << ")') value
      else
        write(part4(1:last4),'(" = >>",E19.10," << ")') value
      end if

      ! print response function
      line(1:79) = '>============================================='//
     &             '===============================<'
      write(luout,*)
      write(luout,*) line(1:79)
      if (order.eq.1) then
        write(luout,*) part2(1:last2)//part4(1:last4)
      else if (order.ge.2) then
        write(luout,*) part2(1:last2)//part3(1:last3)//part4(1:last4)
      end if

      ! try to recognize property and print as such
      property = '    '
      if (all(pop(cmp(ifreq(1:order))%pop_idx)%name.eq.'d')) then
        if (mod(order,2).eq.0.and.order.gt.0) value = -value
        property = 'elec'
      end if
      if (all(pop(cmp(ifreq(1:order))%pop_idx)%name.eq.'r').and.
     &    order.gt.0) then
        value = -value
        property = 'elec'
      end if
      if (order.eq.0) property = 'ener'

      if (order.eq.1) then
        got_nucmom = .false.
        nucmom = 0d0
        inquire(file='4traf.out',exist=file_exists)
        if (file_exists) then
          ! read nuclear dipole moment from file 4traf.out if existing
          open(1,file='4traf.out')
          do while (iostatus.ge.0)
            read(1,'(a25)',iostat=iostatus) phrase(1:25)
            if (phrase(1:25).eq.' dipole moment of nuclei:') then
              if (pop(cmp(ifreq(1))%pop_idx)%comp.ne.'X')
     &                   read(1,*) phrase(1:25)
              if (pop(cmp(ifreq(1))%pop_idx)%comp.eq.'Z')
     &                   read(1,*) phrase(1:25)
              read(1,*) phrase, nucmom
              got_nucmom = .true.
              write(luout,'(" adding nuclear dipole moment <<",a1,
     &                       "nuc>> = ",f18.12)')
     &                       pop(cmp(ifreq(1))%pop_idx)%comp, nucmom
              exit
            end if
          end do
          close(1)
        end if
        value = value + nucmom
      end if

      last4 = 28
      if (abs(value).lt.100000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(": >>>",f19.12," <<<")') value
      else if (abs(value).lt.10000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(": >>>",f19.10," <<<")') value
      else if (abs(value).lt.1000000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(": >>>",f19.8," <<<")') value
      else if (abs(value).lt.100000000000d0.and.abs(value).gt.1d-1) then
        write(part4(1:last4),'(": >>>",f19.6," <<<")') value
      else
        write(part4(1:last4),'(": >>>",E19.10," <<<")') value
      end if

      select case(property)
      case('ener')
        write(luout,*) 'Energy'//part4(1:last4)
      case('elec')
      if (order.eq.1) then
        if (got_nucmom) then
          last1 = 22
          part1(1:last1) = 'Dipole moment (total) '
        else
          last1 = 27
          part1(1:last1) = 'Dipole moment (electronic) '
        end if
        write(luout,*) part1(1:last1)//part5(1:last5)//part4(1:last4)
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
      if (order.ge.2) then
        write(luout,*) part1(1:last1)//part5(1:last5)//part3(1:last3)//
     &                 part4(1:last4)
      end if
      end select
          

      write(luout,*) line(1:79)
      write(luout,*)

      return
      end
