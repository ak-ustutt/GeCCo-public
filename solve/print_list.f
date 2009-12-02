*----------------------------------------------------------------------*
      subroutine print_list(message,mel,mode,orb_info,str_info)
*----------------------------------------------------------------------*
*     print list + message
*     andreas, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'
      include 'stdunit.h'
      include 'def_me_list.h'
      include 'def_filinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_input.h'

      type(me_list), intent(inout) ::
     &     mel
      character(len=*), intent(in) ::
     &     message, mode
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info

      logical ::
     &     openit
      real(8) ::
     &     value

      real(8), external ::
     &     xnormop

      ! open list (if necessary)
      openit = mel%fhand%unit.lt.0
      if (openit) call file_open(mel%fhand)

      select case(mode(1:4))

      case('NORM','SCAL')
        if (mode(1:4).eq.'SCAL') then
          call get_vec(mel%fhand,value,1,1)
        else
          value = xnormop(mel)
        end if

        write(luout,'(x,a,'//trim(mode(5:))//')')
     &       trim(message),value

      case('LIST')

        write(luout,'(x,a)') trim(message)

        call wrt_mel_file(luout,5,
     &       mel,
     &       1,mel%op%n_occ_cls,
     &       str_info,orb_info)

      case default
        call quit(0,'print_list','unknown mode: '//trim(mode))

      end select

      if (openit) call file_close_keep(mel%fhand)

      return
      end

