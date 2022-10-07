*----------------------------------------------------------------------*
      subroutine print_list(message,mel,mode,check,expected,
     &                      orb_info,str_info)
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

      real(8), parameter ::
     &     thresh_warn = 1d-10

      type(me_list), intent(inout) ::
     &     mel
      character(len=*), intent(in) ::
     &     message, mode
      real(8), intent(in) ::
     &     check, expected
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info

      logical ::
     &     openit
      integer ::
     &     nblk, idisc_off
      real(8) ::
     &     value
      type(filinf) ::
     &     ffout

      real(8), pointer ::
     &     xnrm(:)

      real(8), external ::
     &     xnormop

      ! open list (if necessary)
      openit = mel%fhand%unit.lt.0.and..not.mel%fhand%buffered
      if (openit) call file_open(mel%fhand)

      select case(mode(1:4))

      case('NORM','SCAL')
        if (mode(1:4).eq.'SCAL') then
          idisc_off = mel%fhand%length_of_record*
     &               (mel%fhand%current_record-1)
          call get_vec(mel%fhand,value,idisc_off+1,idisc_off+1)
        else
          value = xnormop(mel)
        end if

        write(luout,'(x,a,'//trim(mode(5:))//',f20.12)')
     &       trim(message),value
!       write(luout,*)
!    &       trim(message),value

        ! check result (if requested)
        if (check.gt.0) then
          if (abs(value-expected).gt.check) then
            write(luout,'(x,a,'//trim(mode(5:))//')')
     &       'We had expected the result: ',expected
            call quit(1,'print_list','Not the result we want!')
          else if (abs(value-expected).gt.thresh_warn) then
             write(luout,'(x,a,'//trim(mode(5:))//')')
     &       'We had expected the result: ',expected
            call warn('print_list','Deviation from expected result.')
          end if
        end if

      case('LIST')

        write(luout,'(x,a)') trim(message)

        call wrt_mel_file(luout,5,
     &       mel,
     &       1,mel%op%n_occ_cls,
     &       str_info,orb_info)

      case('BLKS')

        write(luout,'(x,a)') trim(message)

        nblk = mel%op%n_occ_cls
        allocate(xnrm(nblk))
        call blk_norm_for_list(xnrm,nblk,mel,mel)
        call print_mel_contribs(luout,mel,'e10.4',xnrm,nblk)
        deallocate(xnrm)

c        call wrt_mel_file(luout,1,
c     &       mel,
c     &       1,mel%op%n_occ_cls,
c     &       str_info,orb_info)

      case default

        write(luout,'(x,a)') trim(message)
        write(luout,'(x,a,a)') 'Written to file named ',trim(mode)

        call file_init(ffout,trim(mode),ftyp_sq_frm,0)
        call file_open(ffout)
        call wrt_mel_to_file(ffout%unit,
     &       mel,
     &       1,mel%op%n_occ_cls,
     &       str_info,orb_info)
        call file_close_keep(ffout)

c        call quit(0,'print_list','unknown mode: '//trim(mode))

      end select

      if (openit) call file_close_keep(mel%fhand)

      return
      end

