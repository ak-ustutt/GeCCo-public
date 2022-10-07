*----------------------------------------------------------------------*
      subroutine push_result_to_file(message,mel,mode,
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
      integer ::
     &     nblk, idoff
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
          idoff =mel%fhand%length_of_record*(mel%fhand%current_record-1)
          call get_vec(mel%fhand,value,idoff+1,idoff+1)
        else
          value = xnormop(mel)
        end if

        write(lures,'(a,X,'//trim(mode(5:))//')')
     &       trim(message),value

      case('LIST')

        write(lures,'(a)') trim(message)

        call wrt_mel_file(lures,5,
     &       mel,
     &       1,mel%op%n_occ_cls,
     &       str_info,orb_info)

      case('BLKS')

        write(lures,'(a)') trim(message)

        nblk = mel%op%n_occ_cls
        allocate(xnrm(nblk))
        call blk_norm_for_list(xnrm,nblk,mel,mel)
        call print_mel_contribs(lures,mel,'e10.4',xnrm,nblk)
        deallocate(xnrm)

      case default

        call quit(0,'push_result','unknown mode: '//trim(mode))

      end select

      if (openit) call file_close_keep(mel%fhand)

      return
      end

