*----------------------------------------------------------------------*
      subroutine switch_mel_record(mel,rec)
*----------------------------------------------------------------------*
*     safeguarded change of record on ME-list file
*     use only this routine and do not change contents of mel
*     directly
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'

      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     rec
      
      type(filinf), pointer ::
     &     ffop

      ffop => mel%fhand

      if (rec.lt.ffop%active_records(1).or.
     &    rec.gt.ffop%active_records(2)) then
        write(lulog,*) 'list, operator, file: ',
     &       trim(mel%label),
     &       trim(mel%op%name),
     &       trim(ffop%name)
        write(lulog,*) 'bounds:  ',
     &       ffop%active_records(1:2)
        write(lulog,*) 'request: ',rec
        call quit(1,'switch_mel_record',
     &       'requested record out of bounds')
      end if

      if (ffop%buffered)
     &     call quit(1,'switch_mel_record',
     &     'switching and buffering: did you take care of that?')

c dbg
      print *,'file: ',trim(ffop%name)
      print *,'switching to rec. ',rec
c dbg
      ffop%current_record = rec

      return
      end
