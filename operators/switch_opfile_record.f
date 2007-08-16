*----------------------------------------------------------------------*
      subroutine switch_opfile_record(idxop,rec,op_info)
*----------------------------------------------------------------------*
*     safeguarded change of record on operator file
*     use only this routine and do not change contents of op_info
*     directly
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'

      integer ::
     &     idxop, rec
      type(operator_info), intent(inout) ::
     &     op_info
      
      type(filinf), pointer ::
     &     ffop

      ffop => op_info%opfil_arr(idxop)%fhand

      if (rec.lt.ffop%active_records(1).or.
     &    rec.gt.ffop%active_records(2)) then
        write(luout,*) 'operator, file: ',
     &       trim(op_info%op_arr(idxop)%op%name),
     &       trim(ffop%name)
        write(luout,*) 'bounds:  ',
     &       ffop%active_records(1:2)
        write(luout,*) 'request: ',rec
        call quit(1,'switch_opfile_record',
     &       'requested record out of bounds')
      end if

      if (ffop%buffered)
     &     call quit(1,'switch_opfile_record',
     &     'switching and buffering: did you take care of that?')

c dbg
      print *,'file: ',trim(ffop%name)
      print *,'switching to rec. ',rec
c dbg
      ffop%current_record = rec

      return
      end
