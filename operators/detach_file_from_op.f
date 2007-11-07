*----------------------------------------------------------------------*
      subroutine detach_file_from_op(idxop,remove,op_info)
*----------------------------------------------------------------------*
*     remove association of a file with operator #idxop
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'par_globalmarks.h'
      include 'par_opnames_gen.h'
      include 'ifc_memman.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'

      integer, intent(in) ::
     &     idxop
      logical, intent(in) ::
     &     remove
      type(operator_info), intent(inout) ::
     &     op_info

      type(filinf), pointer ::
     &     ffop
      type(operator), pointer ::
     &     op

      op => op_info%op_arr(idxop)%op
      ffop => op_info%opfil_arr(idxop)%fhand

      if (.not.associated(ffop))
     &     call quit(1,'detach_file_from_op',
     &     'operator is not assigned a file')

      if (iprlvl.ge.10) then
        write(luout,'(3x,4a)')
     &       'detaching file: ',trim(ffop%name),' from operator ',
     &       trim(op%name)
      end if

      if (remove) then
        call file_delete(ffop)

        if (ffop%buffered)
     &       call quit(1,'detach_file_from_op',
     &       'buffering not yet considered!')

        deallocate(ffop)
      end if

      op_info%opfil_arr(idxop)%fhand => null()

      return

      end
