*----------------------------------------------------------------------*
      subroutine detach_mel_file(mel,remove)
*----------------------------------------------------------------------*
*     suspend association of a file with ME-list mel
*     remove: the file will be removed
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'par_globalmarks.h'
      include 'par_opnames_gen.h'
      include 'ifc_memman.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'

      type(me_list), intent(inout) ::
     &     mel
      logical, intent(in) ::
     &     remove

      type(filinf), pointer ::
     &     ffop
      type(operator), pointer ::
     &     op

      op => mel%op
      ffop => mel%fhand

      if (.not.associated(ffop))
     &     call quit(1,'detach_mel_file',
     &     'list is not assigned a file')

      if (iprlvl.ge.20) then
        write(luout,'(3x,7a)')
     &       'detaching file: ',trim(ffop%name),' from list ',
     &       trim(mel%label),' (operator: ',trim(op%name),')'
      end if

      if (remove) then
        call file_delete(ffop)

        if (ffop%buffered)
     &       call quit(1,'detach_file_from_op',
     &       'buffering not yet considered!')

        if (associated(ffop%last_mod))
     &       deallocate(ffop%last_mod)

        deallocate(ffop)
      end if

      mel%fhand => null()

      return

      end
