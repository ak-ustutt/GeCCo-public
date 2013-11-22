*----------------------------------------------------------------------*
      subroutine reset_file_rec(fhand)
*----------------------------------------------------------------------*
*     mark the "current_record" in file fhand as not up to date
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      type(filinf), intent(inout) ::
     &     fhand
      
      if (.not.associated(fhand%last_mod))
     &     call quit(1,'reset_file_rec',
     &     'last_mod array not initialized')

      if (fhand%current_record.lt.fhand%active_records(1).or.
     &    fhand%current_record.gt.fhand%active_records(2)) then
        write(lulog,*) 'file: ',trim(fhand%name)
        write(lulog,*) 'low, current, high: ',
     &       fhand%active_records(1),
     &       fhand%current_record,
     &       fhand%active_records(2)
        call quit(1,'reset_file_rec','record info is inconsistent')
      end if

      fhand%last_mod(fhand%current_record) = -1

      return
      end
