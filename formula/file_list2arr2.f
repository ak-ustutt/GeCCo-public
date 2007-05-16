*----------------------------------------------------------------------*
      subroutine file_list2arr2(fl_list,fl_arr,nfiles)
*----------------------------------------------------------------------*
*     set up an array of pointers that point to the entries of the
*     list on fl_list
*     new version with a real pointer array
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_list.h'
      include 'def_file_array.h'

      type(file_list), intent(in), target ::
     &     fl_list
      integer, intent(in) ::
     &     nfiles
      type(file_array), intent(out) ::
     &     fl_arr(nfiles)

      type(file_list), pointer ::
     &     current

      integer ::
     &     ifile

      current => fl_list

      do ifile = 1, nfiles
        if (.not.associated(current%fhand))
     &       call quit(1,'file_list2arr',
     &                   'unallocated file handle on list')
        fl_arr(ifile)%fhand => current%fhand
        if (ifile.lt.nfiles.and..not.(associated(current%next)))
     &       call quit(1,'file_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
