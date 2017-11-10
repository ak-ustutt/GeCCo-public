*----------------------------------------------------------------------*
!>    deallocates all resources owned by a vectorspace
!!
!!    important the me_lists the vectorspace keeps are not disposed
*----------------------------------------------------------------------*
      subroutine vecsp_del(vecspace)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vecsp_init"

      type(vector_space_t),intent(inout)::
     &     vecspace
      
      type(filinf),pointer::
     &     fhand
      integer::
     &     ii
      vecspace%me_lists=>null()
      
      do ii=1,vecspace%nlists
         fhand=> vecspace%vectors(ii)%fhand
         if (ntest.ge.100)  then     !keep files in debugging mode
            if (fhand%unit.gt.0)
     &           call file_close_keep(fhand)
            !else continue
            
         else
            if (fhand%unit.gt.0)then
               call file_close_delete(fhand)
            else
               call file_delete(fhand)
            end if
         end if
         deallocate (vecspace%vectors(ii)%fhand)
      end do
      vecspace%maxvec=0
      vecspace%nvec=0
      vecspace%nlists=0
      vecspace%ivec=0
      end subroutine
