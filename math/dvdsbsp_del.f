*----------------------------------------------------------------------*
!>    disposes of a davidson subspace 
*----------------------------------------------------------------------*
      subroutine dvdsbsp_del(dvdsbsp)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      include 'ifc_memman.h'

      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="dvdsbsp_del"

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer::
     &     ifree

      
      call vecsp_del(dvdsbsp%vspace)
      call vecsp_del(dvdsbsp%Mvspace)
      if(dvdsbsp%with_metric)then
         call vecsp_del(dvdsbsp%Svspace)
      end if
      dvdsbsp%with_metric=.false.
      dvdsbsp%nmaxsub=0
      dvdsbsp%ncursub=0
      dvdsbsp%icursub=0
      deallocate(dvdsbsp%use_metric)
      end subroutine
