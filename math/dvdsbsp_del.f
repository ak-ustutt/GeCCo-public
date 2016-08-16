*----------------------------------------------------------------------*
!>    disposes of a davidson subspace without memory leaks
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

      ifree=mem_dealloc("davidson subspace vMv")
      call vecspace_del(dvdsbsp%vspace)
      call vecspace_del(dvdsbsp%Mvspace)
      dvdsbsp%nmaxsub=0
      dvdsbsp%ncursub=0
      dvdsbsp%icursub=0

      end subroutine
