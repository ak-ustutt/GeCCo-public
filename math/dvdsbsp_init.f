*----------------------------------------------------------------------*
!>    initializes and allocates all fields of a vectorspace
!!
!!    @param dvdsp a davidson subspace (will be reset if already initialized)
!!    @param maxsub maximum dimension of the subspaces
!!    @param me_lists(nlists) me_lists for shape information of the lists a reference will be kept, will not be deleted upon deletion of the subspace
*----------------------------------------------------------------------*
      subroutine dvdsbsp_init(dvdsbsp, maxsub, me_lists, nlists)
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
     &     i_am="dvdsbsp_init"

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      integer, intent(in) ::
     &     maxsub,              !maximum of vectors possible in this subspace
     &     nlists               !number of lists a vector consists of
      real(8),dimension(:),pointer::
     &     mymat
      integer::
     &     ii,jj,
     &     ifree
      intrinsic::
     &     reshape

      type(me_list_array),dimension(nlists), intent(in)::
     &     me_lists

      if(associated(dvdsbsp%vMv_mat))then
         call warn(i_am,
     &        "initiating an already initiated davidson subspace")

         call dvdsbsp_del(dvdsbsp)
      end if
      
      call vecsp_init(dvdsbsp%vspace, maxsub, me_lists, nlists, 
     &     "dvd_vspace")
      call vecsp_init(dvdsbsp%Mvspace, maxsub, me_lists, nlists, 
     &     "dvd_Mvspace")

      dvdsbsp%nmaxsub=maxsub
      dvdsbsp%ncursub=0
      dvdsbsp%icursub=0

      ifree=mem_alloc_real(mymat, maxsub*maxsub, 
     &     "davidson subspace vMv")
      do ii=1,maxsub*maxsub
         mymat(ii)=0
      end do
      dvdsbsp%vMv_mat(1:maxsub, 1:maxsub)=> mymat
      mymat=>null()



      end subroutine
