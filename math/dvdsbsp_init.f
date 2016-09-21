*----------------------------------------------------------------------*
!>    initializes and allocates all fields of a vectorspace
!!
!!    @param dvdsp a davidson subspace (will be reset if already initialized)
!!    @param maxsub maximum dimension of the subspaces
!!    @param me_lists(nlists) me_lists for shape information of the lists a reference will be kept, will not be deleted upon deletion of the subspace
*----------------------------------------------------------------------*
      subroutine dvdsbsp_init(dvdsbsp, maxsub, me_lists, nlists, use_s)
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
      logical::
     &     use_s(nlists)
      
      real(8),dimension(:),pointer::
     &     mymat
      integer::
     &     ii,jj,
     &     ifree
      logical::
     &     use_any_s
      
      type(me_list_array),dimension(nlists), intent(in)::
     &     me_lists


      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.ge.100)then
         write (lulog,*)"me_list template:"
         do ii=1,nlists
            write (lulog,*)me_lists(ii)%mel%label
         end do
      end if
      dvdsbsp%with_metric=.false.
      do ii=1,nlists
         dvdsbsp%with_metric=dvdsbsp%with_metric.or.use_s(ii)
      end do
      if(associated(dvdsbsp%vMv_mat))then
         call warn(i_am,
     &        "initiating an already initiated davidson subspace")

         call dvdsbsp_del(dvdsbsp)
      end if
      
      call vecsp_init(dvdsbsp%vspace, maxsub, me_lists, nlists, 
     &     "dvd_vspace")
      call vecsp_init(dvdsbsp%Mvspace, maxsub, me_lists, nlists, 
     &     "dvd_Mvspace")
      if (dvdsbsp%with_metric)then
         call vecsp_init(dvdsbsp%Svspace, maxsub, me_lists, nlists, 
     &        "dvd_Svspace")
      end if
      
      dvdsbsp%nmaxsub=maxsub
      dvdsbsp%ncursub=0
      dvdsbsp%icursub=0
      dvdsbsp%lcursub=0
      ifree=mem_alloc_real(dvdsbsp%vMv_mat, maxsub*maxsub,
     &     "davidson subspace vMv")
      dvdsbsp%vMv_mat(1:maxsub*maxsub) = 0d0
      if (dvdsbsp%with_metric)then
         ifree=mem_alloc_real(dvdsbsp%vSv_mat, maxsub*maxsub,
     &        "davidson subspace vSv")
         dvdsbsp%vSv_mat(1:maxsub*maxsub) = 0d0
      end if 


      end subroutine
