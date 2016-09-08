*----------------------------------------------------------------------*
!>    initializes and allocates all fields of a vectorspace
!!
!!    @param vecsp  the vectorspace
!!    @param maxdim maximum dimensions of the vectorspace
!!    @param me_lists(nlists)  me_lists for shape definition of all lists within a vector
!!    @param name a name for the vectorspace
*----------------------------------------------------------------------*
      subroutine vecsp_init(vecsp, maxdim, me_lists, nlists, name)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      include 'ioparam.h'

      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vecsp_init"


      character(len=*),intent(in)::
     &     name
      
      type(vector_space_t),intent(inout)::
     &     vecsp

      integer,intent(in)::
     &     maxdim,
     &     nlists

      type(me_list_array), target,dimension(nlists), intent(in)::
     &     me_lists
      type(filinf), pointer ::
     &     fhand
      type(me_list),pointer::
     &     mel
      integer::
     &     ii, jj

      jj=0
      vecsp%me_lists=> me_lists
      
      allocate(vecsp%vectors(nlists))
      do ii=1,nlists
         allocate(vecsp%vectors(ii)%fhand)
         fhand => vecsp%vectors(ii)%fhand
         mel=>  vecsp%me_lists(ii)%mel

         call file_init(fhand,
     &        trim(dvd_get_filename(ii, name, jj)),
     &        ftyp_da_unf, lblk_da  )

         fhand%length_of_record=
     &     ((mel%len_op-1)/fhand%reclen+1)*fhand%reclen
         fhand%active_records(1)=1
         fhand%active_records(2)=maxdim
         fhand%current_record=1
         allocate(fhand%last_mod(maxdim))
         fhand%last_mod(1:maxdim) = -1
         call file_open(fhand)
      end do
      
      vecsp%maxvec=maxdim
      vecsp%nvec=0
      vecsp%nlists=nlists
      vecsp%ivec=1
      return
      contains 
*----------------------------------------------------------------------*
!>   builds a non existing filename for a vectorspace
!!
!!
*----------------------------------------------------------------------*
      function dvd_get_filename(ii, name, jj)
*----------------------------------------------------------------------*
      character(len=maxfilnam)::       !Note def_filinf declares maxfilenam this should be smaller
     &     dvd_get_filename  
      integer,intent(in)::
     &     ii

      character(len=*),intent(in)::
     &     name
      integer,intent(inout)::jj

      logical::
     &     file_exists

      if (ii.ge.100) call quit(1,i_am,
     &     "can only generate up to 100 names for a vectorspace.")
      if (jj.ge. 100)call quit(1,i_am,
     &     "can only generate up to 100 names"//
     &     " for similar vectorspaces.")


      if (jj.eq.0)then
         write(dvd_get_filename, fmt='(A,"_",I2.2,".da")') trim(name),ii !try shortened form with name as unique identifier
      else 
         write(dvd_get_filename,
     &        fmt='(A,I2.2,"_",I2.2,".da")') trim(name),jj,ii
      end if

      inquire(file=trim(dvd_get_filename),exist=file_exists)
      do while (file_exists)
         jj=jj+1
         write(dvd_get_filename,
     &        fmt='(A,I2.2,"_",I2.2,".da")') trim(name),jj,ii
        if (jj.ge. 100)call quit(1,i_am,
     &        "can only generate up to 100 names"//
     &        " for similar vectorspaces.")
        inquire(file=trim(dvd_get_filename),exist=file_exists)     
      end do

      end function
      end subroutine 




