*----------------------------------------------------------------------*
!>    a subroutine to orthogonalize a vector with respect to a subspace( assumed to be already orthogonal)
!!
!!    @param vecsp the vectorspace
!!    @param me_lists the lists of the vector
!!    @param nlists the number of lists in the vector     
*----------------------------------------------------------------------*
      subroutine vecsp_orthvec(vecsp, me_lists, nlists,
     $     xbuf1, xbuf2, lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      
      
      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vecsp_orthvec"

      integer,intent(in)::
     &     nlists, lbuf
      type(vector_space_t), intent(inout)::
     &     vecsp
      type(me_list_array),intent(in)::
     &     me_lists(*)
      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*)
      
      real(8)::
     &     xnorm,
     &     fac
      real(8),allocatable::
     &     contrib(:)

      integer::
     &     nvec,
     &     ivec,
     &     ilist,
     &     lenlist,
     &     irec,
     &     nblock,
     &     njoined

      type(filinf),pointer::
     &     ffme,ffvec
      real(8),external::
     &     da_ddot
      
      nvec=vecsp%nvec
      allocate(contrib(1:nvec))
      contrib(1:nvec)=0d0

      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         ffvec => vecsp%vectors(ilist)%fhand
         lenlist= me_lists(ilist)%mel%len_op
         do ivec=1,nvec
            contrib(ivec)=contrib(ivec)
     &           +da_ddot(ffme, irec,ffvec,  ivec, lenlist,
     &           xbuf1, xbuf2, lbuf)
            
         end do 
      end do
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         ffvec => vecsp%vectors(ilist)%fhand
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf1,lenlist)
         do ivec=1,nvec
            call  vecsp_get_list_buf(vecsp, ivec, ilist, lenlist,
     &           xbuf2, lbuf)
            xbuf1(1:lenlist)=xbuf1(1:lenlist)
     &           -contrib(ivec)*xbuf2(1:lenlist)
         end do
         call vec_to_da(ffme,irec,xbuf1,lenlist)
      end do
      end subroutine
