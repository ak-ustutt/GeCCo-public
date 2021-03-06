*----------------------------------------------------------------------*
!>    assembles the residual for one root
      
*----------------------------------------------------------------------*
      subroutine dvdsbsp_assemble_residual(dvdsbsp,
     &     coeff, rval, ncoeff,
     $     me_lists, nlists, xnrm2,
     $     xbuf1, xbuf2, nincore, lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="dvdsbsp_assemble_residual"

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      integer,intent(in)::
     &     ncoeff,nlists,nincore,lbuf
      
      real(8),intent(in)::
     &     coeff(*),            !Vector with coefficients
     &     rval                 !Ritz-value
      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*)
      real(8),intent(out)::
     &     xnrm2(nlists)
      real(8),external::
     &     me_ddot
      integer::
     &     ilist
      type(me_list_array),intent(inout)::
     &     me_lists(*)              ! me-lists the vector should be assembled at

  

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      if (dvdsbsp%with_metric)then
         call vecsp_assemble_vector(dvdsbsp%Svspace,
     &        me_lists, nlists,
     &        coeff, ncoeff,
     &        0d0,
     &        xbuf1,xbuf2,nincore,lbuf)
      else
         call vecsp_assemble_vector(dvdsbsp%vspace,
     &        me_lists, nlists,
     &        coeff, ncoeff,
     &        0d0,
     &        xbuf1,xbuf2,nincore,lbuf)
      end if
      call vecsp_assemble_vector(dvdsbsp%Mvspace,
     &     me_lists, nlists,
     &     coeff, ncoeff,
     &     -rval,
     &     xbuf1,xbuf2,nincore,lbuf)
      do ilist=1,nlists
         xnrm2(ilist)=me_ddot(me_lists(ilist)%mel,me_lists(ilist)%mel ,
     &        xbuf1, xbuf2, nincore, lbuf)
      end do
*----------------------------------------------------------------------*
      end subroutine

      
      

      


      
