*----------------------------------------------------------------------*
!>    assembles the residual for one root
      
*----------------------------------------------------------------------*
      subroutine dvdsbsp_assemble_residual(dvdsbsp,
     &     coeff, rval, ncoeff,
     $     me_lists, nlists, xnrm,
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
     &     xnrm
 
      type(me_list_array),intent(inout)::
     &     me_lists(*)              ! me-lists the vector should be assembled at

  

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      
      call vecsp_assemble_vector(dvdsbsp%vspace,
     &     me_lists, nlists,
     &     coeff, ncoeff,
     &     0d0,
     &     xbuf1,xbuf2,nincore,lbuf)
      
      call vecsp_assemble_vector(dvdsbsp%Mvspace,
     &     me_lists, nlists,
     &     coeff, ncoeff,
     &     -rval,
     &     xbuf1,xbuf2,nincore,lbuf)

      xnrm=vec_calculate_norm(me_lists, nlists, xbuf1,lbuf)

      contains
      function vec_calculate_norm(mels,nmels, xbuf, lbuf)
      implicit none
      include 'stdunit.h'

      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vec_calculate_norm"
      
      real(8)::
     &     vec_calculate_norm
      
      real(8),intent(inout)::
     &     xbuf(*)

      
      type(me_list_array), intent(in)::
     &     mels(*)

      integer,intent(in)::
     &     nmels,
     &     lbuf                 !> length of buffer
      
      integer::
     &     imel, lenmel, irec
      real(8),external::
     &     ddot
     
      do imel=1,nmels
         if (lbuf.ge. mels(imel)%mel%len_op)then
            lenmel= mels(imel)%mel%len_op
            irec=mels(imel)%mel%fhand%current_record
            call vec_from_da(mels(imel)%mel%fhand,
     &           irec,xbuf,lenmel)
            vec_calculate_norm=vec_calculate_norm+
     &           ddot(lenmel, xbuf,1, xbuf,1)
         else
            call quit(0,i_am, "list to long")
         end if
      end do
      end function
      end subroutine

      
      

      


      
