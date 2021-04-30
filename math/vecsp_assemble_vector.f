*----------------------------------------------------------------------*
!>     assembles a vector from a vectorspace
!!     @param vecsp the vectorspace
!!     @param me_lists array of me_lists that make up the output vector
!!     @param nmels number of me-lists in me_lists
!!     @param coeff vector with coefficients for the linear combination
!!     @param lcoeff lengh of the coefficent vector
!!     @param xfac the elements of the initial me-lists contribute scaled by xfac to the final vector
!!     @param xbuf1,xbuf2 scratch space status on input and output undefined
!!     @param nincore how many scratch buffer are usable
!!     @param lbuf length of the scratch buffer
*----------------------------------------------------------------------*
      subroutine vecsp_assemble_vector(vecsp,
     &     me_lists, nlists,
     &     coeff, lcoeff,
     &     xfac,
     &     xbuf1, xbuf2, nincore, lbuf)
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
     &     i_am="vecsp_assemble_vector"



      integer,intent(inout)::
     &     nlists,               !> number of me-lists in a vector
     &     lcoeff,                !> length of coeff
     &     nincore, lbuf        !> number of buffers and length of buffers

      type(vector_space_t),intent(inout)::
     &     vecsp                !>

      type(me_list_array), intent(in)::
     &     me_lists(*)          !> me-lists the vector is assembled at

      real(8),intent(in)::
     &     xfac,                !> factor the values from the existing me-lists contribute to the vector
     &     coeff(*)               !> coefficients for the vectors of the vectorspace

      real(8),intent(inout)::
     &     xbuf1(*),            !> accumulator
     &     xbuf2(*)             !>scratch memory


      type(filinf),pointer::
     &     ffvec

      logical::
     &     file_close

      integer::
     &     ilist, icoeff, ii,
     &     lenlist

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      if (lcoeff .gt. vecsp%nvec) call quit(0,i_am,
     &     "more coefficients than vectors")

      if (nincore.ge.2) then
         do ilist=1,nlists
            if (.not. mel_equal_format(me_lists(ilist)%mel,
     &           vecsp%me_lists(ilist)%mel))
     &           call quit(1,i_am,"inconsistent format of me-lists")

            lenlist=me_lists(ilist)%mel%len_op

            if (lenlist.gt.lbuf)
     &           call quit(0,i_am,"not prepared for long lists")

            ffvec=> me_lists(ilist)%mel%fhand
            if (ffvec%unit.lt.0)then
               call file_open(ffvec)
               file_close=.True.
            else
               file_close=.False.
            end if

            call prepare_accumulator(xbuf1, lenlist, xfac, ffvec)

            do icoeff=1,lcoeff
               if (coeff(icoeff).eq.0d0) cycle
               call vecsp_get_list_buf(vecsp, icoeff, ilist, lenlist,
     &              xbuf2, lbuf)
               xbuf1(1:lenlist)=xbuf1(1:lenlist)+
     &              xbuf2(1:lenlist)*coeff(icoeff)
            end do
            call put_accumulator(xbuf1,lenlist,ffvec)

            if (file_close) call file_close_keep(ffvec)
         end do
      else
         call quit(1,i_am, "not prepared for only 1 incore buffer")

      end if

      return
      contains

*----------------------------------------------------------------------*
!>    tests if two me-lists have compatible formats
!!    @param me_list1,me_list2 the me-lists to be compared
!!    @TODO implement fully
*----------------------------------------------------------------------*
      pure function mel_equal_format(me_list1, me_list2)
*----------------------------------------------------------------------*
      logical::
     &     mel_equal_format

      type(me_list),intent(in)::
     &     me_list1,me_list2

      mel_equal_format=me_list1%len_op .eq.me_list2%len_op !> @TODO implement more tests

      end function
*----------------------------------------------------------------------*
!>    prepares the accumulator buffer
!!    @param xbuf buffer for the accumulator
!!    @param llist length of the used part of the accumulator(must be smaller than length of buffer)
!!    @param xfac the accumulator is initalized with the element of the current list from ffhand scaled by xfac
!!    @param ffhand filinf object for the initializaion of the buffer (current_record is used)
*----------------------------------------------------------------------*
      subroutine prepare_accumulator(xbuf, llist, xfac, ffhand)
*----------------------------------------------------------------------*
      real(8),intent(inout)::
     &     xbuf(*)
      real(8),intent(in)::
     &     xfac

      integer,intent(in)::
     &     llist

      type(filinf),intent(inout)::
     &     ffhand
      integer::
     &     irec


      if(xfac.eq.0d0) then
         xbuf1(1:llist) = 0d0
      else
         irec=ffvec%current_record
         call vec_from_da(ffvec,irec,xbuf1,llist)

         if (xfac.eq.1d0)then
            continue
         else
            xbuf1(1:llist)=xbuf1(1:llist)*xfac
         end if
      end if
      return
      end subroutine

*----------------------------------------------------------------------*
!>    saves the accumulator buffer to file
!!    @param xbuf1 accumulator buffer
!!    @param llist lenght of the used part of the accumulator buffer
!!    @param ffhand filinf for the file the buffer is saved to
*----------------------------------------------------------------------*
      subroutine put_accumulator(xbuf1,llist,ffhand)
*----------------------------------------------------------------------*
      implicit none
      real(8),intent(inout)::
     &     xbuf1(*)
      integer,intent(in)::
     &     llist

      type(filinf),intent(inout)::
     &     ffhand

      integer::
     &     irec

      irec=ffhand%current_record
      call vec_to_da(ffvec,irec,xbuf1,llist)

      end subroutine
      end subroutine
