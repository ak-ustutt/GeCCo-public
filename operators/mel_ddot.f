
*----------------------------------------------------------------------*
!>      multiplies two me_lists
!
!!      @param[in] me1 -
!!      @param[in] me2 -
!!      @param[inout] buf1,buf2 buffer for the
!!      @param nincore without function
!!      @param lbuf lenght of the buffers
*----------------------------------------------------------------------*
      function me_ddot(me1,me2, buf1, buf2, nincore, lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'ioparam.h'
      
      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="me_ddot"

      real(8)::
     &     me_ddot
      
      type(me_list), intent(in)::
     &     me1,me2       
      real(8),intent(inout)::
     &     buf1(*),buf2(*)


      integer,intent(in)::
     &     lbuf,nincore

      integer::
     &     llist, irec1, irec2

      real(8),external ::
     &     da_ddot
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if

      if (me1%len_op.ne.me2%len_op)
     &     call quit(0,i_am,"list have to be of identical length")
      llist=me1%len_op

      irec1=me1%fhand%current_record
      irec2=me2%fhand%current_record
      me_ddot=da_ddot(me1%fhand, irec1, me2%fhand, irec2,
     &     llist,buf1,buf2,lbuf)
      end function
