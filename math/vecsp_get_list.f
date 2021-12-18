
*----------------------------------------------------------------------*
!>    retrieves the ilist list from ivec vector  from a vectorspace onto a buffer
!!
!!    @param[in] vecsp the vectorspace
!!    @param[out] listlen  length of the list
!!    @param[inout] buf contains the elements of the list on exit
!!    @param[in] lbuf absolute length of the buffer
*----------------------------------------------------------------------*
      subroutine vecsp_get_list_buf(vecsp, ivec, ilist, listlen,
     &     buf, lbuf)
*----------------------------------------------------------------------*

      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer,parameter ::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="vecsp_get_list_buf"


      type(vector_space_t), intent(inout)::
     &     vecsp

      integer,intent(in)::
     &     ivec,                !vector the list is to be replaced of
     &     ilist,               !list this is for
     &     lbuf                 !len of buf
      integer,intent(out)::
     &     listlen

      real(8),dimension(lbuf),intent(inout)::
     &     buf

      listlen=vecsp%me_lists(ilist)%mel%len_op
      if ( listlen.gt.lbuf)
     &     call quit(1,i_am,"buffer too short")

      call vec_from_da(vecsp%vectors(ilist)%fhand, ivec, buf,
     &     listlen)

      return
      end subroutine
