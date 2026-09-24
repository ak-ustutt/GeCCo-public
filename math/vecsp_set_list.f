
      subroutine vecsp_set_list_mel(vecsp,mel, ivec, ilist, buf, lbuf)
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      integer,parameter ::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="vecsp_set_list"


      type(vector_space_t), intent(inout)::
     &     vecsp
      
      type(me_list),intent(in)::
     &     mel

      integer,intent(in)::
     &     ivec,                !vector the list is to be replaced of
     &     ilist,               !list this is for
     &     lbuf                 !len of buf

      real(8),dimension(lbuf),intent(inout):: 
     &     buf                 
      
      integer::
     &     islice, nslice,
     &     listlen,              !length of list to be copied
     &     irecst

      
      listlen=mel%len_op
      if (listlen .gt. lbuf)
     &     call quit(1,i_am,"not prepared for operator"//
     &     " longer than buffer.")

      call vec_from_da(mel%fhand,mel%fhand%current_record, buf,listlen)
      call vec_to_da(vecsp%vectors(ilist)%fhand, ivec, buf, listlen)
      vecsp%nvec=max(vecsp%nvec,ivec)
      end subroutine 


      subroutine vecsp_set_list_buf(vecsp, buf, ivec, ilist, lbuf)
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      integer,parameter ::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="vecsp_set_list"


      type(vector_space_t), intent(inout)::
     &     vecsp
      

      integer,intent(in)::
     &     ivec,                !vector the list is to be replaced of
     &     ilist,               !list this is for
     &     lbuf                 !len of buf


      real(8),dimension(lbuf),intent(inout):: 
     &     buf                  ! scratch space
      
      integer::
     &     islice, nslice,
     &     listlen,              !length of list to be copied
     &     irecst

      
      listlen=vecsp%me_lists(ilist)%mel%len_op
      if (listlen .gt. lbuf)
     &     call quit(1,i_am,"not prepared for operator"//
     &     " longer than buffer.")

      call vec_to_da(vecsp%vectors(ilist)%fhand, ivec, buf, listlen)
      vecsp%nvec=max(vecsp%nvec,ivec)
      end subroutine 
