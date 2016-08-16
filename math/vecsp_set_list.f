
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
     &     buf                  ! scratch space
      
      integer::
     &     islice, nslice,
     &     listlen,              !length of list to be copied
     &     irecst
      type(filinf)::
     &     mel_fhand, vec_fhand

      
      listlen=mel%len_op
      if (lenlist .gt. lbuf)
     &     call quit(1,i_am,"not prepared for operator"//
     &     " longer than buffer.")

      mel_fhand=mel%fhand
      vec_fhand=vecsp%vectors(ilist)%fhand
      
      call vec_from_da(mel_fhand, mel_fhand%current_record, buf,listlen)
      call vec_to_da(vec_fhand, ivec, buf, listlen)

      end subroutine 

