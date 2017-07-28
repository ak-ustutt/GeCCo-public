

      pure function dvdsbsp_get_curlen(dvdsbsp)
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      
      type(davidson_subspace_t),intent(in)::
     &     dvdsbsp
      integer ::
     &     dvdsbsp_get_curlen

      dvdsbsp_get_curlen=dvdsbsp%ncursub



      end function
*----------------------------------------------------------------------*
!>   returns the number of new vectors that do not yet have a matching mvproduct
*----------------------------------------------------------------------*
      pure function dvdsbsp_get_nnew_vvec(dvdsbsp)
      implicit none
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'
      type(davidson_subspace_t),intent(in)::
     &     dvdsbsp
      integer::
     &     dvdsbsp_get_nnew_vvec
      
     
      dvdsbsp_get_nnew_vvec=dvdsbsp%icursub-dvdsbsp%lcursub

      if (dvdsbsp_get_nnew_vvec.lt.0)
     &     dvdsbsp_get_nnew_vvec=dvdsbsp_get_nnew_vvec
     &     +dvdsbsp%ncursub
      end function
