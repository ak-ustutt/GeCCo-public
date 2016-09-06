

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
