*----------------------------------------------------------------------*
      subroutine clean_formula_dependencies(depend)
*----------------------------------------------------------------------*
*     deallocate subarrays in depend
*----------------------------------------------------------------------*

      implicit none

      include 'def_dependency_info.h'

      type(dependency_info), intent(inout) ::
     &     depend

      if (associated(depend%idxlist))
     &     deallocate(depend%idxlist)
      if (associated(depend%depends_on_idxlist))
     &     deallocate(depend%depends_on_idxlist)

      return
      end
