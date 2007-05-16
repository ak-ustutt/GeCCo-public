*----------------------------------------------------------------------*
      subroutine form_list2arr(form_list,form_arr,nform)
*----------------------------------------------------------------------*
*     set up an array of pointers that point to the entries of the
*     list on form_list
*     new version with a real pointer array
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_formula.h'
      include 'def_formula_list.h'
      include 'def_formula_array.h'

      type(formula_list), intent(in), target ::
     &     form_list
      integer, intent(in) ::
     &     nform
      type(formula_array), intent(out) ::
     &     form_arr(nform)

      type(formula_list), pointer ::
     &     current

      integer ::
     &     iform

      current => form_list

      do iform = 1, nform
        if (.not.associated(current%form))
     &       call quit(1,'form_list2arr',
     &       'unallocated formula on list')
        form_arr(iform)%form => current%form
        if (iform.lt.nform.and..not.(associated(current%next)))
     &       call quit(1,'form_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
