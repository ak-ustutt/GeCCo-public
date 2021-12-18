*----------------------------------------------------------------------*
      subroutine relink_formula_list(form,fpl_reo)
*----------------------------------------------------------------------*
*     re-link the linked list stored on form using the order of the
*     pointer list fpl_reo
*     currently, the first item must not be different
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      type(formula_item), intent(inout), target ::
     &     form
      type(formula_item_list), intent(inout), target ::
     &     fpl_reo

      type(formula_item_list), pointer ::
     &     fpl_pnt
      type(formula_item), pointer ::
     &     form_pnt

      form_pnt => form
      if (.not.associated(form_pnt,fpl_reo%item))
     &     call quit(1,'relink_formula_list',
     &     'first item must remain first')

      fpl_pnt => fpl_reo
      do
        if (associated(fpl_pnt%prev)) then
          fpl_pnt%item%prev => fpl_pnt%prev%item
        else
          fpl_pnt%item%prev => null()
        end if
        if (associated(fpl_pnt%next)) then
          fpl_pnt%item%next => fpl_pnt%next%item
        else
          fpl_pnt%item%next => null()
        end if

        if (.not.associated(fpl_pnt%next)) exit
        fpl_pnt => fpl_pnt%next
      end do

      return
      end
