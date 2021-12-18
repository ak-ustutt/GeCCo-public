*----------------------------------------------------------------------*
      subroutine init_formula_plist(form_plist)
*----------------------------------------------------------------------*
*     initialize all pointers in first formula item pointer
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      type(formula_item_list), intent(inout) ::
     &     form_plist

      form_plist%prev => null()
      form_plist%next => null()
      form_plist%item => null()

      return
      end
