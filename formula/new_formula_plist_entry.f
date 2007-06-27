*----------------------------------------------------------------------*
      subroutine new_formula_plist_entry(fpl_pnt)
*----------------------------------------------------------------------*
*     allocate new entry for linked list
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'
      
      type(formula_item_list), intent(inout), target ::
     &     fpl_pnt

      allocate(fpl_pnt%next)
      fpl_pnt%next%prev => fpl_pnt
      fpl_pnt%next%next => null()
      
      return
      end
