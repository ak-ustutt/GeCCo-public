*----------------------------------------------------------------------*
      subroutine dealloc_formula_plist(formula_plist)
*----------------------------------------------------------------------*
*     deallocate pointer list (from tail to head, head will not be
*     deallocated due to ifort bug)
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'

      type(formula_item_list), intent(inout), target ::
     &     formula_plist

      type(formula_item_list), pointer ::
     &     current
      
      current => formula_plist
      ! go to end of list
      do while(associated(current%next))
        current => current%next
      end do
     
      ! loop backward through list and deallocate
      do
        if (.not.associated(current%prev)) exit
        current => current%prev
        deallocate(current%next)
      end do

      return
      end
