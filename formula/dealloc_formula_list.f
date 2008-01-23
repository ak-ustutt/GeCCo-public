*----------------------------------------------------------------------*
      subroutine dealloc_formula_list(formula_list)
*----------------------------------------------------------------------*
*     deallocate pointer list (from tail to head, head will not be
*     deallocated due to ifort bug)
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout), target ::
     &     formula_list

      type(formula_item), pointer ::
     &     current
      
      current => formula_list
      ! go to end of list
      do while(associated(current%next))
        current => current%next
      end do
     
      ! loop backward through list and deallocate
      do
        if (associated(current%interm)) deallocate(current%interm)
        if (associated(current%contr)) deallocate(current%contr)
        if (.not.associated(current%prev)) exit
        current => current%prev
        deallocate(current%next)
      end do

      return
      end
