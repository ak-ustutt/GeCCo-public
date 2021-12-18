*----------------------------------------------------------------------*
      subroutine dealloc_fl_item(fl_item)
*----------------------------------------------------------------------*
*     deallocate all subarrays of fl_item
*     as the fl_item is not passed as pointer (ifort problems)
*     we cannot deallocate it in this routine, it must be done
*     AFTER calling this routine
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout) ::
     &     fl_item

      if (associated(fl_item%contr)) then
        call dealloc_contr(fl_item%contr)
        deallocate(fl_item%contr)
      end if
      if (associated(fl_item%interm)) then
        call dealloc_operator(fl_item%interm)
        deallocate(fl_item%interm,fl_item%parent1,fl_item%parent2)
      end if
      if (associated(fl_item%bcontr)) then
        call dealloc_bcontr(fl_item%bcontr)
        deallocate(fl_item%bcontr)
      end if
      if (associated(fl_item%reo)) then
        call dealloc_reorder(fl_item%reo)
        deallocate(fl_item%reo)
      end if
      if (associated(fl_item%label))
     &  deallocate(fl_item%label)
      
      return
      end
