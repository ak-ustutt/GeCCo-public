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

      integer ::
     &     isave
      
      current => formula_list
      ! go to end of list
      isave = 0
      do while(associated(current%next))
        current => current%next
        isave = isave+1
        if (isave.gt.10 000 000)
     &       call quit(1,'dealloc_formula_list','infinite loop (1)?')
      end do
     
      isave = 0
      ! loop backward through list and deallocate
      do
        call dealloc_fl_item(current)
        if (.not.associated(current%prev)) exit
        current => current%prev
        deallocate(current%next)
        isave = isave+1
        if (isave.gt.10 000 000)
     &       call quit(1,'dealloc_formula_list','infinite loop (2)?')
      end do

      return
      end
