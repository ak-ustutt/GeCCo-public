*----------------------------------------------------------------------*
      subroutine delete_fl_node(node)
*----------------------------------------------------------------------*
*     remove node from linked list
*     as the node is not passed as pointer (ifort problems)
*     we cannot deallocate it in this routine, it must be done
*     AFTER calling this routine
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout) ::
     &     node

      if (associated(node%prev)) node%prev%next => node%next
      if (associated(node%contr)) then
        call dealloc_contr(node%contr)
        deallocate(node%contr)
      end if
      if (associated(node%interm)) then
c        call dealloc_op()
        deallocate(node%interm)
      end if
      if (associated(node%next)) node%next%prev => node%prev
      
      return
      end
