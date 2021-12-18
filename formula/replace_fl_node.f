*----------------------------------------------------------------------*
      subroutine replace_fl_node(node,nodes_new)
*----------------------------------------------------------------------*
*     remove node from linked list and replace it be new nodes
*     as the node is not passed as pointer (ifort problems)
*     we cannot deallocate it in this routine, it must be done
*     AFTER calling this routine
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout), target ::
     &     node, nodes_new

      type(formula_item), pointer ::
     &     node_pnt

      ! if first node to insert is an [END]: do nothing
      if (nodes_new%command.eq.command_end_of_formula)
     &     return

      node_pnt => nodes_new
      if (associated(node%prev)) node%prev%next => node_pnt
      node_pnt%prev => node%prev
      do while(associated(node_pnt%next))
        node_pnt => node_pnt%next
      end do
      ! make sure that last inserted node is not and [END]
      if (node_pnt%command.eq.command_end_of_formula) then
        node_pnt => node_pnt%prev
        deallocate(node_pnt%next)
      end if
      node_pnt%next => node%next
      if (associated(node%next)) node%next%prev => node_pnt

      if (associated(node%contr)) then
        call dealloc_contr(node%contr)
        deallocate(node%contr)
      end if
      if (associated(node%interm)) then
c        call dealloc_op()
        deallocate(node%interm)
      end if
      
      return
      end
