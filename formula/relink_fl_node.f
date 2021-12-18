*----------------------------------------------------------------------*
      subroutine relink_fl_node(node_move,node_append)
*----------------------------------------------------------------------*
*     relink node in linked list
*     node_move is removed from prev. position and inserted
*     after node_append
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout), target ::
     &     node_move,
     &     node_append

      ! link the nodes around node_move
      if (associated(node_move%prev)) 
     &               node_move%prev%next => node_move%next
      if (associated(node_move%next)) 
     &               node_move%next%prev => node_move%prev

      ! and relink
      node_move%prev => node_append
      node_move%next => node_append%next
      node_append%next => node_move
      if (associated(node_move%next))
     &               node_move%next%prev => node_move

      
      return
      end
