*----------------------------------------------------------------------*
      subroutine set_primitive_formula(form,idx_op,op_info)
*----------------------------------------------------------------------*
*     store blocks of operator in formula format, corresponds to
*       F = Op(block1) + Op(block2) + ....
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout), target ::
     &     form
      integer, intent(in) ::
     &     idx_op
      type(operator_info), intent(in), target ::
     &     op_info

      type(operator), pointer ::
     &     op
      type(formula_item), pointer ::
     &     form_pnt
      integer ::
     &     iblk_op
      
      op => op_info%op_arr(idx_op)%op
      form_pnt => form

      if (form_pnt%command.ne.command_end_of_formula)
     &     call quit(1,'set_primitive_formula',
     &     'formula seems positioned incorrectly on input')
      if (associated(form_pnt%contr))
     &     call quit(1,'set_primitive_formula',
     &     'contr is already associated on entry')

      do iblk_op = 1, op%n_occ_cls
        ! new entry
        call new_formula_item(form_pnt,command_add_contribution,idx_op)
        
        ! set "contraction"
        call resize_contr(form_pnt%contr,1,0,0,0)
        form_pnt%contr%idx_res = idx_op
        form_pnt%contr%iblk_res = iblk_op
        form_pnt%contr%nvtx = 1
        form_pnt%contr%vertex(1)%idx_op = idx_op
        form_pnt%contr%vertex(1)%iblk_op = iblk_op
        
        form_pnt => form_pnt%next
      end do

      return
      end
