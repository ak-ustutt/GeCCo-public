*----------------------------------------------------------------------*
      subroutine print_form_count(lulog,label,form_head,op_info)
*----------------------------------------------------------------------*
*     print info about formula on linked list to unit lulog
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog
      character(len=form_maxlen_label) ::
     &     label
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx_contr, idx_bcontr

      type(formula_item), pointer ::
     &     form_ptr

      ! only one of the two should become non-zero:
      idx_contr = 0
      idx_bcontr = 0
      form_ptr => form_head
      do
        select case(form_ptr%command)
        case(command_add_contribution)
          idx_contr = idx_contr+1
        case(command_add_intm,command_cp_intm,command_add_bc,
     &       command_add_bc_reo,command_bc,command_bc_reo,
     &       command_add_reo)
          idx_bcontr = idx_bcontr+1
        end select
        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      if (idx_contr>0) then
         write(lulog,'(1x,a,": ",i8," diagrams")') trim(label),idx_contr
      end if
      if (idx_bcontr>0) then
         write(lulog,'(1x,a,": ",i8," binary tensor contractions")')
     &     trim(label), idx_bcontr
      end if

      return
      end

