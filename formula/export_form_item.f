*----------------------------------------------------------------------*
      subroutine export_form_item(lulog,idx,fl_item,op_info)
*----------------------------------------------------------------------*
*     export formula item to unit lulog
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lulog
      integer, intent(inout) ::
     &     idx
      type(formula_item), intent(inout), target ::
     &     fl_item
      type(operator_info), intent(inout) ::
     &     op_info


      select case(fl_item%command)
      case(command_end_of_formula)
        write(lulog,'(a)') '[END]'
      case(command_set_target_init)
        !write(lulog,'(a)') '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
        !write(lulog,'(a)') '[SET TARGET]',fl_item%target
      case(command_add_contribution)
        idx = idx+1
        write(lulog,'(a,i8)') '[CONTR] # ',idx
        call contr_set_indices(fl_item%contr,op_info)
        call contr_get_eqv_line_factor(fl_item%contr,op_info)
        call prt_contr_export(lulog,fl_item%contr,op_info)
      case(command_symmetrise)
        write(lulog,'(a)') '[SYMMETRISE]',fl_item%target
      case default
        write(lulog,*) 'unknown command ',fl_item%command,
     &       fl_item%target
      end select
      
      end
