*----------------------------------------------------------------------*
      subroutine print_form_item(luout,idx,fl_item,op_info)
*----------------------------------------------------------------------*
*     print formula item to unit luout
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     luout
      integer, intent(inout) ::
     &     idx
      type(formula_item), intent(in), target ::
     &     fl_item
      type(operator_info), intent(in) ::
     &     op_info

      select case(fl_item%command)
      case(command_end_of_formula)
        write(luout,*) '[END]'
      case(command_set_target_init)
        write(luout,*) '[INIT TARGET]',fl_item%target
      case(command_set_target_update)
        write(luout,*) '[SET TARGET]',fl_item%target
      case(command_new_intermediate)
        write(luout,*) '[NEW INTERMEDIATE]',fl_item%target
        write(luout,'(2x,a)') trim(fl_item%interm%name)
        write(luout,'(2x,"attribute parentage: ",a," ",a)')
     &                        trim(fl_item%parent1),
     &                        trim(fl_item%parent2)
        write(luout,'(2x,"incore: ",i2)') fl_item%incore
        call print_op_occ(luout,fl_item%interm)
      case(command_del_intermediate)
        write(luout,*) '[DELETE INTERMEDIATE]',fl_item%target
        write(luout,'(2x,a)') trim(fl_item%label)
      case(command_add_contribution)
        idx = idx+1
        write(luout,*) '[CONTR]',fl_item%target,'( term #',idx,')'
        call prt_contr2(luout,fl_item%contr,op_info)
      case(command_add_intm)
        idx = idx+1
        write(luout,*) '[ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
      case(command_cp_intm)
        idx = idx+1
        write(luout,*) '[COPY]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
      case(command_add_bc)
        idx = idx+1
        write(luout,*) '[CONTRACT][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
      case(command_add_bc_reo)
        idx = idx+1
        write(luout,*) '[CONTRACT][REORDER][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
        call prt_reorder(luout,fl_item%reo)
      case(command_bc)
        idx = idx+1
        write(luout,*) '[CONTRACT]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
      case(command_bc_reo)
        idx = idx+1
        write(luout,*) '[CONTRACT][REORDER]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
        call prt_reorder(luout,fl_item%reo)
      case(command_reorder)
        write(luout,*) '[REORDER]',fl_item%target
        call prt_reorder(luout,fl_item%reo)
      case(command_add_reo)
        idx = idx+1
        write(luout,*) '[REORDER][ADD]',
     &       fl_item%target,'( term #',idx,')'
        call prt_bcontr(luout,fl_item%bcontr)
        call prt_reorder(luout,fl_item%reo)
      case(command_symmetrise)
        write(luout,*) '[SYMMETRISE]',fl_item%target
      case default
        write(luout,*) 'unknown command ',fl_item%command,
     &       fl_item%target
      end select
      
      end
