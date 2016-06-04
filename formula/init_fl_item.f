*----------------------------------------------------------------------*
      subroutine init_fl_item(form_pnt,command,target)
*----------------------------------------------------------------------*
*     prepare new entry to take the information coming
*     along with command "command" (i.e. contr or interm is alloc.d)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      
      type(formula_item), intent(inout), target ::
     &     form_pnt
      integer, intent(in) ::
     &     command, target

      select case(command)
      case(command_set_target_init,
     &     command_set_target_update,
     &     command_symmetrise,
     &     command_internal)
        ! nothing to do
      case(command_add_contribution)
        allocate(form_pnt%contr)
        call init_contr(form_pnt%contr)
      case(command_del_intermediate)
        allocate(form_pnt%label)
      case(command_new_intermediate)
        allocate(form_pnt%interm,form_pnt%parent1,form_pnt%parent2)
        call init_operator_0(form_pnt%interm)
        form_pnt%parent1(1:len_opname) = ' '
        form_pnt%parent2(1:len_opname) = ' '
        form_pnt%incore = 0
        form_pnt%tra = .false.
        form_pnt%tra1 = .false.
        form_pnt%tra2 = .false.
      case(command_reorder)
        allocate(form_pnt%reo)
        call init_reorder(form_pnt%reo)
      case(command_add_bc_reo,command_bc_reo,command_add_reo)
        allocate(form_pnt%bcontr)
        call init_bcontr(form_pnt%bcontr)
        allocate(form_pnt%reo)
        call init_reorder(form_pnt%reo)
      case(command_add_intm,command_cp_intm,command_add_bc,command_bc)
        allocate(form_pnt%bcontr)
        call init_bcontr(form_pnt%bcontr)
      case default
        call quit(1,'init_fl_item','unknown mode')
      end select

      form_pnt%command = command
      form_pnt%target = target
      
      return
      end
