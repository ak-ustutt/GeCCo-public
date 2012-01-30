*----------------------------------------------------------------------*
      subroutine extract_rhs(fl_rhs,fl_traf,fl_raw,
     &     idx_traf,idx_rhs,idx_x,nx,op_info)
*----------------------------------------------------------------------*
*     sort formula into rhs and transformation contributions
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'

      type(formula_item), intent(inout), target ::
     &     fl_rhs, fl_traf, fl_raw
      integer, intent(in) ::
     &     idx_traf, idx_rhs, idx_x(nx), nx
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     ix
      logical ::
     &     found

      type(formula_item), pointer ::
     &     fl_rhs_pnt, fl_traf_pnt, fl_raw_pnt

      integer, external ::
     &     vtx_in_contr

      fl_rhs_pnt => fl_rhs
      fl_traf_pnt => fl_traf
      fl_raw_pnt => fl_raw

      if (fl_raw_pnt%command.ne.command_add_contribution)
     &     fl_raw_pnt => fl_raw_pnt%next

      if (fl_traf_pnt%command.ne.command_end_of_formula)
     &     call quit(1,'extract_rhs','improper positioned fl_traf')
      if (fl_rhs_pnt%command.ne.command_end_of_formula)
     &     call quit(1,'extract_rhs','improper positioned fl_rhs')

      ! set [INIT] items
      call new_formula_item(fl_traf_pnt,
     &         command_set_target_init,idx_traf)
      fl_traf_pnt => fl_traf_pnt%next
      call new_formula_item(fl_rhs_pnt,
     &         command_set_target_init,idx_rhs)
      fl_rhs_pnt => fl_rhs_pnt%next

      do while(fl_raw_pnt%command.ne.command_end_of_formula)
        if (fl_raw_pnt%command.ne.command_add_contribution)
     &       call quit(1,'extract_rhs',
     &       'raw formula must not contain other commands than [ADD]')
        
        found = .true.
        do ix = 1, nx
          if (vtx_in_contr(idx_x(ix),.false.,1,fl_raw_pnt%contr).gt.0)
     &       exit
          found = ix.ne.nx
        end do
        if (found) then
          call new_formula_item(fl_traf_pnt,
     &         command_add_contribution,idx_traf)
          call copy_contr(fl_raw_pnt%contr,fl_traf_pnt%contr)
          ! we assume that the blocks of the raw operator
          ! and that of idx_trad/idx_rhs coincide
          fl_traf_pnt%contr%idx_res = idx_traf
          fl_traf_pnt => fl_traf_pnt%next
        else
          call new_formula_item(fl_rhs_pnt,
     &         command_add_contribution,idx_rhs)
          call copy_contr(fl_raw_pnt%contr,fl_rhs_pnt%contr)
          fl_rhs_pnt%contr%idx_res = idx_rhs
          fl_rhs_pnt => fl_rhs_pnt%next
        end if
        
        fl_raw_pnt => fl_raw_pnt%next
      end do

      return
      end
