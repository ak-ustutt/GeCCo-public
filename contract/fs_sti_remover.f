*----------------------------------------------------------------------*
      subroutine fs_sti_remover(fl_item,op_info)
*----------------------------------------------------------------------*
*     remove operators and list that are marked as 
*     short term intermediates (names start with _STIN)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(in) ::
     &     fl_item

      type(operator_info), intent(inout) ::
     &     op_info

      logical ::
     &     is_bc, is_reo
      integer ::
     &     n_operands, iop, idx
      character(maxlen_bc_label) ::
     &     label_op
      character(mxlen_melabel) ::
     &     label_mel
      type(binary_contr), pointer ::
     &     bc
      type(reorder), pointer ::
     &     reo

      integer, external ::
     &     idx_oplist2

      select case(fl_item%command)
      case(command_add_intm,command_cp_intm,command_bc,command_add_bc,
     &       command_bc_reo,command_add_bc_reo)
        if (.not.associated(fl_item%bcontr))
     &     call quit(1,'fs_sti_remover','no bcontr?')
        is_bc  = .true.
        is_reo = .false.
      case(command_reorder,command_add_reo)
        if (.not.associated(fl_item%reo))
     &     call quit(1,'fs_sti_remover','no reo?')
        is_bc  = .false.
        is_reo = .true.
      case default
        return
      end select

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'fs_sti_remover active')
        write(lulog,*) 'is_bc, is_reo: ',is_bc,is_reo
      end if

      if (is_bc) then
        bc => fl_item%bcontr
        n_operands = bc%n_operands
      else if (is_reo) then
        reo => fl_item%reo
        n_operands = 1
      else
        call quit(1,'fs_sti_remover','unknown route')
      end if

      ! check whether any of the input operators is a 
      ! short-term intermediate (STIN)
      
      do iop = 1, n_operands

        if (is_reo.and.iop.eq.1) label_op = reo%label_in
        if (is_bc.and.iop.eq.1) label_op = bc%label_op1
        if (is_bc.and.iop.eq.2) label_op = bc%label_op2

        if (ntest.ge.100)
     &       write(lulog,*) 'current op label: ',trim(label_op)

        if (label_op(1:5).ne.'_STIN') cycle

        idx = idx_oplist2(label_op,op_info)
        label_mel = op_info%op_arr(idx)%op%assoc_list

        if (ntest.ge.100)
     &       write(lulog,*) 'removing: ',trim(label_op),trim(label_mel)
        
        call del_me_list(label_mel,op_info)
        call del_operator(label_op,op_info)

      end do

      return
      end
