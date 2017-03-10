*----------------------------------------------------------------------*
      subroutine transpose_formula(form_head,op_info,multi)
*----------------------------------------------------------------------*
*     transpose all entries on a formula list
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      !include 'ifc_formula.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in), optional ::
     &     multi

      type(formula_item), pointer ::
     &     form_ptr
      type(contraction), pointer ::
     &     contr
      logical ::
     &     scal, use_multi

      if(present(multi)) then
        use_multi = multi
      else
        use_multi = .false. 
      endif

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'transpose_formula')
        write(lulog,*) 'formula on entry'
        call print_form_list(lulog,form_head,op_info)
      end if

      ! A new logical 'scal' has been used as the change of sign of
      ! of the target operator is not necessary while transposing 
      ! the formula for operator corresponding to CI coefficients
      if(use_multi) then 
        scal = op_info%op_arr(form_head%target)%op%n_occ_cls.gt.1
      else 
        scal = .true.
      end if

      form_ptr => form_head
      do
        
        if (scal) form_ptr%target = -form_ptr%target

        if (form_ptr%command.eq.command_add_contribution)
     &       call transpose_contr(form_ptr%contr,op_info, use_multi)

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      if (ntest.ge.100) then
        write(lulog,*) 'formula on exit'
        call print_form_list(lulog,form_head,op_info)
      end if

      return
      end

