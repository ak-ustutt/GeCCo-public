*----------------------------------------------------------------------*
      subroutine transpose_formula_wrap(label_in,label_res,
     &      label_op,init,multi,
     &      op_info,form_info)
*----------------------------------------------------------------------*
*     wrapper for transpose_formula
*----------------------------------------------------------------------*

      implicit none 

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'mdef_formula_info.h'
      include 'def_formula_item.h'
!     include 'def_formula.h'
      include 'ifc_formula.h'
      
      character(*), intent(in) ::
     &     label_in, label_res, label_op
      logical, intent(in) ::
     &     init, multi
      type(operator_info), intent(inout) ::
     &     op_info      
      type(formula_info), intent(inout) ::
     &     form_info

      integer ::
     &     idx_in,idx_res,idx_op
      type(formula), pointer ::
     &     form_in, form_res      
      type(formula_item), target ::
     &     flist_in
      type(formula_item), pointer ::
     &     flist_res

      integer, external ::
     &     idx_formlist, idx_oplist2

      idx_op = idx_oplist2(trim(label_op),op_info)

      if (idx_op.le.0)
     &       call quit(0,'transpose_formula_wrap',
     &       'operator label does not exist: '//trim(label_op))

      idx_in = idx_formlist(trim(label_in),form_info)

      if (idx_in.le.0)
     &       call quit(0,'transpose_formula_wrap',
     &       'input formula does not exist: '//trim(label_in))

      form_in => form_info%form_arr(idx_in)%form

      idx_res = idx_formlist(trim(label_res),form_info)

      if (init) then
        if (idx_res.gt.0)
     &       call quit(0,'transpose_formula_wrap',
     &       'formula does already exist: '//trim(label_res))
        call add_formula(form_info,trim(label_res))
        idx_res = idx_formlist(trim(label_res),form_info)
      else 
        if (idx_res.le.0)
     &       call quit(0,'transpose_formula_wrap',
     &       'formula does not exist: '//trim(label_res))
      end if 
 
      form_res => form_info%form_arr(idx_res)%form

      if (init)
     &  call new_formula_item(flist_res,
     &       command_set_target_init,idx_op)

      call init_formula(flist_in)
      call read_form_list(form_in%fhand,flist_in,.true.)

      call transpose_formula(flist_in,
     &       op_info, multi)

!     flist_in%target = - flist_in%target

!     call copy_fml(flist_in,flist_res,idx_op)

      flist_res => flist_in
!
!     flist_res%target = idx_op
!     flist_res%contr%idx_res = idx_op

!     fl_loop: do

!       if (flist_res%command.eq.command_end_of_formula)
!    &       exit fl_loop


!       flist_res%target = idx_op

!       if (associated(flist_res%contr)) then

!       flist_res%contr%idx_res = idx_op

!      endif

!       flist_res => flist_res%next

!       if (.not.associated(flist_res))
!    &       call quit(1,'copy_fml',
!    &       'unexpected end of formula list')

!     end do fl_loop
 
      call write_form_list(form_res%fhand,flist_in,'XXX')

!     call read_form_list(form_res%fhand,flist,.true.)

!     flist%target = idx_op

!     call write_form_list(form_res%fhand,flist,'XXX')

      return
      end
