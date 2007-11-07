*----------------------------------------------------------------------*
      subroutine leq_post(form_traf,form_rhs,form_raw,
     &                    split,
     &                    idx_traf,idx_rhs,idx_x,
     &                    title_traf,title_rhs,
     &                    op_info)
*----------------------------------------------------------------------*
*     some post-processing for linear equations
*     split into transformation and right-hand side part
*          (if split==.true.)
*     reorder formulae
*     idx_traf :  index of operator associated with Ax-trafo
*     idx_rhs  :  index of operator associated with rhs
*     idx_x    :  index of operator associated with variable x
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      type(formula), intent(inout), target ::
     &     form_traf, form_rhs, form_raw
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     split
      integer, intent(in) ::
     &     idx_traf, idx_rhs, idx_x
      character, intent(in) ::
     &     title_traf*(*), title_rhs*(*)

      type(formula_item) ::
     &     fl_traf, fl_rhs, fl_raw
      character ::
     &     name*256


      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'leq_post at work')

      ! initialize lists:
      if (split) call init_formula(fl_raw)
      if (split) call init_formula(fl_rhs)
      call init_formula(fl_traf)

      ! read raw list
      if (split) then
        call read_form_list(form_raw%fhand,fl_raw)

        ! sort into traf and rhs list
        call extract_rhs(fl_rhs,fl_traf,fl_raw,
     &       idx_traf,idx_rhs,idx_x,op_info)

      else
        call read_form_list(form_raw%fhand,fl_traf)
      end if
      
      ! reorder new lists
      if (split) call reorder_formula(fl_rhs,op_info)
      call reorder_formula(fl_traf,op_info)

      ! write lists to file
      ! check whether raw formula was identical to either
      ! of the new formulae
      if (form_raw%label.ne.form_traf%label) then
        write(name,'(a,".fml")') trim(form_traf%label)
        call file_init(form_traf%fhand,name,ftyp_sq_unf,0)
      end if
      form_traf%comment = title_traf
      call write_form_list(form_traf%fhand,fl_traf,title_traf)

      if (split) then
        if (form_raw%label.ne.form_rhs%label) then
          write(name,'(a,".fml")') trim(form_rhs%label)
          call file_init(form_rhs%fhand,name,ftyp_sq_unf,0)
        end if
        form_rhs%comment = title_rhs
        call write_form_list(form_rhs%fhand,fl_rhs,title_rhs)
      end if

      if (split.and.ntest.ge.100) then
        call write_title(luout,wst_around_single,'Extracted lists')
        call write_title(luout,wst_uline_single,'RHS')
        call print_form_list(luout,fl_rhs,op_info)
        call write_title(luout,wst_uline_single,'Trafo')
        call print_form_list(luout,fl_traf,op_info)
      end if

      ! tidy up again
      if (split) call dealloc_formula_list(fl_raw)
      if (split) call dealloc_formula_list(fl_rhs)
      call dealloc_formula_list(fl_traf)

      return
      end
