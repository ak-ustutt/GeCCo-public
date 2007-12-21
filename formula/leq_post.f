*----------------------------------------------------------------------*
      subroutine leq_post(label_f_traf,label_f_rhs,label_f_raw,
     &                    label_traf,label_rhs,label_x,
     &                    title_traf,title_rhs,
     &                    op_info,form_info)
*----------------------------------------------------------------------*
*     some post-processing for linear equations
*     split formula labelled "label_f_raw" into transformation
*     ("label_f_traf") and right-hand side part ("label_f_rhs")
*     reorder formulae
*     label_traf :  operator associated with Ax-trafo
*     label_rhs  :  operator associated with rhs
*     label_x    :  operator associated with variable x
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'mdef_formula_info.h'

      integer, parameter ::
     &     ntest = 100

      character(*), intent(in) ::
     &     label_f_raw, label_f_traf, label_f_rhs,
     &     label_traf, label_rhs, label_x
      type(operator_info), intent(in) ::
     &     op_info
      type(formula_info), intent(in) ::
     &     form_info
      character, intent(in) ::
     &     title_traf*(*), title_rhs*(*)

      integer ::
     &     idx, idx_traf, idx_rhs, idx_x
      type(formula), pointer ::
     &     form_traf, form_rhs, form_raw
      type(formula_item) ::
     &     fl_traf, fl_rhs, fl_raw

      integer, external ::
     &     idx_formlist, idx_oplist2

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'leq_post at work')

      ! check labels
      ! this one must exist:
      idx = idx_formlist(label_f_raw,form_info)
      if (idx.le.0)
     &     call quit(1,'leq_post',
     &     'label not on list '//trim(label_f_raw))
      form_raw => form_info%form_arr(idx)%form

      idx = idx_formlist(trim(label_f_rhs),form_info)
      if (idx.le.0) then
        call add_formula(form_info,trim(label_f_rhs))
        idx = idx_formlist(trim(label_f_rhs),form_info)
      end if
      form_rhs => form_info%form_arr(idx)%form

      idx = idx_formlist(trim(label_f_traf),form_info)
      if (idx.le.0) then
        call add_formula(form_info,trim(label_f_traf))
        idx = idx_formlist(trim(label_f_traf),form_info)
      end if
      form_traf => form_info%form_arr(idx)%form

      ! operator labels (must exist)
      idx_traf = idx_oplist2(label_traf,op_info)
      idx_rhs  = idx_oplist2(label_rhs,op_info)
      idx_x    = idx_oplist2(label_x,op_info)

      if (idx_traf.le.0.or.idx_rhs.le.0.or.idx_x.le.0) then
        write(luout,*) 'labels: "',trim(label_traf),'" "',
     &                             trim(label_rhs), '" "',
     &                             trim(label_x),   '"'
        write(luout,*) 'indices: ',idx_traf,idx_rhs,idx_x
        call quit(1,'leq_post',
     &       'one or more of the above labels was not defined')
      end if

      ! initialize lists:
      call init_formula(fl_raw)
      call init_formula(fl_rhs)
      call init_formula(fl_traf)

      ! read raw list
      call read_form_list(form_raw%fhand,fl_raw)

      ! sort into traf and rhs list
      call extract_rhs(fl_rhs,fl_traf,fl_raw,
     &       idx_traf,idx_rhs,idx_x,op_info)

      ! reorder new lists
      call reorder_formula(fl_rhs,op_info)
      call reorder_formula(fl_traf,op_info)

      ! write lists to file
      ! check whether raw formula was identical to either
      ! of the new formulae
      form_traf%comment = title_traf
      call write_form_list(form_traf%fhand,fl_traf,title_traf)

      form_rhs%comment = title_rhs
      call write_form_list(form_rhs%fhand,fl_rhs,title_rhs)

      if (ntest.ge.100) then
        call write_title(luout,wst_around_single,'Extracted lists')
        call write_title(luout,wst_uline_single,'RHS')
        call print_form_list(luout,fl_rhs,op_info)
        call write_title(luout,wst_uline_single,'Trafo')
        call print_form_list(luout,fl_traf,op_info)
      end if

      ! tidy up again
      call dealloc_formula_list(fl_raw)
      call dealloc_formula_list(fl_rhs)
      call dealloc_formula_list(fl_traf)

      return
      end
