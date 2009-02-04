*----------------------------------------------------------------------*
      subroutine form_extract_order(f_output,f_input,
     &                      title, label, order, op_info)
*----------------------------------------------------------------------*
*
*     driver for extracting terms of specific order
*     if order<0, terms not obeying (2n+1) and (2n+2) rules are deleted
*
*     f_input and f_output may be identical
*
*     andreas, jan, matthias 2008
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'mdef_formula_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     order
      character(*), intent(in) ::
     &     title, label
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     same, transpose
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item) ::
     &     flist
      integer ::
     &     idx_tgt,iprint

      integer, external ::
     &     idx_oplist2

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'form_extract_order reports')
        write(luout,*) ' f_input = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        write(luout,*) ' order = ',order
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist)

      ! extract terms of order
      idx_tgt = idx_oplist2(label, op_info)
      if (idx_tgt.le.0) call quit(1,'form_extract_order',
     &    'Operator '//label//' not found')
      call pert_order_truncation(flist,order,idx_tgt,op_info)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (iprint.ge.100) then
        call write_title(luout,wst_around_double,'Extracted formula:')
        call print_form_list(luout,flist,op_info)
      end if

      call dealloc_formula_list(flist)
      
      return
      end
