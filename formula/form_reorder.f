*----------------------------------------------------------------------*
      subroutine form_reorder(f_output,f_input,
     &                      title,op_info)
*----------------------------------------------------------------------*
*
*     driver for call to reorder_formula (command REORDER_FORMULA)
*
*     andreas, october 2012
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
c      include 'def_contraction_list.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      character(*), intent(in) ::
     &     title
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     ok, same, transpose
      integer ::
     &     nterms, idum, idxinp, idx, icmpnd, idxres, len
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'here speaks form_indep')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_reorder',
     &     'empty formula list? something is buggy')

      call reorder_formula(flist,op_info)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      call dealloc_formula_list(flist)

      return
      end
