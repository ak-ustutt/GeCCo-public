*----------------------------------------------------------------------*
      subroutine form_sum_hermite(f_output,f_input,
     &                      title,label_opres,
     &                      op_info)
*----------------------------------------------------------------------*
*     sum up terms, which are hermitian conjugates
*     driver routine for sum_hermite()
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      character(*), intent(in) ::
     &     title,
     &     label_opres
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     same
      integer ::
     &     idxres
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks form_sum_hermite')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
        write(lulog,*) ' op_res  = ',trim(label_opres)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      if (idxres.lt.0)
     &     call quit(1,'form_sum_hermite',
     &     'required operators are not yet defined? '//
     &       trim(label_opres))

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      call sum_hermite(flist,.true.,.false.,op_info)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
        f_output%comment = trim(title)
      end if
      call write_form_list(f_output%fhand,flist,title)

      call dealloc_formula_list(flist)
      
      return
      end
