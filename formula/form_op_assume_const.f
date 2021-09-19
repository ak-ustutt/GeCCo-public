*----------------------------------------------------------------------*
      subroutine form_op_assume_const(f_output,f_input,
     &                      title,
     &                      nreplace,label_replace,val_replace,
     &                      op_info)
*----------------------------------------------------------------------*
*
*     driver for replacing (scalar) operators in input formula by constant values
*
*     label_replace contains nreplace labels, val_replace the corresp. values
*
*     f_input and f_output may be identical
*
*     andreas, august 2021
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
     &     nreplace
      character(*), intent(in) ::
     &     title,
     &     label_replace(nreplace)
      real(8), intent(in) ::
     &     val_replace(nreplace)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     same
      integer ::
     &     irepl, idx, len
      character ::
     &     name*(form_maxlen_label*2)

      type(filinf), pointer ::
     &     ffintm
      type(formula_item) ::
     &     flist, fl_intm

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &     'form_op_assume_const reports')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
        do irepl = 1, nreplace
          write(lulog,*) irepl,trim(label_replace(irepl)),' ',
     &         val_replace(irepl)
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! loop over intermediates
      do irepl = 1, nreplace

        ! get index
        idx = idx_oplist2(trim(label_replace(irepl)),op_info)
        if (idx.le.0)
     &       call quit(1,'form_op_assume_const',
     &    'operator label 1 not on list: '//trim(label_replace(irepl)))
        
        
        if (ntest.ge.100)
     &       write(lulog,*)
     &       'now replacing: ',
     &       trim(label_replace(irepl)),
     &       ' by ',
     &       val_replace(irepl)

        call op_assume_const(idx,val_replace(irepl),flist,op_info)

      end do

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (ntest.ge.100) then
        call write_title(lulog,wst_around_double,'Modified formula:')
        call print_form_list(lulog,flist,op_info)
      end if

      call dealloc_formula_list(flist)

      return
      end
