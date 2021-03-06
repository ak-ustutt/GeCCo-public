*----------------------------------------------------------------------*
      subroutine form_select_special(f_output,f_input,
     &                      labels,nlabels,
     &                      type,mode,
     &                      op_info)
*----------------------------------------------------------------------*
*     driver routine for various special purpose selection routines
*     switch to routine by type; mode can contain additional information
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

      integer ::
     &     nlabels
      character(len=*) ::
     &     type, mode, labels(nlabels)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     on_list, same, delete
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next


      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks form_select_special')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      if (ntest.ge.150) then
        write(lulog,*) 'initial formula:'
        call print_form_list(lulog,flist,op_info)
      end if

      select case(trim(type))
      case('F12x','f12x')
        call select_f12x(flist,labels,nlabels,mode,op_info)
      case('MRCC_F12','mrcc_f12')
        call select_mrcc_f12(flist,labels,nlabels,mode,op_info)
      case('MRCC3','mrcc3')
        call select_mrcc_lag3(flist,labels,nlabels,mode,op_info)
      case('MRCC2','mrcc2')
        call select_mrcc_lag2(flist,labels,nlabels,mode,op_info)
      case('MRCC','mrcc')
        call select_mrcc_lag(flist,labels,nlabels,mode,op_info)
      case('MRCCtrunc','mrcctrunc')
        call select_mrcc_trunc(flist,labels,nlabels,mode,op_info)
      case('MRCCrem0res','mrccrem0res')
        call select_mrcc_rem0res(flist,labels,nlabels,mode,op_info)
      case('CONNECTED','connected')
          call select_connected(flist,labels,nlabels,mode,op_info)
      case('OPT1','opt1')
        call select_xsp_opt1(flist,labels,nlabels,mode,op_info)
      case('SAME','same')
        call select_same_blk(flist,labels,nlabels,mode,op_info)
      case('HTT','htt')
        call select_repl_htt(flist,labels,nlabels,mode,op_info)
      case('FORMAL','formal')
        call select_formal_blk(flist,mode,op_info)
      case('NONZERO','nonzero')
        call del_zero_terms(flist,mode,op_info,1d-12)
      case('RANK','rank')
        call select_rank(flist,labels,nlabels,mode,op_info)
      case default
        call quit(1,'form_select_special','unknown type: "'
     &       //trim(type)//'"')
      end select

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = 'XXX'
      call write_form_list(f_output%fhand,flist,'XXX')

      if (ntest.ge.100) then
        write(lulog,*) 'final formula:'
        call print_form_list(lulog,flist,op_info)
      end if

      call dealloc_formula_list(flist)

      return
      end

