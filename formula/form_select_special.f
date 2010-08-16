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
        call write_title(luout,wst_dbg_subr,
     &       'here speaks form_select_special')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist)

      if (ntest.ge.150) then
        write(luout,*) 'initial formula:'
        call print_form_list(luout,flist,op_info)
      end if

      select case(trim(type))
      case('F12x','f12x')
        call select_f12x(flist,labels,nlabels,mode,op_info)
      case('MRCC')
        call select_mrcc_lag(flist,labels,nlabels,mode,op_info)
      case('OPT1','opt1')
        call select_xsp_opt1(flist,labels,nlabels,mode,op_info)
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
        write(luout,*) 'final formula:'
        call print_form_list(luout,flist,op_info)
      end if

      call dealloc_formula_list(flist)

      return
      end

