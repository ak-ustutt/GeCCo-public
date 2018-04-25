*----------------------------------------------------------------------*
      subroutine form_itf(f_input,name_output,op_info)
*----------------------------------------------------------------------*
*     Driver for outputing ITF algo code
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

      character(len=*), intent(in) ::
     &     name_output
      type(formula), intent(inout) ::
     &     f_input
      type(operator_info), intent(in) ::
     &     op_info

      type(filinf) ::
     &     ffitf
      type(formula_item) ::
     &     flist

      call write_title(lulog,wst_section,'Translating formula to ITF')

      call file_init(ffitf,name_output,ftyp_sq_frm,0)
      call file_open(ffitf)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      if (iprlvl.gt.0) then
        write(lulog,*) 'Translating ',trim(f_input%label),
     &                 ' formulae into ITF algo code'
        write(lulog,*) 'Code written to file: ',trim(ffitf%name)
      end if

      call print_itf(ffitf%unit,flist,op_info)

      write(ffitf%unit,*) "Hello world"

      call file_close_keep(ffitf)

      call dealloc_formula_list(flist)
      
      return
      end
