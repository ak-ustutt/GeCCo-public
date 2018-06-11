*----------------------------------------------------------------------*
      subroutine form_itf(f_input,name_out,form_out,op_info)
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
     &     name_out,
     &     form_out
      type(formula), intent(inout) ::
     &     f_input
      type(operator_info), intent(in) ::
     &     op_info

      type(filinf) ::
     &     fitf,
     &     fform
      type(formula_item) ::
     &     flist
      logical ::
     &     print_form
      integer ::
     &     e      ! Exit satatus from python

      print_form=form_out.ne.'##not_set##'

      call write_title(lulog,wst_section,'Translating formula to ITF')

      ! Open file to write to
      call file_init(fitf,name_out,ftyp_sq_frm,0)
      call file_open(fitf)

      if (print_form) then
        ! Open optional formulae file
        call file_init(fform,form_out,ftyp_sq_frm,0)
        call file_open(fform)
        write(fform%unit,*) "Generated formulae"
      end if

      ! Read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      if (iprlvl.gt.0) then
        write(lulog,*) 'Translating ',trim(f_input%label),
     &                 ' formulae into ITF algo code'
        write(lulog,*) 'Code written to file: ',trim(fitf%name)
      end if

      ! Translate formula list into ITF code
      call print_itf(fitf%unit,flist,op_info,print_form,fform%unit)

      ! Close file
      call file_close_keep(fitf)
      call dealloc_formula_list(flist)

      if (print_form) then
        write(fform%unit,*) "End formulae"
        call file_close_keep(fform)
      end if

      ! Process ITF lines with external python script
!      call execute_command_line("python3 ../test_python/process.py",
!     &                          exitstat=e)
      
      return
      end
