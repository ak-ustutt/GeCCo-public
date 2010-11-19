*----------------------------------------------------------------------*
      subroutine tex_formula_drv(f_input,name_output,op_info)
*----------------------------------------------------------------------*
*
*     driver for printing formula
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
c      include 'def_contraction_list.h'
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
     &     fftex

      type(formula_item) ::
     &     flist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'tex_formula_drv')
        write(luout,*) ' f_input  = ',trim(f_input%label)
      end if

      call file_init(fftex,name_output,ftyp_sq_frm,0)
      call file_open(fftex)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! TeX it ...
      if (iprlvl.gt.0) then
        write(luout,*) 'A TeX formatted version of formula ',
     &       trim(f_input%label)
        write(luout,*) 'is written to file: ',trim(fftex%name)
      end if

      call tex_form_list(fftex%unit,flist,op_info)

      call file_close_keep(fftex)

      call dealloc_formula_list(flist)
      
      return
      end
