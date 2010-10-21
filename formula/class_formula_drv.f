*----------------------------------------------------------------------*
      subroutine class_formula_drv(f_input,name_output,op_info,mode)
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
      integer, intent(in) ::
     &     mode

      type(filinf) ::
     &     ffprint
      integer ::
     &     luprint

      type(formula_item) ::
     &     flist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'class_formula_drv')
        write(luout,*) ' f_input  = ',trim(f_input%label)
      end if

      if (name_output.ne.'stdout') then
        call file_init(ffprint,name_output,ftyp_sq_frm,0)
        call file_open(ffprint)
        luprint = ffprint%unit
      else
        luprint = luout
      end if

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! print it ...
      call class_form_list(luprint,flist,op_info,mode)

      if (name_output.ne.'stdout')
     &     call file_close_keep(ffprint)

      call dealloc_formula_list(flist)
      
      return
      end
