*----------------------------------------------------------------------*
      subroutine print_formula_drv(f_input,name_output,mode,op_info)
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
     &     name_output, mode
      type(formula), intent(inout) ::
     &     f_input
      type(operator_info), intent(in) ::
     &     op_info

      type(filinf) ::
     &     ffprint
      integer ::
     &     luprint

      type(formula_item) ::
     &     flist

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &     'tex_formula_drv')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
      end if

      if (name_output.eq.'stdout') then
        luprint = luout
      else if (name_output.eq.'stderr') then
        luprint = lulog
      else
        call file_init(ffprint,name_output,ftyp_sq_frm,0)
        call file_open(ffprint)
        luprint = ffprint%unit
      end if

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! print it ...
      if (mode.eq.'long'.or.mode.eq.'LONG') then
        call print_form_list(luprint,flist,op_info)
      else if  (mode.eq.'short'.or.mode.eq.'SHORT') then
        call print_form_list_short(luprint,flist,op_info)
      else
        call warn('PRINT_FORMULA','unknown mode: '//trim(mode))
      end if

      if (name_output.ne.'stdout')
     &     call file_close_keep(ffprint)

      call dealloc_formula_list(flist)
      
      return
      end
