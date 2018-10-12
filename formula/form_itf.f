*----------------------------------------------------------------------*
      subroutine form_itf(f_input,name_out,form_out,multi,op_info)
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

      character(len=*), intent(in) ::
     &     name_out,    ! Output file name for ITF algo code
     &     form_out     ! Output file name for GeCco formulae code
      type(formula), intent(inout) ::
     &     f_input
      type(operator_info), intent(in) ::
     &     op_info

      type(filinf) ::
     &     fline,       ! Temporary file which contrains ITF binary contractions
     &     fitf,        ! File for ITF algo code
     &     fform        ! File for GeCco formulae
      type(formula_item) ::
     &     flist        ! Linked list of binary contractions
      logical ::
     &     print_form,  ! If true, outputs GeCco formulae to file
     &     multi        ! Flag which is passed to python processer, false if a single-ref calculation
      integer ::
     &     e            ! Exit satatus from python
      character(len=100) ::
     &     exe_line     ! Line for shell to execute
      character(len=1) ::
     &     singlr       ! 0: single-reference; >0: multireference

      print_form=form_out.ne.'##not_set##'

      call write_title(lulog,wst_section,'Translating formula to ITF')

      ! Open tmp binary contraction file
      call file_init(fline,'bcontr.tmp',ftyp_sq_frm,0)
      call file_open(fline)

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
        write(lulog,*) ' Binary contractions written to: ',
     &                 trim(fline%name)
      end if

      ! Translate formula list into ITF binary contractions
      call print_itf(fline%unit,flist,op_info,print_form,fform%unit)

      ! Is this a single-refenece or a multireferecne calculation?
      if (.not. multi) then
         singlr = '0'
      else
         singlr = '1'
      endif

      ! Process ITF lines with python script
      ! Don't need to inialise the file, the python script takes care of
      ! that
      write(lulog,*) ' ITF algo file written to:       ',
     &               trim(name_out)

      exe_line='python3 $GECCO_DIR/itf_python/process.py -i '
     &          //trim(fline%name)//' -o '//trim(name_out)
     &          //' -s '//singlr
      write(lulog,*) "Executing: ", exe_line
      call execute_command_line(trim(exe_line),exitstat=e)

      ! Close binary contraction file
      call file_close_keep(fline)
      ! Deallocat link list
      call dealloc_formula_list(flist)

      if (print_form) then
        write(fform%unit,*) "End formulae"
        call file_close_keep(fform)
      end if

      return
      end
