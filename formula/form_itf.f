*----------------------------------------------------------------------*
      subroutine form_itf(f_input,name_out,form_out,multi,process,kext,
     &                    tasks,init_res,itin,op_info)
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
      logical, intent(in) ::
     &     multi,       ! Flag which is passed to python processer, false if a single-ref calculation
     &     process,     ! Process bcontr.tmp file to create .itfaa file
     &     kext,        ! True if constructing INTpp tensor to contract in Kext
     &     tasks,       ! True if using ITF tasks to generate algos
     &     init_res,    ! Produce Init_residual algo code
     &     itin         ! Produce ITIN lines, or symmetrise residual at end

      type(filinf) ::
     &     fline,       ! Temporary file which contrains ITF binary contractions
     &     fitf,        ! File for ITF algo code
     &     fform        ! File for GeCco formulae
      type(formula_item) ::
     &     flist        ! Linked list of binary contractions
      logical ::
     &     print_form   ! If true, outputs GeCco formulae to file
      integer ::
     &     e            ! Exit satatus from python
      character(len=200) ::
     &     exe_line     ! Line for shell to execute
      character(len=100) ::
     &     flags        ! set the flag options for the python script

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0   ! Timing variables

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

      ! Set time variables
      call atim_csw(cpu0,sys0,wall0)

      ! Translate formula list into ITF binary contractions
      call print_itf(fline%unit,flist,itin,op_info,print_form,
     &               fform%unit)

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'Time to process formulae',
     &           cpu-cpu0,sys-sys0,wall-wall0)

      ! Is this a single-refenece or a multireferecne calculation?
      if (.not. multi) then
         flags = '--no-multi'
      else
         flags = '--multi'
      endif

      ! Is this a single-refenece or a multireferecne calculation?
      if (kext) then
         flags = trim(flags)//' --kext'
      else
         flags = trim(flags)//' --no-kext'
      endif

      ! Initalise the residual with only integral terms (T=0)
      if (init_res) then
         flags = trim(flags)//' --init-res'
      else
         flags = trim(flags)//' --no-init-res'
      endif

      ! Process ITF lines with python script
      ! Don't need to inialise the file, the python script takes care of
      ! that
      write(lulog,*) ' ITF algo file written to:       ',
     &               trim(name_out)

      if (process) then
         ! Quick and dirty scan through tmp binary contraction file to
         ! find if any repeated lines can be replaced by a factor
         call atim_csw(cpu0,sys0,wall0)

         exe_line='python3 $GECCO_DIR/itf_python/simplify.py -i '
     &             //trim(fline%name)//' -o bcontr2.tmp'
         write(lulog,*) "Executing: ", exe_line
         call execute_command_line(trim(exe_line),exitstat=e)

         call atim_csw(cpu,sys,wall)
         call prtim(lulog,'Time to simplify ITF code',
     &              cpu-cpu0,sys-sys0,wall-wall0)

         if (e > 0) then
            write(lulog,*) "Error in executing simplify.py"
            call quit(1,'Please check the bcontr.tmp file')
         end if


         exe_line='python3 $GECCO_DIR/itf_python/remove_extra_lines.py'
     &             //' bcontr2.tmp'
         write(lulog,*) "Executing: ", exe_line
         call execute_command_line(trim(exe_line),exitstat=e)

         call atim_csw(cpu,sys,wall)
         call prtim(lulog,'Time to remove exta lines in ITF code',
     &              cpu-cpu0,sys-sys0,wall-wall0)

         if (e > 0) then
            write(lulog,*) "Error in executing remove_exta_lines.py"
            call quit(1,'Please check the bcontr2.tmp file')
         end if


         ! Main ITF algo file processing
         call atim_csw(cpu0,sys0,wall0)

         exe_line='python3 $GECCO_DIR/itf_python/process.py -i '
     &             //'bcontr3.tmp -o '//trim(name_out)
     &             //' '//trim(flags)
         write(lulog,*) "Executing: ", exe_line
         call execute_command_line(trim(exe_line),exitstat=e)

         call atim_csw(cpu,sys,wall)
         call prtim(lulog,'Time to create ITF algo file',
     &              cpu-cpu0,sys-sys0,wall-wall0)

         if (e > 0) then
            write(lulog,*) "Error in executing process.py"
            call quit(1,'Please check the bcontr3.tmp file')
         end if

!         ! Quick and dirty scan to remove unnecessary load and drop lines
!         call atim_csw(cpu0,sys0,wall0)
!
!         exe_line='python3 $GECCO_DIR/itf_python/reduce.py -i '
!     &             //trim(name_out)
!         write(lulog,*) "Executing: ", exe_line
!         call execute_command_line(trim(exe_line),exitstat=e)
!
!         call atim_csw(cpu,sys,wall)
!         call prtim(lulog,'Time to reduce ITF algo file',
!     &              cpu-cpu0,sys-sys0,wall-wall0)
!
!         if (e > 0) then
!            write(lulog,*) "Error in executing reduce.py"
!            call quit(1,'Please check the ITF algo file')
!         end if
!
!         exe_line='mv tmp.itfaa '//trim(name_out)
!         write(lulog,*) "Executing: ", exe_line
!         call execute_command_line(trim(exe_line),exitstat=e)
!
!         if (e > 0) then
!            write(lulog,*) "Error in renaming files"
!            call quit(1,'Please check the ITF algo file')
!         end if
      end if

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
