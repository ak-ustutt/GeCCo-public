*----------------------------------------------------------------------*
      subroutine form_expand_subexpr(f_output,f_input,
     &                      title,
     &                      nintm,label_f_intm,
     &                      op_info,form_info)
*----------------------------------------------------------------------*
*
*     driver for expanding intermediates in input formula
*
*     f_input and f_output may be identical
*
*     andreas, jan 2008
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

      integer, intent(in) ::
     &     nintm
      character(*), intent(in) ::
     &     title,
     &     label_f_intm(nintm)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info), intent(in) ::
     &     op_info
      type(formula_info), intent(in) ::
     &     form_info

      logical ::
     &     same, transpose
      integer ::
     &     iintm, idx, len
      character ::
     &     name*(form_maxlen_label*2)

      type(filinf), pointer ::
     &     ffintm
      type(formula_item) ::
     &     flist, fl_intm

      integer, external ::
     &     idx_formlist

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'form_expand_subexpr reports')
        write(luout,*) ' f_input = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        do iintm = 1, nintm
          write(luout,*) iintm,label_f_intm(iintm)
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! loop over intermediates
      do iintm = 1, nintm
        call atim_csw(cpu0,sys0,wall0)

        ! look for transposition label
        len = len_trim(label_f_intm(iintm))
        transpose = (label_f_intm(iintm)(len-1:len).eq.'^+') 
        if (transpose) len = len-2
         
        ! get index
        idx = idx_formlist(label_f_intm(iintm)(1:len),form_info)

        if (idx.le.0)
     &       call quit(1,'form_expand_subexpr',
     &       'formula label not on list: '//trim(label_f_intm(iintm)))
        
        if (ntest.ge.100)
     &       write(luout,*)
     &       'now expanding: ',trim(label_f_intm(iintm)),
     &       ' transpose: ',transpose

        ffintm => form_info%form_arr(idx)%form%fhand
        if (.not.associated(ffintm))
     &       call quit(1,'form_factor_out',
     &       'formula file does not exist for '//
     &       trim(form_info%form_arr(idx)%form%label))

        call init_formula(fl_intm)
        call read_form_list(ffintm,fl_intm,.true.)

        if (transpose)
     &       call transpose_formula(fl_intm,op_info)

        call expand_subexpr(flist,fl_intm,.false.,op_info)

        call sum_terms(flist,op_info)

        call atim_csw(cpu,sys,wall)
        call prtim(luout,'expand',cpu-cpu0,sys-sys0,wall-wall0)

        call dealloc_formula_list(fl_intm)

      end do

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (ntest.ge.100) then
        call write_title(luout,wst_around_double,'Expanded formula:')
        call print_form_list(luout,flist,op_info)
      end if

      call dealloc_formula_list(flist)

      return
      end
