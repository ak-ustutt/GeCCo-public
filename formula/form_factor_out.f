*----------------------------------------------------------------------*
      subroutine form_factor_out(f_output,f_input,
     &                      title,
     &                      nintm,label_f_intm,split,
     &                      op_info,form_info)
*----------------------------------------------------------------------*
*
*     driver for factorizing intermediates from input formula
*
*     split: allow for factoring out by splitting certain terms
*
*     f_input and f_output may be identical
*
*     andreas, jan 2008 (update: 2021)
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
      include 'ifc_formula.h' 

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nintm
      logical, intent(in) ::
     &     split
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
     &     same, transpose, debug
      integer ::
     &     iintm, idx, len, nrpl, nspl
      character ::
     &     name*(form_maxlen_label*2)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      type(filinf), pointer ::
     &     ffintm
      type(formula_item) ::
     &     flist, fl_intm

      integer, external ::
     &     idx_formlist, form_count

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &     'form_factor_out reports')
        write(lulog,*) ' f_input = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
        do iintm = 1, nintm
          write(lulog,*) iintm,label_f_intm(iintm)
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      ! loop over intermediates
      do iintm = 1, nintm

        ! look for transposition label
        len = len_trim(label_f_intm(iintm))

        ! trigger debug - uncomment (and set "ntest" in subroutines)
        debug = .false. !label_f_intm(iintm)(1:6) == 'F_T1T2'

        transpose = (label_f_intm(iintm)(len-1:len).eq.'^+') 
        if (transpose) len = len-2
         
        ! get index
        idx = idx_formlist(label_f_intm(iintm)(1:len),form_info)

        if (idx.le.0)
     &       call quit(1,'form_factor_out',
     &       'formula label not on list: '//trim(label_f_intm(iintm)))
        
        if (ntest.ge.100)
     &       write(lulog,*)
     &       'now factorizing: ',trim(label_f_intm(iintm)),
     &       ' transpose: ',transpose

        ffintm => form_info%form_arr(idx)%form%fhand
        if (.not.associated(ffintm))
     &       call quit(1,'form_factor_out',
     &       'formula file does not exist for '//
     &       trim(form_info%form_arr(idx)%form%label))

        call init_formula(fl_intm)
        call read_form_list(ffintm,fl_intm,.true.)
c dbg
c        ! Can be used to bypass an error in find_possible_subexpr
c        call reorder_formula(fl_intm,op_info)
c        print *,'reordered intermediate'
c        call print_form_list(lulog,fl_intm,op_info)
c dbgend

        if (transpose)
     &       call transpose_formula(fl_intm,op_info)

        call atim_csw(cpu0,sys0,wall0)
        call factor_out_subexpr2(flist,fl_intm,split,nrpl,nspl,
     &                           op_info,debug)
        call atim_csw(cpu,sys,wall)

        if (iprlvl.ge.2) write(lulog,
     &       '(x,a40,": ",i6," replacements with ",i6," splits")')
     &                  trim(label_f_intm(iintm)),nrpl,nspl
        if (iprlvl.ge.2) then
          len = form_count(flist)
          write(lulog,'(x,"formula reduced to ",i10," items")') len
        end if

        ! consider sum_terms here, if there were splits
        if (nspl.gt.0) call del_zero_terms(flist,'sum',op_info,1d-12)
        
        if (iprlvl.ge.2)
     &     call prtim(lulog,'factoring out',
     &     cpu-cpu0,sys-sys0,wall-wall0)

        call dealloc_formula_list(fl_intm)

      end do

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

      if (ntest.ge.1000) then
        call write_title(lulog,wst_around_double,'Factored formula:')
        call print_form_list(lulog,flist,op_info)
      end if

      call dealloc_formula_list(flist)
      
      return
      end
