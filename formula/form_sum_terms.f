*----------------------------------------------------------------------*
      subroutine form_sum_terms(f_output,f_input,
     &                      title,thresh,
     &                      op_info)
*----------------------------------------------------------------------*
*     sum up terms
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest =  00

      character(len=*), intent(in) ::
     &     title
      real(8), intent(in) ::
     &     thresh
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     same
      integer ::
     &     idxres
      character(len=form_maxlen_label*2) ::
     &     name

      type(formula_item), target ::
     &     flist

      integer, external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'here speaks form_sum_terms')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
        write(lulog,*) ' thresh   = ',thresh
      end if

      call atim_csw(cpu0,sys0,wall0)

      same = trim(f_input%label).eq.trim(f_output%label)


      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      !call sum_terms(flist,op_info) <- is called below by option "sum"
      call del_zero_terms(flist,'sum',op_info,thresh)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
        f_output%comment = trim(title)
      end if
      call write_form_list(f_output%fhand,flist,title)

      call dealloc_formula_list(flist)

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'Sum terms',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      
      return
      end
