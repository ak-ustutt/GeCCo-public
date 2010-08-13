*----------------------------------------------------------------------*
      subroutine set_formula(form,
     &     form_str,title,
     &     op_info)
*----------------------------------------------------------------------*
*
*     set up formula as described by form_str
*
*     very initial version for setting up sums of operators
*
*     basic syntax
*
*     <result_op> = <operator> [+ <operator> [...]]
* 
*     e.g.:
*
*     OP0 = OP1 + OP2 + OP3
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula), intent(inout), target ::
     &     form

      character*(*), intent(in) ::
     &     title, form_str

      type(operator_info), intent(in) ::
     &     op_info

      character(len=len_opname) ::
     &     label

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag
      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     ipos, jpos, idx_res, idx_op, len_form_str
      real(8) ::
     &     fac

      integer, external ::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'Setting formula')
        write(luout,*) ' form_str  = "',trim(form_str),'"'
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! advanced: call here a lexer to analyze the string

      ! for the moment very hand-made:
      ipos = index(trim(form_str),'=')
      if (ipos.eq.0) call quit(0,'set_formula','no "=" in str: '//
     &     trim(form_str))
      label(1:len_opname) = ' '
      label = trim(form_str(1:ipos-1))
      ipos = ipos+1

      if (ntest.ge.100)
     &     write(luout,*) 'result label = "',trim(label),'"'

      idx_res = idx_oplist2(label,op_info)
      if (idx_res.le.0)
     &     call quit(1,'set_formula',
     &                 'unknown label: "'//trim(label)//'"')

      ! initialize formula
      call init_formula(flist_lag)
      flist_pnt => flist_lag
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idx_res)
      flist_pnt => flist_pnt%next

      ! loop over contributions
      len_form_str = len_trim(form_str)
      do while(ipos.le.len_form_str)

        jpos = index(trim(form_str(ipos:)),'+')+ipos-1
        if (jpos.eq.ipos-1)
     &       jpos=len_form_str+1

        label(1:len_opname) = ' '
        label = trim(form_str(ipos:jpos-1))
        ipos = jpos+1

        if (ntest.ge.100)
     &       write(luout,*) 'next label = "',trim(label),'"'

        idx_op = idx_oplist2(label,op_info)
        if (idx_op.le.0)
     &       call quit(1,'set_formula',
     &                   'unknown label: "'//trim(label)//'"')

        fac = 1d0
        do while(associated(flist_pnt%next))
          flist_pnt => flist_pnt%next
        end do
        call set_primitive_formula(flist_pnt,idx_op,
     &       fac,idx_res,.false.,op_info) 

      end do
      ! reorder according to result blocks
      call reorder_formula(flist_lag,op_info)

      ! assign comment
      form%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form%label)
      call file_init(form%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form%fhand,flist_lag,
     &     form%comment)

      if (max(ntest,iprlvl).ge.50) then
        call write_title(luout,wst_around_double,'Generated formula:')
        call print_form_list(luout,flist_lag,op_info)
      end if

      call dealloc_formula_list(flist_lag)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'set formula',cpu-cpu0,sys-sys0,wall-wall0)

      end
