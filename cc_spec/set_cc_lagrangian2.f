*----------------------------------------------------------------------*
      subroutine set_cc_lagrangian2(form_cclag,
     &     title,label_oplcc,label_opham,label_optbar,label_opt,
     &     op_info)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines a CC-Lagrangian within the chosen operator space 
*
*     new version using gen_contr()
*     written by andreas, june 2007
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
c      include 'par_formnames_gen.h'
c prelim
      include 'ifc_input.h'
c prelim

      type(formula), intent(inout), target ::
     &     form_cclag

      character*(*), intent(in) ::
     &     title,
     &     label_opham,label_optbar,label_opt,label_oplcc

      type(operator_info), intent(in) ::
     &     op_info

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag
      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     nterms,
     &     idxham,idxtbar,idxtop,idxlcc

      integer, external ::
     &     idx_oplist2

c prelim
      character(len=8) ::
     &     trmode
c prelim

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'Setting up CC-Lagrangian')
        write(luout,*) ' op_ham  = ',trim(label_opham)
        write(luout,*) ' op_tbar = ',trim(label_optbar)
        write(luout,*) ' op_t    = ',trim(label_opt)
        write(luout,*) ' op_lcc  = ',trim(label_oplcc)
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! get indices
      idxlcc  = idx_oplist2(label_oplcc,op_info)
      idxham  = idx_oplist2(label_opham,op_info)
      idxtbar = idx_oplist2(label_optbar,op_info)
      idxtop  = idx_oplist2(label_opt,op_info)
      if (idxlcc.lt.0.or.idxham.lt.0.or.idxtbar.lt.0.or.idxtop.lt.0)
     &     call quit(1,'set_cc_lagrangian2',
     &     'required operators are not yet defined')

      ! initialize formula
      call init_formula(flist_lag)
      flist_pnt => flist_lag
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxlcc)
      flist_pnt => flist_pnt%next

      ! expand <0|(1+Tbar) e^{-T} H e^T|0> =
      ! <0| e^{-T} H e^T|0> +
      call expand_op_bch(flist_pnt,2,idxlcc,
     &     1d0,-1,idxham,1d0,idxtop,1,2,op_info)

      ! advance pointer
      do while(associated(flist_pnt%next))
        flist_pnt => flist_pnt%next
      end do
      ! <0|Tbar e^{-T} H e^T|0>
      call expand_op_bch(flist_pnt,4,idxlcc,
     &     1d0,idxtbar,idxham,1d0,idxtop,1,-1,op_info)

      ! insert here procedure to produce approx. expansions      
c prelim
      trmode = '        '
      call get_argument_value('method.CC','truncate',str=trmode)
      if (trim(trmode).ne.'no')
     &     call pert_truncation(flist_lag,trmode,
     &     idxtbar,idxham,idxtop,op_info)
c prelim

      ! post_processing and term counting:
      call cc_form_post(flist_lag,nterms,idxtbar,idxham,idxtop,iprlvl,
     &     op_info)

      ! assign comment
      form_cclag%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_cclag%label)
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cclag%fhand,flist_lag,
     &     form_cclag%comment)

      if (iprlvl.ge.50
c prelim
     &     .or.trim(trmode).ne.'no'
c prelim
     &     ) then
        call write_title(luout,wst_around_double,'Generated formula:')
        call print_form_list(luout,flist_lag,op_info)
      end if

      call dealloc_formula_list(flist_lag)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      end
