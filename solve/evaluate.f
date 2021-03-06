
*----------------------------------------------------*
!>     just evaluate given formula (must be factorized)
!!
!!     this subroutine is to be called from userspace
!!     @param label_form label of the formula
!!     @param init  if the result operator should be set to zero before evaluation
!!     @param op_info  operator definitions and files
!!     @param str_info string information (to be passed to subroutines)
!!     @param orb_info orbital space information (to be passed)
*
*     andreas, Aug. 2007
*
*----------------------------------------------------------------------*
      subroutine evaluate(
     &     label_form,init,
     &     op_info,form_info,str_info,strmap_info,orb_info)

*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'
      include 'ifc_input.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     ntest = 000

      character(*) ::
     &     label_form
      logical ::
     &     init
      type(formula_info) ::
     &     form_info
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      integer ::
     &     ifree, nout, iout, idx, iprint
      logical ::
     &     pz_eval

      type(formula), pointer ::
     &     f_eval
      type(formula_item) ::
     &     fl_eval
      type(dependency_info) ::
     &     depend

      real(8), pointer ::
     &     xret(:)

      integer, external ::
     &     idx_formlist

      iprint=max(iprlvl,ntest)

      ifree = mem_setmark('evaluate')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'entered evaluate')
      end if

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'evaluate',
     &     'did not find formula '//trim(label_form))
      f_eval => form_info%form_arr(idx)%form

      ! read formula
      call read_form_list(f_eval%fhand,fl_eval,.true.)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_eval,op_info)

      nout = depend%ntargets
      allocate(xret(nout))

      ! call the scheduler - I have now set NO_SKIP to true, as skipping of targets is
      ! nearly always an unwanted effect for the EVALUATE statement
      call frm_sched(xret,fl_eval,depend,0,0,
     &               init,.true.,op_info,str_info,strmap_info,orb_info)

      if (iprint.ge.10) then
        call write_title(lulog,wst_title,
     &       'norms/values of output operators')
        do iout = 1, nout
          write(lulog,'(">>>",1x,i4," --> ",g20.14,x,"(",a,")")')
     &        iout, xret(iout), trim(label_form)

        end do
      end if

      if (iprint.ge.30) then
        do iout = 1, nout
          idx = depend%idxlist(iout)
          write(lulog,*) 'dump of result for ',
     &         trim(op_info%mel_arr(idx)%mel%label)
          call wrt_mel_file(lulog,5,
     &       op_info%mel_arr(idx)%mel,
     &       1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        end do
      end if

      deallocate(xret)
      call dealloc_formula_list(fl_eval)
      call clean_formula_dependencies(depend)

      ifree = mem_flushmark()

      return
      end

