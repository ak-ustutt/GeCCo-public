*----------------------------------------------------------------------*
*
!>     just evaluate given formula (must be factorized)
!!
!!     difference to evaluate(): formula head passed, not formula label
!!     @param fl_eval Head of the formula
!!     @param init  if the result operator should be set to zero before evaluation
!!     @param force if evaluation should occur even if the result op seems up to date.
!!     @param op_info  operator definitions and files
!!     @param str_info string information (to be passed to subroutines)
!!     @param orb_info orbital space information (to be passed)
*     andreas, Aug. 2007, matthias, Nov. 2010
**----------------------------------------------------------------------*
      subroutine evaluate2(fl_eval,init,force,
     &                     op_info,str_info,strmap_info,orb_info,
     &                     xret_out,get_xret)
*----------------------------------------------------------------------*      implicit none

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
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 000

      type(formula_item) ::
     &     fl_eval
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info
      logical ::
     &     get_xret, init, force
      real(8) ::
     &     xret_out(*)

      type(dependency_info) ::
     &     depend

      integer ::
     &     ifree, nout, iout, idx

      real(8), pointer ::
     &     xret(:)

      ifree = mem_setmark('evaluate2')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'entered evaluate2')
      end if

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_eval,op_info)

      nout = depend%ntargets
      allocate(xret(nout))

      ! call the scheduler
      call frm_sched(xret,fl_eval,depend,0,0,
     &            init,force,op_info,str_info,strmap_info,orb_info)

      if (iprlvl.ge.20) then
        call write_title(lulog,wst_title,
     &       'norms/values of output operators')
        do iout = 1, nout
          write(lulog,'(">>>",1x,i4," --> ",g20.14)') iout, xret(iout)
        end do
      end if

      if (ntest.ge.1000) then
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

      if (get_xret) then
        do iout = 1, nout
          xret_out(iout) = xret(iout)
        end do
      end if

      deallocate(xret)
      call clean_formula_dependencies(depend)

      ifree = mem_flushmark()

      return
      end

