*----------------------------------------------------------------------*
      subroutine evaluate(
     &     label_form,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     just evaluate given formula (must be factorized)
*
*     op_info:  operator definitions and files
*     str_info: string information (to be passed to subroutines)
*     orb_info: orbital space information (to be passed)
*
*     andreas, Aug. 2007
*
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
     &     ifree, nout, iout, idx
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

      ifree = mem_setmark('evaluate')

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'entered evaluate')
      end if

      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,'evaluate',
     &     'did not find formula '//trim(label_form))
      f_eval => form_info%form_arr(idx)%form

      ! read formula
      call read_form_list(f_eval%fhand,fl_eval)

      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_eval,op_info)

      nout = depend%ntargets
      allocate(xret(nout))

      ! call the scheduler
      call frm_sched(xret,fl_eval,depend,0,0,
     &               op_info,str_info,strmap_info,orb_info)

c dbg
      call get_argument_value('method.R12','pz_eval',lval=pz_eval)
      if(pz_eval)then
        do iout = 1, nout
          idx = depend%idxlist(iout)
          if(trim(op_info%mel_arr(idx)%mel%op%name).eq.op_z_inter.or.
c          if(trim(op_info%mel_arr(idx)%mel%op%name).eq.op_z_test.or.
     &     trim(op_info%mel_arr(idx)%mel%op%name).eq.op_p_inter)then
            print *,'EVALUATION OF P/Z'
            call wrt_mel_seq(luout,
     &           op_info%mel_arr(idx)%mel,
     &           1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
     &           str_info,orb_info)
          endif
        enddo
      endif
c dbg

      if (iprlvl.ge.5) then
        call write_title(luout,wst_title,
     &       'norms/values of output operators')
        do iout = 1, nout
          write(luout,'(">>>",1x,i4," --> ",g16.10)') iout, xret(iout)
        end do
      end if

      if (ntest.ge.1000) then
        do iout = 1, nout
          idx = depend%idxlist(iout)
          write(luout,*) 'dump of result for ',
     &         trim(op_info%mel_arr(idx)%mel%label)
          call wrt_mel_file(luout,5,
     &       op_info%mel_arr(idx)%mel,
     &       1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
     &       str_info,orb_info)
        end do
      end if

      call clean_formula_dependencies(depend)

      ifree = mem_flushmark()

      return
      end

