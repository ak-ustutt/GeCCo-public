*----------------------------------------------------------------------*
      subroutine transform(
     &     label_me_in,label_me_out,label_form,init,
     &     op_info,form_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*
*     transform each record on list with label "label_me_in" 
*     to a corresp. record on list with label "label_me_out"
*     using given formula (must be factorized)
*
*     op_info:  operator definitions and files
*     str_info: string information (to be passed to subroutines)
*     orb_info: orbital space information (to be passed)
*
*     andreas, Oct. 2014
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
      character(len=9), parameter ::
     &     i_am = 'transform'

      character(*) ::
     &     label_form, label_me_in, label_me_out
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
     &     ifree, nout, iout, idx,
     &     idxmel, recst_in, recnd_in, recst_out, recnd_out,
     &     idx_in, idx_out

      type(formula), pointer ::
     &     f_eval
      type(formula_item) ::
     &     fl_eval
      type(dependency_info) ::
     &     depend

      type(me_list), pointer ::
     &     me_in, me_out

      real(8), pointer ::
     &     xret(:)

      integer, external ::
     &     idx_formlist, idx_mel_list

      ifree = mem_setmark('transform')

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'entered transform')
      end if

      ! get lists
      idxmel = idx_mel_list(label_me_in,op_info)
      if (idxmel.le.0)
     &     call quit(1,i_am,'list not found: '//trim(label_me_in))
      me_in => op_info%mel_arr(idxmel)%mel
      idxmel = idx_mel_list(label_me_out,op_info)
      if (idxmel.le.0)
     &     call quit(1,i_am,'list not found: '//trim(label_me_out))
      me_out => op_info%mel_arr(idxmel)%mel

      ! get formula
      idx = idx_formlist(label_form,form_info)
      if (idx.le.0)
     &     call quit(1,i_am,
     &     'did not find formula '//trim(label_form))
      f_eval => form_info%form_arr(idx)%form

      ! read formula
      call read_form_list(f_eval%fhand,fl_eval,.true.)

c dbg
c     call print_form_list_short(lulog,fl_eval,op_info)
c dbg
      ! set dependency info for submitted formula list
      call set_formula_dependencies(depend,fl_eval,op_info)

      nout = depend%ntargets
      allocate(xret(nout))

      recst_in = me_in%fhand%active_records(1)
      recnd_in = me_in%fhand%active_records(2)
      recst_out = me_out%fhand%active_records(1)
      recnd_out = me_out%fhand%active_records(2)

      if (recnd_in-recst_in.ne.recnd_out-recst_out)
     &     call quit(1,i_am,'incompatible number of records')

      idx_out = recst_out-1
      do idx_in = recst_in, recnd_in
        idx_out = idx_out+1
        call switch_mel_record(me_in,idx_in)
        call switch_mel_record(me_out,idx_out)

        ! call the scheduler
        call frm_sched(xret,fl_eval,depend,0,0,
     &               init,.false.,op_info,str_info,strmap_info,orb_info)

        if (iprlvl.ge.10) then
          call write_title(lulog,wst_title,
     &         'norms/values of output operators')
          do iout = 1, nout
            write(lulog,'(">>>",1x,i4," --> ",g20.14)') iout, xret(iout)
          end do
        end if

        ! this was take from "evaluate" without adaption ...
        if (ntest.ge.1000) then
          do iout = 1, nout
            idx = depend%idxlist(iout)
            write(lulog,*) 'dump of result for ',
     &           trim(op_info%mel_arr(idx)%mel%label)
            call wrt_mel_file(lulog,5,
     &           op_info%mel_arr(idx)%mel,
     &           1,op_info%mel_arr(idx)%mel%op%n_occ_cls,
     &           str_info,orb_info)
          end do
        end if

      end do

      deallocate(xret)
      call dealloc_formula_list(fl_eval)
      call clean_formula_dependencies(depend)

      ifree = mem_flushmark()

      return
      end

