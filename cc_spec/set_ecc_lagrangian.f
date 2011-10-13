*----------------------------------------------------------------------*
      subroutine set_ecc_lagrangian(form_lag,
     &     title,label_op,nlabels,ansatz,mode,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines a ECC-Lagrangian within the chosen operator space 
*
*     written by andreas, june 2008
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_lag

      character*(*), intent(in) ::
     &     title,
     &     label_op(*)
      integer, intent(in) ::
     &     nlabels

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local constants
      character(3), parameter ::
     &     op_hb_temp = '_HB'

      ! local variables

      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_lag, flist_hbar
      type(formula_item), pointer ::
     &     fl_pnt, fl_hbar_pnt
      integer, intent(in) ::
     &     ansatz
      character*(*) ::
     &     mode

      logical ::
     &     l_h0d
      integer ::
     &     nterms, ilabel, idx, 
     &     idxham,idxtbar,idxtop,idxtpt,idxtptbar,idxlag,
     &     idx_hb_temp, max_rank_t, max_ext_t, t1xmode

      type(operator), pointer::
     &     hb_temp_pnt

      integer, external ::
     &     idx_oplist2, maxxlvl_op

c prelim
      character(len=8) ::
     &     trmode
c prelim

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'Setting up ECC-Lagrangian')
      end if

      call atim_csw(cpu0,sys0,wall0)

      do ilabel = 1, nlabels
        idx = idx_oplist2(label_op(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_ecc_lagrangian',
     &       'label not on list: '//trim(label_op(ilabel)))
        if (ilabel.eq.1) idxlag = idx
        if (ilabel.eq.2) idxham = idx
        if (ilabel.eq.3) idxtbar = idx
        if (ilabel.eq.4) idxtop = idx
      end do

      ! initialize formula
      call init_formula(flist_lag)
      fl_pnt => flist_lag
      ! put [INIT] at the beginning
      call new_formula_item(fl_pnt,command_set_target_init,idxlag)
      fl_pnt => fl_pnt%next

c prelim
      trmode = '        '
      call get_argument_value('method.ECC','truncate',str=trmode)
      max_ext_t = 1
c prelim
      ! -------------------
      ! define and set Hbar
      ! -------------------
      call add_operator(op_hb_temp,op_info)
      idx_hb_temp = idx_oplist2(op_hb_temp,op_info)
      hb_temp_pnt => op_info%op_arr(idx_hb_temp)%op
      max_rank_t = maxxlvl_op(op_info%op_arr(idxtop)%op)
      if (trim(trmode).eq.'no') then
        call set_xop(hb_temp_pnt,op_hb_temp,.false.,
     &     0,4*max_rank_t-2,2*max_ext_t,0,0,orb_info)
      else
        call set_xop(hb_temp_pnt,op_hb_temp,.false.,
     &     0,max_rank_t+1,0,0,0,orb_info)
      end if

      ! expand e^{-T} H e^T 
      call init_formula(flist_hbar)
      fl_hbar_pnt => flist_hbar
      call new_formula_item(fl_hbar_pnt,command_set_target_init,
     &                      idx_hb_temp)
      fl_hbar_pnt => fl_hbar_pnt%next
      call expand_op_bch(fl_hbar_pnt,4,idx_hb_temp,
     &     1d0,-1,idxham,1d0,idxtop,-1,-1,op_info)
      call reorder_formula(flist_hbar,op_info)

      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'Hbar:')
        call print_form_list(luout,flist_hbar,op_info)
      end if

      ! -----------------------------
      ! expand <0| e^{Tbar} Hbar e^{-Tbar} |0>
      ! -----------------------------
      call expand_op_bch(fl_pnt,4,idxlag,
     &     1d0,-1,idx_hb_temp,-1d0,idxtbar,-1,-1,op_info)
      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'Initial Lagrangian:')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! replace Hbar by e^-T H e^T
      call expand_subexpr(flist_lag,flist_hbar,0,op_info)
      call sum_terms(flist_lag,op_info)
      if (ntest.ge.1000) then
        call write_title(luout,wst_title,'Expanded Lagrangian:')
        call print_form_list(luout,flist_lag,op_info)
        call print_form_list_p(luout,flist_lag,op_info)
      end if


      ! reform to doubly linked expression
      call make_doubly_connected(flist_lag,
     &     idxtbar,idxham,idxtop,op_info)
      

      if (trim(trmode).ne.'no') then
        call pert_truncation(flist_lag,trmode,
     &     idxtbar,idxham,idxtop,op_info)
      else
c        call quit(1,'set_ecc_lagrangian',
c     &       'only valid in truncation mode')
      end if
c quick'n'dirty:
      call get_argument_value('method.ECC','T1ext',ival=t1xmode)
      if (t1xmode.gt.0) then
        write(trmode,'("ord",i1)') t1xmode
        call get_argument_value('method.ECC','H0_T1ext',ival=t1xmode)
        write(trmode(6:),'(i1)') t1xmode
        call get_argument_value('method.ECC','H0d',lval=l_h0d)
        if (l_h0d) trmode(7:7) = 'd'
        call t1x_truncation(flist_lag,trmode,
     &       idxtbar,idxham,idxtop,op_info)
      end if

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_lag,op_info)
      end if

      ! assign comment
      form_lag%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_lag%label)
      call file_init(form_lag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_lag%fhand,flist_lag,
     &     form_lag%comment)

      call dealloc_formula_list(flist_lag)
      call dealloc_formula_list(flist_hbar)

      call del_operator(op_hb_temp,op_info)

      call atim_csw(cpu,sys,wall)
c      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'ECC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

c dbg
      if (ntest.ge.100) stop
c dbg

      end
