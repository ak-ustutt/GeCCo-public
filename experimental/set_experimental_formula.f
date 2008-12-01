*----------------------------------------------------------------------*
      subroutine set_experimental_formula(form_exp,
     &     title,label,nlabels,ansatz,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     expands the following operator product given the operators L,T,Op:
*     Formula = <0| (1+L) exp(-T) Op exp(T) |0>
*     ansatz, approx can be freely used
*
*     matthias, 2008
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

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
     &     form_exp
      integer, intent(in) ::
     &     ansatz, nlabels
      character(*), intent(in) ::
     &     label(nlabels), title, approx

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local parameters
      integer, parameter ::
     &     maxterm = 4

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist_exp
      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     ilabel, idx, idx_l, idx_t, idx_resplag, idx_op,
     &     iprint, iterm, op_idx(6+maxterm), op_num(6+maxterm),
     &     connect(2*maxterm+2), switch, no_l
 
      integer, external::
     &     idx_oplist2

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'output from set_experimental')
        write(luout,*) ' ansatz = ',ansatz
        write(luout,*) ' approx = ',trim(approx)
      end if

      call atim_csw(cpu0,sys0,wall0)
      ! transform labels into indices in this way
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_exp',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idx_resplag = idx
        if (ilabel.eq.2) idx_l = idx
        if (ilabel.eq.3) idx_op = idx
        if (ilabel.eq.4) idx_t = idx
      end do

      ! initialize formula
      call init_formula(flist_exp)
      flist_pnt => flist_exp
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,
     &    command_set_target_init,idx_resplag)

      ! expand <0|exp(-T) Op exp(T) |0> (no_l=1)
      ! +  <0|L exp(-T) Op exp(T) |0> (no_l=0)
      do no_l = 1,0,-1
        op_idx(1) = idx_resplag
        op_idx(2) = idx_l
        op_idx(3-no_l) = idx_op
        op_idx(4-no_l) = idx_resplag
        op_num(1:4) = (/1,2,3,1/)
        if (no_l.eq.1) op_num(3) = 1
        connect(1) = 0
        switch = 1      ! needed for terms without T: connect(1:1)=0
        do iterm = 1,maxterm+1
          do while(associated(flist_pnt%next))
            flist_pnt => flist_pnt%next
          end do
          call expand_op_product2(flist_pnt,idx_resplag,
     &         1.d0,iterm+3-no_l,iterm+2-no_l,
     &         op_idx(1:iterm+3-no_l),
     &         op_num(1:iterm+3-no_l),
     &         -1,-1,
     &         connect(1:2*iterm-switch),iterm-1,
     &         0,0,
     &         0,0,
     &         op_info)
          op_idx(iterm+3-no_l) = idx_t
          op_idx(iterm+4-no_l) = idx_resplag
          op_num(iterm+3-no_l) = iterm+3-no_l
          op_num(iterm+4-no_l) = 1
          connect(2*iterm-1) = 3-no_l
          connect(2*iterm) = iterm+3-no_l
          switch = 2
        end do
      end do

      if (iprint.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist_exp,op_info)
      end if

      ! assign comment
      form_exp%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_exp%label)
      call file_init(form_exp%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_exp%fhand,flist_exp,
     &     form_exp%comment)

      call dealloc_formula_list(flist_exp)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'set_exp',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
