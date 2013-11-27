      subroutine select_repl_htt(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     Replace H1bar by H in terms like [[H,T],T] (or more Ts)
*     Works only if H1bar and H have identical blocks
*
*     matthias, Jun 2011
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_input.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     nlabels
      character(len=*), intent(in) ::
     &     labels(nlabels), mode

      logical ::
     &     error
      integer ::
     &     idxtop, idxham, idxhbar, ntop, iterm, ii, ivtx, nvtx, idx_op,
     &     ivtxham
      integer ::
     &     idxop(nlabels)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next
      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_repl_htt')
c        write(lulog,*) 'mode = ',trim(mode)
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.3
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.3)
     &       call quit(1,'select_repl_htt','need 3 labels')
        call quit(1,'select_repl_htt','Labels not on list!')
      end if

      idxhbar  = idxop(1)
      idxham  = idxop(2)
      idxtop  = idxop(3)

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)

          iterm = iterm + 1
c dbg
c          print *,'current term: ',iterm
c dbgend
          contr => form_pnt%contr
          nvtx = contr%nvtx
          vertex => contr%vertex

          ! find out:
          ! - number of T operators
          ntop  = 0
          ivtxham = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) ntop = ntop+1
            if (idx_op.eq.idxhbar) ivtxham = ivtx
          end do

          ! if two or more T: replace H1bar by H
          if (ntop.ge.2.and.ivtxham.gt.0) then
            vertex(ivtxham)%idx_op = idxham
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_repl_htt','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
