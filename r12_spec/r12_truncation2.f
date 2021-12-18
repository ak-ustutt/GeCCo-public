      subroutine r12_truncation2(flist,trunc_type,trunc_t1x,RGRc,
     &     idxham,idxsbar,idxsop,op_info)
*----------------------------------------------------------------------*
*     new routine for truncated R12 expansions
*     we work with the raw formula that still contains "S"
*     we can much more easily spot T2', T3' etc. pp.
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     trunc_type, idxsbar, idxsop, idxham, 
     &     trunc_t1x, RGRc

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op, rank,
     &     ord_t, ord_ham, ord_tx, ord_tbx,
     &     ntop, ntx, nt1x, nham, ntbar, ntbx, nt1bx
      character*64 ::
     &     op_name
      logical ::
     &     xtern, diag_xx, offd_xx, t1x_involved

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

c      integer, external ::
c     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'pert_trunction')
      endif

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

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! find out:
          ! - number and X-indices of T operators
          ! - perturbation order of H
          ntop    = 0
          ntx     = 0
          nt1x    = 0
          ntbar   = 0
          ntbx    = 0
          nt1bx   = 0
          nham    = 0
          ord_ham = 0
          ord_tx  = 0
          ord_tbx = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxsop) then
              rank = op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)
              xtern=op_info%op_arr(idx_op)%
     &              op%ihpvca_occ(IEXTR,1,iblk_op).gt.0 
              if (.not.xtern) ntop = ntop+1
              if (xtern)      ntx  = ntx +1
              if (xtern.and.rank.eq.1) nt1x  = nt1x +1
              if (xtern)      ord_tx = ord_tx + rank-1
            end if
            if (idx_op.eq.idxsbar) then
              rank = op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)
              xtern=op_info%op_arr(idx_op)%
     &              op%ihpvca_occ(IEXTR,2,iblk_op).gt.0 
              if (.not.xtern) ntbar = ntbar+1
              if (xtern)      ntbx  = ntbx +1
              if (xtern.and.rank.eq.1) nt1bx  = nt1bx +1
              if (xtern)      ord_tbx = ord_tbx + rank-1
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              ord_ham = ord_ham
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'r12_truncation2','strange: nham.ne.1')

          ! well, old behaviour only for trunc_t1x.eq.-1
          ! else, we forget about these operators and treat them later:
          if (trunc_t1x.ne.-1) then
            ntx = ntx-nt1x
            ntbx = ntbx-nt1bx
            nt1x = 0
            nt1bx = 0
          end if

          ! always linear in T1X
          delete = nt1x.gt.1

          if (trunc_type.ne.1) then
            ! (F12):
            ! always linear in external terms
            delete = delete.or.ntx.gt.1
            ! external projection:
            ! diagonal coupling case (e.g. TBAR2' with T2')
            diag_xx = ntbx.gt.0.and.ntx.gt.0.and.ord_tbx.eq.ord_tx
            ! off-diagonal coupling case (e.g. TBAR2' with T3')
            offd_xx = ntbx.gt.0.and.ntx.gt.0.and.ord_tbx.ne.ord_tx
            t1x_involved = nt1bx.gt.0.or.nt1x.gt.0
            ! on the diagonal, we only consider the F matrix
            if (diag_xx) then
              delete = delete.or.ord_ham+ntop.gt.0
            end if
            ! for the off-diagonal, it looks that we better keep
            ! all the terms; however, for T1' we previously defined
            ! to omit these terms, so for compatibility:
            ! (well, for RGRc==0 we drop the offd_xx terms)
            if (offd_xx.and.(t1x_involved.or.RGRc.eq.0)) then
              delete = delete.or.ord_ham+ntop.gt.0
            end if
            
          else if (trunc_type.eq.1) then
            ! linearized R12:
            delete = delete.or.ntx.gt.1
          end if
          if (trunc_type.eq.2) then
            ! no R12 at all:
            delete = delete.or.ntx.gt.0
            delete = delete.or.ntbx.gt.0
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'r12_truncation2','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
