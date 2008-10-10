      subroutine r12_truncation(flist,trunc_type,
     &     idxr12,idxham,idxtbar,idxtop,op_info)
*----------------------------------------------------------------------*
*     preliminary start-up for truncated R12 expansions
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 1000

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
     &     trunc_type, idxr12, idxtbar, idxtop, idxham

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op,
     &     ntop, ntx, nrdag, nr12, nham, ntbar, ntbx,
     &     ord_t, ord_ham,
     &     max_pert, max_comm, max_t1
      character*64 ::
     &     op_name
      logical ::
     &     dagger

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

c      integer, external ::
c     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'pert_trunction')
      endif

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(luout,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(luout,*) '[INIT_TARGET]'
        case(command_add_contribution)

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! find out:
          ! - number and X-indices of T operators
          ! - perturbation order of H
          ! - number of R+
          ! - number of R
          ntop  = 0
          ntx   = 0
          ntbar = 0
          ntbx  = 0
          nrdag = 0
          nr12  = 0
          nham    = 0
          ord_ham = 0
          ord_t   = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              ord_t = ord_t
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
              if (op_info%op_arr(idx_op)%
     &            op%ihpvca_occ(IEXTR,1,iblk_op).gt.0)
     &             ntx = ntx+1
            end if
            if (idx_op.eq.idxtbar) then
              ntbar = ntbar+1
              if (op_info%op_arr(idx_op)%
     &            op%ihpvca_occ(IEXTR,2,iblk_op).gt.0)
     &             ntbx = ntbx+1
            end if
            if (idx_op.eq.idxr12) then
              if (vertex(ivtx)%dagger) then
                nrdag = nrdag+1
              else
                nr12 = nr12+1
              end if
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              ord_ham = ord_ham
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'r12_truncation','strange: nham.ne.1')

          ! always linear in T1X
          delete = ntx.gt.1

          if (trunc_type.eq.0) then
            ! (R12):
            ! always linear in R12
            delete = delete.or.nr12+ntx.gt.1
            ! R12 projection:
            delete = delete.or.nrdag+ntbx.gt.0.and.
     &           ord_ham+ord_t.gt.0.and.nr12+ntx.gt.0
            
          else if (trunc_type.eq.1) then
            ! linearized R12:
            delete = delete.or.nr12.gt.1
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*) 'Deleted formula item:'
              call prt_contr2(luout,form_pnt%contr,op_info)
              write(luout,*) 'nrdag,nr12,ntbx,ntx,ord_t,ord_ham: ',
     &             nrdag,nr12,ntbx,ntx,ord_t,ord_ham
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'r12_truncation','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
