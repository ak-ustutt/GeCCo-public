      subroutine r12_hhat_truncation(flist,idxham,idxtop,op_info)
*----------------------------------------------------------------------*
*     truncate formula for Hhat
*
*     matthias, 2009
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
     &     idxtop, idxham

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op,
     &     ntop, ntx, nham, ord_ham

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

c      integer, external ::
c     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'Hhat_truncation')
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
          ! - number of R+
          ! - number of R
          ntop  = 0
          ntx   = 0
          nham    = 0
          ord_ham = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              if (op_info%op_arr(idx_op)%
     &            op%ihpvca_occ(IEXTR,1,iblk_op).gt.0)
     &             ntx = ntx+1
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              ord_ham = ord_ham
     &           + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
     &           + op_info%op_arr(idx_op)%op%ihpvca_occ(IEXTR,1,iblk_op)
     &           + op_info%op_arr(idx_op)%op%ihpvca_occ(IEXTR,2,iblk_op)
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'r12_hhat_truncation','strange: nham.ne.1')
          if (ntop-ntx.lt.0)
     &         call quit(1,'r12_hhat_truncation',
     &                     'strange: ntop-ntx.lt.0')

          ! linear in T1X and no "formal" blocks of hamiltonian
          delete = ntx.gt.1.or.ord_ham.ge.3

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
              write(lulog,*) 'ntx,ord_ham: ',ntx,ord_ham
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'r12_hhat_truncation','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
