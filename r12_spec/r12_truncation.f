      subroutine r12_truncation(flist,trunc_type,trunc_t1x,
     &     idxr12,idxham,idxtbar,idxtop,idxcbar,idxc12,op_info)
*----------------------------------------------------------------------*
*     preliminary start-up for truncated R12 expansions
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
     &     trunc_type, idxr12, idxtbar, idxtop, idxham, idxc12, idxcbar,
     &     trunc_t1x

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op,
     &     ntop, ntx, nrdag, nr12, nham, ntbar, ntbx, nc12, ncbar,
     &     ord_t, ord_ham, ord_c12, ord_cbar,
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
          ! - number of R+
          ! - number of R
          ! - T (TBAR) p connected to R (R+): count as 
          !   special T (TBAR) operator
          ntop  = 0
          ntx   = 0
          ntbar = 0
          ntbx  = 0
          nrdag = 0
          nr12  = 0
          ncbar = 0
          nc12  = 0
          nham    = 0
          ord_ham = 0
          ord_t   = 0
          ord_c12 = 0
          ord_cbar= 0
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
            if (idx_op.eq.idxcbar) then
              ncbar = ncbar+1
              ! no need to add (there should be only one)
              ord_cbar = op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
            if (idx_op.eq.idxc12) then
              nc12 = nc12+1
              ord_c12 = ord_c12+
     &             op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'r12_truncation','strange: nham.ne.1')
          if (ntop-ntx.lt.0)
     &         call quit(1,'r12_truncation','strange: ntop-ntx.lt.0')

          ! well, old behaviour only for trunc_t1x.eq.-1
          ! else, we forget about these operators and treat them later:
          if (trunc_t1x.ne.-1) ntx = 0
          if (trunc_t1x.ne.-1) ntbx = 0

          ! always linear in T1X
          delete = ntx.gt.1

          if (trunc_type.ne.1) then
            ! (R12):
            ! always linear in R12
            delete = delete.or.nr12+ntx.gt.1
            ! R12 projection:
            delete = delete.or.nrdag+ntbx.gt.0.and.
     &           ord_ham+(ntop-ntx).gt.0.and.nr12+ntx.gt.0
     &          .and.(ntop-ntx.gt.0 .or. ord_cbar.eq.ord_c12)            
          else if (trunc_type.eq.1) then
            ! linearized R12:
            delete = delete.or.nr12.gt.1
          end if
          if (trunc_type.eq.2) then
            ! no R12 at all:
            delete = delete.or.nr12.gt.0
            delete = delete.or.nrdag.gt.0
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
              write(lulog,*) 'nrdag,nr12,ntbx,ntx,ord_t,ord_ham: ',
     &             nrdag,nr12,ntbx,ntx,ord_t,ord_ham
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'r12_truncation','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
