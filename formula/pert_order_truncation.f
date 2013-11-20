*----------------------------------------------------------------------*
      subroutine pert_order_truncation(flist,order,idx_tgt,op_info)
*----------------------------------------------------------------------*
*     truncate formula to terms of specific total perturbational order
*     matthias, 2008
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
     &     order, idx_tgt

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx, iord, idx_op, t_max_ord, l_max_ord, op_spec

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'pert_order_trunction')
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
          form_pnt%target = idx_tgt
        case(command_add_contribution)

          form_pnt%target = idx_tgt
          ! assign new result index
          ! comment: operator block should also be changed
          form_pnt%contr%idx_res = idx_tgt

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex


          ! find out: total perturbation order of contraction
          iord = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (op_info%op_arr(idx_op)%op%order.lt.0)
     &        call quit(1,'pert_order_truncation',
     &           'no order assigned to '//
     &           trim(op_info%op_arr(idx_op)%op%name))
            iord = iord + op_info%op_arr(idx_op)%op%order
          end do

          ! order >= 0: delete term if total order does not equal order
          ! order < 0: delete term if zero by (2n+1) and (2n+2) rules
          if (order.ge.0) then
            delete = (iord.ne.order)
          else
            t_max_ord = int((real(iord)-1)/2+0.6)
            l_max_ord = int((real(iord)-2)/2+0.6)
            if (iord.eq.0) l_max_ord = -1
            delete = .false.
            do ivtx = 1, nvtx
              idx_op  = vertex(ivtx)%idx_op
              iord = op_info%op_arr(idx_op)%op%order
              op_spec = op_info%op_arr(idx_op)%op%species
              if ( ( (op_spec.eq.1) .and.
     &                 (iord.gt.t_max_ord) ) .or.
     &                 ( (op_spec.eq.2) .and.
     &                 (iord.gt.l_max_ord) ) ) delete = .true.
            end do
          end if

          if (delete) then
            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'delete_non_fact','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
