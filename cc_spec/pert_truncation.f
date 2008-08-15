      subroutine pert_truncation(flist,mode,
     &     idxtbar,idxham,idxtop,op_info)
*----------------------------------------------------------------------*
*     preliminary start-up for truncated CC expansions
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

      character(len=*) ::
     &     mode
      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     idxtbar, idxtop, idxham

      logical ::
     &     delete
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op,
     &     ntop, nt1, nham, ntbar,
     &     ord_t, ord_tbar, ord_ham,
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

      ! experimental version for CC2 project
      select case(trim(mode))
      case('CC2')
        max_pert = 2
        max_comm = 4
        max_t1   = 4
      case('CC2-l')
        max_pert = 2
        max_comm = 1
        max_t1   = 4
      case('CC2-t1l')
        max_pert = 2
        max_comm = 4
        max_t1   = 1
      case('CC2-q')
        max_pert = 2
        max_comm = 2
        max_t1   = 4
      case('CC2-c')
        max_pert = 2
        max_comm = 3
        max_t1   = 4
      case default
        call quit(1,'pert_truncation','what do you mean: '//trim(mode))
      end select

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
          ! - number and perturbation order of T operators
          ! - perturbation order of H
          ! - perturbation order of TBAR
          ntop  = 0
          nt1   = 0
          ord_t = 0
          ntbar = 0
          ord_tbar = 0
          nham    = 0
          ord_ham = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              ord_t = ord_t
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
              if (op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op).eq.1)
     &             nt1 = nt1+1
            end if
            if (idx_op.eq.idxtbar) then
              ntbar = ntbar+1
              ord_tbar = ord_tbar
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              ord_ham = ord_ham
     &              + op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'pert_truncation','strange: nham.ne.1')
          if (ntbar.gt.1)
     &         call quit(1,'pert_truncation','strange: ntbar.ne.1')
          ! restrict to second order (T1 counts 0 here)
          delete = (ord_ham+ord_t+ord_tbar).gt.max_pert
          ! avoid <0|TBAR2 [[F,T1],T2]|0>
          delete = delete.or.
     &         (ord_ham.eq.0.and.ord_tbar.gt.0.and.nt1.gt.0)
          ! restrict max. commutators
          delete = delete.or.
     &         (ntop.gt.max_comm)
          ! restrict max. T1
          delete = delete.or.
     &         (nt1.gt.max_t1)

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*)'Deleted formula item:'
              call prt_contr2(luout,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'delete_non_fact','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
