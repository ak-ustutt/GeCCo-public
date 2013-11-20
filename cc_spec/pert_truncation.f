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
     &     ntop, nt1, nham, ntbar, ntmax,
     &     ord_t, ord_this_t, ord_tbar, ord_ham,
     &     max_pert, max_comm, max_t1, tb_trunc
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

      ! experimental version for CC2 project
      select case(trim(mode))
      case('CC2')
        tb_trunc = 1
        max_pert = 2
        max_comm = 4
        max_t1   = 4
      case('CC2-l')
        tb_trunc = 1
        max_pert = 2
        max_comm = 1
        max_t1   = 4
      case('CC2-t1l')
        tb_trunc = 1
        max_pert = 2
        max_comm = 4
        max_t1   = 1
      case('CC2-q')
        tb_trunc = 1
        max_pert = 2
        max_comm = 2
        max_t1   = 4
      case('CC2-c')
        tb_trunc = 1
        max_pert = 2
        max_comm = 3
        max_t1   = 4
      case('CC3')
        tb_trunc = 2
        max_pert = 4
        max_comm = 4
        max_t1   = 4
      case('CCSDT-3')
        tb_trunc = 2
        max_pert = 100
        max_comm = 4
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
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)

          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! find out:
          ! - number and perturbation order of T operators
          ! - perturbation order of H
          ! - perturbation order of TBAR
          ntop  = 0
          ntmax = 0
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
              ord_this_t= op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
              ord_t = ord_t + ord_this_t
              if (ord_this_t.eq.0)
     &             nt1 = nt1+1
              if (ord_this_t.eq.tb_trunc)
     &             ntmax = ntmax+1
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
c          if (ntbar.gt.1)
c     &         call quit(1,'pert_truncation','strange: ntbar.ne.1')
          ! restrict to second order (T1 counts 0 here)
          delete =  ord_tbar.eq.tb_trunc.and.
     &             (ord_ham+ord_t+ord_tbar).gt.max_pert
C this one was syntactically incorrect; does it still work?? 
C          delete = (ord_tbar.eq.tb_trunc.and.
C     &              ord_ham+ord_t+ord_tbar).gt.max_pert
          ! avoid <0|TBARmax [[F,T1],Tmax]|0>
          delete = delete.or.
C     &         (ord_ham.eq.0.and.ord_tbar.gt.0.and.nt1.gt.0)
     &        (ord_tbar.eq.tb_trunc.and.
     &      ord_ham.eq.0.and.ntmax.gt.0.and.nt1.gt.0)
          ! avoid <0|TBARmax [H,Tmax] |0> and higher (CCSDT-3 etc.)
          delete = delete.or.
     &        (ord_tbar.eq.tb_trunc.and.ord_ham.ne.0.and.ntmax.gt.0)
          ! restrict max. commutators
          delete = delete.or.
     &         (ntop.gt.max_comm)
          ! restrict max. T1
          delete = delete.or.
     &         (nt1.gt.max_t1)

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*)'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

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
      
      
