      subroutine t1x_truncation(flist,mode,
     &     idxtbar,idxham,idxtop,op_info)
*----------------------------------------------------------------------*
*     preliminary start-up for truncated CC expansions with 
*     T1 into CABS
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
     &     delete, diag, omit_fpx, special
      integer ::
     &     nvtx, ivtx,
     &     idx_op, iblk_op, iblk_t1x, iblk_l1x, iblk_t1,
     &     nt1, nt1x, nl1x, nham, nblk_h, iblk,
     &     ord_ham, n_ext, rank, xrank, h0_def,
     &     max_pert, max_comm, max_t1x
      character*64 ::
     &     op_name
      logical ::
     &     dagger

      integer ::
     &     occ(ngastp,2)

      integer, pointer ::
     &     po_h(:)

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next
      type(operator), pointer ::
     &     op_h

      integer, external ::
     &     rank_occ, iblk_occ
      logical, external ::
     &     occ_is_diag_blk

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'t1x_trunction')
        write(lulog,*) 'mode = ',trim(mode)
      endif

      ! experimental version for CCS project
      special = .false.
      select case(mode(1:4))
      case('ord1')
        special = .true.
        max_pert = 2
        max_t1x  = 40
      case('ord2')
        max_pert = 2
        max_t1x  = 40
      case('ord3')
        max_pert = 3
        max_t1x  = 40
      case('ord4')
        max_pert = 4
        max_t1x  = 40
      case('ord5')
        max_pert = 5
        max_t1x  = 40
      case('ord6')
        max_pert = 6
        max_t1x  = 40
      case('ord7')
        max_pert = 7
        max_t1x  = 40
      case('ord8')
        max_pert = 8
        max_t1x  = 40
      case('line')
        max_pert = 5
        max_t1x  = 1
      case('quad')
        max_pert = 5
        max_t1x  = 2
      case default
        call quit(1,'t1x_truncation','what do you mean: '//trim(mode))
      end select
      h0_def = 0
      if (mode(6:6).eq.'1') h0_def=1
      if (mode(6:6).eq.'2') h0_def=2
      omit_fpx = (mode(7:7).eq.'d') 

      ! obtain blocks of T1, T1x, and L1x
      occ = 0
      occ(IPART,1) = 1
      occ(IHOLE,2) = 1
      iblk_t1  = iblk_occ(occ,.false.,op_info%op_arr(idxtop)%op,1)
      occ = 0
      occ(IEXTR,1) = 1
      occ(IHOLE,2) = 1
      iblk_t1x = iblk_occ(occ,.false.,op_info%op_arr(idxtop)%op,1)
      iblk_l1x = iblk_occ(occ,.true., op_info%op_arr(idxtbar)%op,1)

      ! set up perturbation order of Hamiltonian elements
      op_h => op_info%op_arr(idxham)%op
      nblk_h = op_h%n_occ_cls
      allocate(po_h(nblk_h))

      do iblk = 1, nblk_h
        n_ext = sum(op_h%ihpvca_occ(IEXTR,1:2,iblk))
        rank  = rank_occ('C',op_h%ihpvca_occ(1:,1:,iblk),1)
        xrank  = rank_occ('X',op_h%ihpvca_occ(1:,1:,iblk),1)
        diag  = occ_is_diag_blk(op_h%ihpvca_occ(1:,1:,iblk),1)
        po_h(iblk) = 0
        if (n_ext.gt.0) po_h(iblk) = 1
        if (h0_def.eq.0.and.rank.gt.1) po_h(iblk) = 1
        if (h0_def.eq.2.and.rank.gt.1) po_h(iblk) = po_h(iblk)+1
        if (rank.eq.1.and.diag) po_h(iblk) = 0
        if (.not.omit_fpx.and.rank.eq.1.and.xrank.eq.0) po_h(iblk) = 0
c fusk
c        call warn('t1x_truncation','PFUSCH!!!!')
c        if (rank.eq.1) po_h(iblk) = 0
c fusk
c test
        if (special) then
          po_h(iblk) = 3
          if (rank.le.1) po_h(iblk) = 0          
        end if
c test
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'PO(H) table:'
        do iblk = 1, nblk_h
          write(lulog,*) iblk,po_h(iblk)
        end do
      end if

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
          nt1   = 0
          nt1x  = 0
          nl1x  = 0
          nham    = 0
          ord_ham = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop.and.iblk_op.eq.iblk_t1 ) nt1  = nt1 +1
            if (idx_op.eq.idxtop.and.iblk_op.eq.iblk_t1x) nt1x = nt1x+1
            if (idx_op.eq.idxtbar.and.iblk_op.eq.iblk_l1x) nl1x = nl1x+1
            if (idx_op.eq.idxham) then
              nham = nham+1
              ord_ham = ord_ham + po_h(iblk_op)
            end if
          end do

          if (nham.ne.1)
     &         call quit(1,'t1x_truncation','strange: nham.ne.1')
          ! restrict order
          delete = (ord_ham+nt1x+nl1x).gt.max_pert
          ! restrict max. T1
          delete = delete.or.(nt1x.gt.max_t1x)
          ! special:
          if (special) then
            delete = delete.or.nt1+nt1x.gt.1
          end if

          if (delete) then
            ! Print the deleted contraction.
            if (ntest.ge.1000) then
              write(lulog,*) 'ord_ham,nt1x,nl1x: ',ord_ham,nt1x,nl1x
              write(lulog,*) 'max_pert = ',max_pert
              write(lulog,*) 'max_t1x = ',max_t1x
              write(lulog,*)'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'t1x_truncation','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      end do

      deallocate(po_h)

      return
      end
      
      
