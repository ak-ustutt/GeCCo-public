      subroutine r12exc_split_terms(flist,
     &     mode,idxop_res,
     &     idx_ham,idx_top,idx_r12,idx_l,idx_r,idx_v,idx_tprime,
     &     op_info)
*----------------------------------------------------------------------*
*     split r12 excitation energy for analysis purposes
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

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     mode,
     &     idxop_res,
     &     idx_ham, idx_r12, idx_top, idx_l, idx_r, idx_v, idx_tprime

      logical ::
     &     delete, is_typ1, is_typ2, is_typ3,
     &     is_r1l1, has_corr, has_vhh
      integer ::
     &     nvtx, ivtx, narc, iarc,
     &     idx_op, iblk_op, idx,
     &     nham, nt_gt_2, nr12, nr12ph, nr12phdg,
     &     nlx, nrx, rank_rx, rank_lx, nother, np_other, rank
      character*64 ::
     &     op_name
      logical ::
     &     found

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     rank_occ

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_terms')
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
          form_pnt%target = idxop_res
        case(command_add_contribution)

          if (ntest.ge.1000) then
            write(lulog,*) 'current item:'
            call prt_contr2(lulog,form_pnt%contr,op_info)
          end if

          nvtx = form_pnt%contr%nvtx
          narc = form_pnt%contr%narc
          vertex => form_pnt%contr%vertex
          arc => form_pnt%contr%arc

          ! count: # of H
          !        # of non-singles T
          !        # of R12
          !        # of non-hh R12 or R12^+
          !        block of L
          !        block of R
          has_vhh = .false.
          nham = 0
          nt_gt_2 = 0
          nr12 = 0
          nr12ph = 0
          nr12phdg = 0
          nlx = 0
          rank_rx = 0
          nrx = 0
          rank_lx = 0
          nother = 0
          np_other = 0
          do ivtx = 1, nvtx
            idx_op = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op

            if (idx_op.eq.idx_ham) then
              nham = nham+1
            else if (idx_op.eq.idx_top) then
              rank = rank_occ('X',op_info%op_arr(idx_top)%op%
     &             ihpvca_occ(1:,1:,iblk_op),1)
              if (rank.gt.1) nt_gt_2 = nt_gt_2+1
            else if (idx_op.eq.idx_tprime) then
              nt_gt_2 = nt_gt_2+1
            else if (idx_op.eq.idx_r12) then
              nr12 = nr12+1
              if (op_info%op_arr(idx_r12)%op%
     &             ihpvca_occ(2,IPART,iblk_op).gt.0) then
                if (.not.vertex(ivtx)%dagger) nr12ph = nr12ph+1
                if (     vertex(ivtx)%dagger) nr12phdg = nr12phdg+1
              end if
            else if (idx_op.eq.idx_l) then
              nlx = nlx+1
              rank_lx = -rank_occ('X',op_info%op_arr(idx_l)%op%
     &             ihpvca_occ(1:,1:,iblk_op),1)
            else if (idx_op.eq.idx_r) then
              nrx = nrx+1
              rank_rx = rank_occ('X',op_info%op_arr(idx_r)%op%
     &             ihpvca_occ(1:,1:,iblk_op),1)
              if (rank_rx.eq.1) then
                do iarc = 1, narc
                  if (arc(iarc)%link(2).eq.ivtx .and.
     &                arc(iarc)%occ_cnt(IPART,2).eq.1) then 
                    if (vertex(arc(iarc)%link(1))%idx_op.eq.idx_r12.or.
     &                  vertex(arc(iarc)%link(1))%idx_op.eq.idx_v)
     &                   rank_rx = rank_rx+1
                  end if
                end do
              end if
            else if (idx_op.eq.idx_v) then
c              nother = nother+1
              has_vhh =
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IPART,1,iblk_op).eq.0 .and.
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IPART,1,iblk_op).eq.0 .and.
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IHOLE,2,iblk_op).eq.1 .and.
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IHOLE,2,iblk_op).eq.1 
            else
              nother = nother+1
              np_other = np_other +
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IPART,1,iblk_op) +
     &             op_info%op_arr(idx_op)%op%
     &             ihpvca_occ(IPART,2,iblk_op)
            end if
          end do

          is_r1l1 = rank_lx.eq.1.and.rank_rx.eq.1

          has_corr = nham.ne.1.or.nt_gt_2.gt.0.or.nr12.gt.0.or.has_vhh

          if (idx_tprime.le.0) then
            is_typ1 = is_r1l1.and..not.has_corr.and.nother.eq.0
            is_typ2 = is_r1l1.and..not.is_typ1
     &         .and.nr12phdg.eq.0.and.
c     &         .and.nr12ph.eq.0.and.nr12phdg.eq.0.and.
     &         (nother.eq.0.or.has_vhh)
c     &         (nother.eq.0.or.(nother.eq.1.and.has_vhh))
          else
            is_typ1 = is_r1l1.and..not.has_corr.and.nother.eq.0
            is_typ2 = is_r1l1.and..not.is_typ1
          end if
          is_typ3 = .not.is_typ1.and..not.is_typ2

          if (ntest.ge.1000) then
            write(lulog,*) 'rank_lx, rank_rx, nt_gt_2, nham:',
     &                      rank_lx, rank_rx, nt_gt_2, nham
            write(lulog,*) 'nr12,nr12ph,nr12phdg,np_other:',
     &                      nr12,nr12ph,nr12phdg,np_other
            write(lulog,*) 'nother,has_vhh:',
     &                      nother,has_vhh
            write(lulog,*) '-> is_r1l1,has_corr: ',
     &                         is_r1l1,has_corr  
            write(lulog,*) '-> is_typ1, is_typ2, is_typ3: ',
     &                         is_typ1, is_typ2, is_typ3
          end if

          select case(mode)
          case(1)
            delete = .not.is_typ1
          case(2)
            delete = .not.is_typ2
          case(3)
            delete = .not.is_typ3
          end select

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*)'---- Deleted formula item ----'
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else
            form_pnt%target = idxop_res
            form_pnt%contr%idx_res =
     &           idxop_res
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'r12exc_split_terms','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
