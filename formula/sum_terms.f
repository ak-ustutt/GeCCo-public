*----------------------------------------------------------------------*
      subroutine sum_terms(fl_tgt,op_info)
*----------------------------------------------------------------------*
*     look for equal terms and sum them
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_array.h'
      include 'def_formula_item_list.h'
      
      type(formula_item), target, intent(inout) ::
     &     fl_tgt
      ! only for debug output:
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idxop_tgt, iblk_tgt, iterm, jterm, ivtx, iarc,
     &     nvtx, nj, narc, nxarc, nsupvtx, sum_op1, sum_op,
     &     sum_blk1, sum_blk, occ1(ngastp,2), occ(ngastp,2)
      logical ::
     &     unique_set, ok
      integer(8) ::
     &     hash

      type(formula_item), pointer ::
     &     fl_tgt_pnt, fl_tgt_pnt_next, fl_tgt_current
      type(cntr_vtx), pointer ::
     &     vtx1(:), vtx2(:)
      type(cntr_arc), pointer ::
     &     arc1(:), arc2(:)

      integer(8), pointer ::
     &     ivtx1(:),topo1(:,:),xlines1(:,:),
     &     ivtx2(:),topo2(:,:),xlines2(:,:)

      integer, external ::
     &     i8list_cmp, njres_contr
      logical, external ::
     &     iocc_zero

      ! if first entry is an [END]: do nothing
      if (fl_tgt%command.eq.command_end_of_formula) return

      if (ntest.ge.100) then
        write(luout,*) '====================='
        write(luout,*) ' info from sum_terms'
        write(luout,*) '====================='
      end if

      if (fl_tgt%command.ne.command_set_target_init) then
        if (fl_tgt%command.ne.command_add_contribution)
     &       call quit(1,'sum_terms',
     &         'must start with [INIT] or [ADD]')
        idxop_tgt = fl_tgt%target
      end if

      iterm = 0
      fl_tgt_current => fl_tgt
      ! loop over target items
      tgt_loop: do

        ! new operator target ?
        if (fl_tgt_current%command.eq.command_set_target_init) then
          idxop_tgt = fl_tgt_current%target
          if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target)')
          if (ntest.ge.100) then
            write(luout,'(70("="))')
            write(luout,*) 'New operator target: ',idxop_tgt
            write(luout,'(70("="))')
          end if
          fl_tgt_current => fl_tgt_current%next
        end if

        ! prepare some integers and checksums for comparison
        iblk_tgt = fl_tgt_current%contr%iblk_res
        nvtx = fl_tgt_current%contr%nvtx
        nj = njres_contr(fl_tgt_current%contr)
        narc = fl_tgt_current%contr%narc
        nxarc = fl_tgt_current%contr%nxarc
        nsupvtx = fl_tgt_current%contr%nsupvtx
        vtx1 => fl_tgt_current%contr%vertex
        arc1 => fl_tgt_current%contr%arc
        sum_op1  = 0
        sum_blk1 = 0
        do ivtx = 1, nvtx
          sum_op1  = sum_op1  + vtx1(ivtx)%idx_op
          sum_blk1 = sum_blk1 + vtx1(ivtx)%iblk_op
        end do
        occ1 = 0
        do iarc = 1, narc
          occ1 = occ1 + arc1(iarc)%occ_cnt
        end do
        unique_set = fl_tgt_current%contr%unique_set
        if (unique_set) then
          ivtx1 => fl_tgt_current%contr%vtx
          topo1 => fl_tgt_current%contr%topo
          xlines1 => fl_tgt_current%contr%xlines
          hash = fl_tgt_current%contr%hash
        end if

        iterm = iterm+1
        jterm = iterm+1
c dbg
c        print *,'term: ',iterm
c dbgend
        if (ntest.ge.100) then
          write(luout,*) 'current term: # ',iterm
          call prt_contr2(luout,fl_tgt_current%contr,op_info)
        end if

        if (.not.associated(fl_tgt_current%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target)')
        fl_tgt_pnt => fl_tgt_current%next
        search_loop: do
          ! we search within a single domain only
          if (fl_tgt_pnt%command.ne.command_add_contribution)
     &       exit search_loop

          ! as the node might get deleted, better save the next pointer
          if (.not.associated(fl_tgt_pnt%next))
     &         call quit(1,'sum_terms',
     &         'unexpected end of list (target, inner loop)')
          fl_tgt_pnt_next => fl_tgt_pnt%next

          ! pre-screening
          ok = .true.
          if (unique_set.and.fl_tgt_pnt%contr%unique_set) then
            ok = hash.eq.fl_tgt_pnt%contr%hash
c dbg
c            if (ok) print*,jterm,' matches with hash (1): ',hash
c dbgend
          else
            ok = fl_tgt_pnt%contr%iblk_res.eq.iblk_tgt.and.
     &           fl_tgt_pnt%contr%nvtx.eq.nvtx.and.
     &           fl_tgt_pnt%contr%narc.eq.narc.and.
     &           fl_tgt_pnt%contr%nxarc.eq.nxarc.and.
     &           fl_tgt_pnt%contr%nsupvtx.eq.nsupvtx

            if (ok) then
              ! pre-screen with checksums
              vtx2 => fl_tgt_pnt%contr%vertex
              sum_op  = sum_op1
              sum_blk = sum_blk1
              do ivtx = 1, fl_tgt_pnt%contr%nvtx
                sum_op  = sum_op  - vtx2(ivtx)%idx_op
                sum_blk = sum_blk - vtx2(ivtx)%iblk_op
              end do
              ok = sum_op.eq.0.and.sum_blk.eq.0
              if (ok) then
                arc2 => fl_tgt_pnt%contr%arc
                occ = occ1
                do iarc = 1, fl_tgt_pnt%contr%narc
                  occ = occ - arc2(iarc)%occ_cnt
                end do
                ok = iocc_zero(occ)
              end if
            end if
            if (ok) then
              ! set unique vtx/topo/xlines if necessary and compare hash
              if (.not.unique_set) then
                call topo_set_unique(fl_tgt_current%contr)
                unique_set = .true.
                ivtx1 => fl_tgt_current%contr%vtx
                topo1 => fl_tgt_current%contr%topo
                xlines1 => fl_tgt_current%contr%xlines
                hash = fl_tgt_current%contr%hash
              end if
              if (.not.fl_tgt_pnt%contr%unique_set)
     &               call topo_set_unique(fl_tgt_pnt%contr)
              ok = hash.eq.fl_tgt_pnt%contr%hash
c dbg
c              if (ok) print*,jterm,' matches with hash (2): ',hash
c dbgend
            end if
          end if

          ! now compare the unique vtx, topo, xlines
          if (ok) then
            ivtx2 => fl_tgt_pnt%contr%vtx
            ok = i8list_cmp(ivtx1,ivtx2,nvtx).eq.0
            if (ok) then
              topo2 => fl_tgt_pnt%contr%topo
              ok = i8list_cmp(topo1,topo2,nvtx*nvtx).eq.0
              if (ok) then
                xlines2 => fl_tgt_pnt%contr%xlines
                ok = i8list_cmp(xlines1,xlines2,nvtx*nj).eq.0
              end if
            end if
          end if

          ! sum terms
          if (ok) then
c dbg
c            print *,' -> summing!'
c dbgend
            if (ntest.ge.100) then
              write(luout,*) 'found equal term: # ',jterm
              call prt_contr2(luout,fl_tgt_pnt%contr,
     &             op_info)
              write(luout,*) 'now summing and deleting'
            end if
            fl_tgt_current%contr%fac =
     &           fl_tgt_current%contr%fac + fl_tgt_pnt%contr%fac
            call delete_fl_node(fl_tgt_pnt)
            deallocate(fl_tgt_pnt)
          end if

          fl_tgt_pnt => fl_tgt_pnt_next
          jterm = jterm+1
        end do search_loop

        ! advance to next item
        ! end of formula list?
        if (.not.associated(fl_tgt_current%next)) exit tgt_loop
        ! go to next item
        fl_tgt_current => fl_tgt_current%next
        if (fl_tgt_current%command.eq.command_end_of_formula)
     &                                             exit tgt_loop

      end do tgt_loop
c dbg
c      print *,' sum_terms yields',iterm,' terms.'
c dbgend
        
      return
      end
