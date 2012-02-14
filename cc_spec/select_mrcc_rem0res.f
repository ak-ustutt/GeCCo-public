      subroutine select_mrcc_rem0res(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     remove residuals for fixed part of T that are assumed to be zero:
*     keep only those terms belonging to such residuals that have
*     at least one T block which is not fixed.
*
*     matthias, nov. 2011
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_input.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     nlabels
      character(len=*), intent(in) ::
     &     labels(nlabels), mode

      logical ::
     &     delete, error
      integer ::
     &     idxtop, idxtfix, iblk, jblk, idx_op,
     &     ivtx, ii, nvtx, iterm, idxl, iblkfix,
     &     jvtx, jblkfix, jdx_op
      integer ::
     &     idxop(nlabels), occ(ngastp,2)

      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(operator), pointer ::
     &     op_l, op_top, op_tfix
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, iblk_occ

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_mrcc_rem0res')
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.3
      
      if (error) then
        write(luout,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(luout,'(a20," - ??")') trim(labels(ii))
          else
            write(luout,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.3)
     &       call quit(1,'select_mrcc_rem0res','wrong number of labels')
        call quit(1,'select_mrcc_rem0res','Labels not on list!')
      end if

      idxl = idxop(1)
      idxtop = idxop(2)
      idxtfix = idxop(3)
      op_l => op_info%op_arr(idxl)%op
      op_top => op_info%op_arr(idxtop)%op
      op_tfix => op_info%op_arr(idxtfix)%op
      iterm = 0

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

          iterm = iterm + 1
c dbg
c          print *,'current term: ',iterm
c dbgend
          contr => form_pnt%contr
          nvtx = contr%nvtx
          vertex => contr%vertex

          delete = .false.

          ! search terms belonging to redundant residuals
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxl) then
              iblk = vertex(ivtx)%iblk_op
              ! is there an equivalent block in Tfix?
              occ(1:ngastp,1:2) = op_l%ihpvca_occ(1:ngastp,1:2,iblk)
              iblkfix = iblk_occ(occ,.true.,op_tfix,
     &                           op_l%blk_version(iblk))
              if (iblkfix.gt.0) then
                delete = .true.
                ! keep only if any T is not part of Tfix
                do jvtx = 1, nvtx
                  jdx_op = vertex(jvtx)%idx_op
                  if (jdx_op.eq.idxtop) then
                    jblk = vertex(jvtx)%iblk_op
                    ! is there an equivalent block in Tfix?
                    occ(1:ngastp,1:2) 
     &                        = op_top%ihpvca_occ(1:ngastp,1:2,jblk)
                    jblkfix = iblk_occ(occ,.false.,op_tfix,
     &                                 op_top%blk_version(jblk))
                    delete = delete.and.jblkfix.gt.0
                    if (.not.delete) exit
                  end if
                end do
              end if
            end if
          end do

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*) 'Deleted formula item:'
              write(luout,*) 'iterm:',iterm
              call prt_contr2(luout,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'select_mrcc_rem0res','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
