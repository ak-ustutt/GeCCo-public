      subroutine select_same_blk(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     given two operators, checks whether they appear with the same
*     occupation class number
*
*     matthias, dec 2010
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
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
     &     ii, idx_op, ivtx, nvtx, idx_op1, idx_op2, iblk
      integer ::
     &     idxop(nlabels)

      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_same_blk')
        write(luout,*) 'mode = ',trim(mode)
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.2
      
      if (error) then
        write(luout,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(luout,'(a20," - ??")') trim(labels(ii))
          else
            write(luout,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.2)
     &       call quit(1,'select_same_blk','need exactly 2 labels')
        call quit(1,'select_same_blk','Labels not on list!')
      end if

      idx_op1  = idxop(1)
      idx_op2  = idxop(2)

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

          contr => form_pnt%contr
          nvtx = contr%nvtx
          vertex => contr%vertex

          ! find out if block numbers are the same:
          iblk = 0
          delete = .false.
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idx_op1.or.idx_op.eq.idx_op2) then
              if (iblk.eq.0) then
                iblk = vertex(ivtx)%iblk_op
                cycle
              else
                delete = iblk.ne.vertex(ivtx)%iblk_op
              end if
            end if
          end do

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*) 'Deleted formula item:'
              call prt_contr2(luout,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'select_same_blk','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
