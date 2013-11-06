      subroutine select_rank(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     only keeps terms in which either (!!!) of the given operators
*     appears with the rank specified in the single-digit
*     integer string "mode".
*     Terms without these operators are also deleted.
*
*     matthias, oct 2013
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
     &     idxop(nlabels), irank(nlabels)

      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_rank')
        write(luout,*) 'mode = ',trim(mode)
      endif

      ! get operator indices and desired ranks
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        read(mode(ii:ii),'(i1)') irank(ii)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.le.0
      
      if (error) then
        write(luout,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(luout,'(a20," - ??")') trim(labels(ii))
          else
            write(luout,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.le.0)
     &       call quit(1,'select_rank','kidding me?')
        call quit(1,'select_rank','Labels not on list!')
      end if

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

          ! find out if any of the operators appears with desired rank
          delete = .true.
          do ii = 1, nlabels
             idx_op = idxop(ii)
            ! does the result operator match?
            if (idx_op.eq.contr%idx_res) then
              iblk = contr%iblk_res
              if (op_info%op_arr(idx_op)%op%ica_occ(1,iblk)
     &            .eq.irank(ii)) then
                delete = .false.
                exit
              end if
            end if
            ! now check the other operators
            do ivtx = 1, nvtx
              if (idx_op.eq.vertex(ivtx)%idx_op) then
                iblk = vertex(ivtx)%iblk_op
                if (op_info%op_arr(idx_op)%op%ica_occ(1,iblk)
     &              .eq.irank(ii)) then
                  delete = .false.
                  exit
                end if
              end if
            end do
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
          call quit(1,'select_rank','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
