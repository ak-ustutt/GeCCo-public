      subroutine select_f12x(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     special selection routine for Stuttgart F12x methods
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
!      include 'def_orbinf.h'
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
     &     idxtbar, idxtop, idxr12, idxham
      integer ::
     &     trunc_e, trunc_t1, trunc_t2, ii,
     &     idx_op, iblk_op, ivtx, nvtx, 
     &     ntop, typ_tbar, ord_ham, nr12, nrdag
      integer ::
     &     idxop(nlabels)


      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_f12x')
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.4
      
      if (error) then
        write(luout,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(luout,'(a20," - ??")') trim(labels(ii))
          else
            write(luout,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.4)
     &       call quit(1,'select_f12x','need exactly 4 labels')
        call quit(1,'select_f12x','Labels not on list!')
      end if

      select case(trim(mode))
      case('F12a','f12a')
        trunc_e = 1
        trunc_t1 = 0
        trunc_t2 = 1
      case('F12b','f12b')
        trunc_e = 999
        trunc_t1 = 0
        trunc_t2 = 1
      case('F12b-S','f12b-S','f12b-s')
        trunc_e = 999
        trunc_t1 = 1
        trunc_t2 = 1
      case default
        call quit(1,'select_f12x','mode?? "'//trim(mode)//'"')
      end select

      idxtbar = idxop(1)
      idxham  = idxop(2)
      idxtop  = idxop(3)
      idxr12  = idxop(4)

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
          ! - number of T operators
          ! - number of R+
          ! - number of R
          ! - type of Tbar
          ntop  = 0
          typ_tbar = 0
          nrdag = 0
          nr12  = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
            end if
            if (idx_op.eq.idxtbar) then
              typ_tbar = op_info%op_arr(idx_op)%
     &             op%ihpvca_occ(IHOLE,1,iblk_op)
            end if
            if (idx_op.eq.idxr12) then
              if (vertex(ivtx)%dagger) then
                nrdag = nrdag+1
              else
                nr12 = nr12+1
              end if
              
            end if
            if (idx_op.eq.idxham) then
              ord_ham = 
     &             op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
            end if
          end do

          if (typ_tbar.eq.0.and.nrdag.eq.0) then
            delete = .false.
          else if (nrdag.eq.1) then
            delete = ntop+nr12+ord_ham.gt.trunc_e
          else if (typ_tbar.eq.1) then
            delete = nr12.gt.trunc_t1
            delete = delete.or.(nr12.gt.0.and.ntop.gt.0)
          else if (typ_tbar.eq.2) then
            delete = nr12.gt.trunc_t2
            delete = delete.or.(nr12.gt.0.and.ntop.gt.0)
          else
            delete = .true.
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(luout,*) 'Deleted formula item:'
              call prt_contr2(luout,form_pnt%contr,op_info)
              write(luout,*) 'nrdag,nr12,ntop,typ_tbar,ord_ham: ',
     &             nrdag,nr12,ntop,typ_tbar,ord_ham
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(luout,*)'command = ',form_pnt%command
          call quit(1,'select_f12x','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo


      return
      end
      
      
