      subroutine select_mrcc_lag3(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     delete terms where the first operator (cumulant) is only
*     connected to the second operator (lambda)
*     (becomes disconnected upon differentiation w.r.t. lambda)
*
*     matthias, april 2012
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
     &     delete, error, pure_con, found, no_tt, no_con
      integer ::
     &     idxtop, idxcum, idxham, idxl
      integer ::
     &     ii, idx_op, ivtx, nvtx, iarc, vtx1, vtx2,
     &     ntop, ncum, svtx1, svtx2, idx_op1, idx_op2, icum, itop
      integer ::
     &     idxop(nlabels)
      integer, allocatable ::
     &     svtx_cum(:), tmp(:), cum_cls(:), cum_cnt(:)
      logical, allocatable ::
     &     t_con(:)


      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, idxlist

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_mrcc_lag3')
c        write(lulog,*) 'mode = ',trim(mode)
      endif

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.2
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.2)
     &       call quit(1,'select_mrcc_lag3','need exactly 2 labels')
        call quit(1,'select_mrcc_lag3','Labels not on list!')
      end if

      idxcum  = idxop(1)
      idxl    = idxop(2)

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

          contr => form_pnt%contr
          nvtx = contr%nvtx
          vertex => contr%vertex

          if (contr%nxarc.gt.0) call quit(1,'select_mrcc_lag3',
     &         'not yet adapted for open diagrams')

          ! find out:
          ! - number of cumulants (and their supervertex number)
          ncum  = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxcum) then
              svtx1 = contr%svertex(ivtx)
              if (ncum.eq.0) then
                ncum = ncum + 1
                allocate(svtx_cum(ncum))
                svtx_cum(ncum) = svtx1
              else if (idxlist(svtx1,svtx_cum,ncum,1)
     &                 .le.0) then
                allocate(tmp(ncum))
                tmp = svtx_cum
                deallocate(svtx_cum,STAT=ii)
                ncum = ncum+1
                allocate(svtx_cum(ncum))
                svtx_cum(1:ncum-1) = tmp(1:ncum-1)
                deallocate(tmp)
                svtx_cum(ncum) = svtx1
              end if
            end if
          end do
c dbg
c          print *,'ncum: ',ncum
c          if (ncum.gt.0) print *,'svtx_cum: ',svtx_cum
c dbgend

          if (ncum.gt.0) then
            ! determine cumulant class (0:default;
            !   1: cnt. with Op. other than L)
            allocate(cum_cls(ncum))
            cum_cls = 0
            do iarc = 1, contr%narc
              vtx1 = contr%arc(iarc)%link(1)
              vtx2 = contr%arc(iarc)%link(2)
              idx_op1 = vertex(vtx1)%idx_op
              idx_op2 = vertex(vtx2)%idx_op
              found = .false.
              if (idx_op1.eq.idxcum) then
                ! no contraction between cumulants
                if (idx_op2.eq.idxcum)
     &            call quit(1,'select_mrcc_lag3','contr. between cum.?')
                idx_op  = idx_op2
                svtx1 = contr%svertex(vtx1)
                found = .true.
              else if (idx_op2.eq.idxcum) then
                idx_op  = idx_op1
                svtx1 = contr%svertex(vtx2)
                found = .true.
              end if
              if (found.and.idx_op.ne.idxl) then
                icum = idxlist(svtx1,svtx_cum,ncum,1)
                cum_cls(icum) = 1
              end if
            end do            
          end if

          delete = .false.
          if (ncum.gt.0) then
            delete = delete.or.any(cum_cls.eq.0)
            deallocate(svtx_cum,cum_cls)
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_mrcc_lag3','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
