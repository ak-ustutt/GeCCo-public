      subroutine select_mrcc_lag(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     selects only those terms of an MRCC lagrangian that meet 3 rules:
*     1.) no self-contractions
*         (if mode='no_tt', contractions involving more than one
*          T operator are also excluded)
*     2.) no pure T-T-contractions
*     only if mode is not 'no_con':
*     3.) connected terms (if mode='pre', contractions involving H
*                          and ops other than T also count as connected)
*
*     matthias, mar 2010
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
     &     idxtop, idxcum, idxham
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
        call write_title(luout,wst_dbg_subr,'select_mrcc_lag')
        write(luout,*) 'mode = ',trim(mode)
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
     &       call quit(1,'select_mrcc_lag','need exactly 3 labels')
        call quit(1,'select_mrcc_lag','Labels not on list!')
      end if

      pure_con = trim(mode).ne.'pre'
      no_tt = trim(mode).eq.'no_tt'
      no_con = trim(mode).eq.'no_con'

      idxham  = idxop(1)
      idxtop  = idxop(2)
      idxcum  = idxop(3)

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

          if (contr%nxarc.gt.0) call quit(1,'select_mrcc_lag',
     &         'not yet adapted for open diagrams')

          ! find out:
          ! - number of T operators
          ! - number of cumulants (and their supervertex number)
          ntop  = 0
          ncum  = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
            end if
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
c          print *,'ntop, ncum: ',ntop,ncum
c          if (ncum.gt.0) print *,'svtx_cum: ',svtx_cum
c dbgend

          if (ncum.gt.0) then
            ! determine cumulant class (0:default;
            !   1: cnt. with Op. other than H,T; 2: cnt. with H)
            ! and whether more than one operator is contracted with cumulant
            allocate(cum_cls(ncum),cum_cnt(ncum))
            cum_cls = 0
            cum_cnt = 0
            do iarc = 1, contr%narc
              vtx1 = contr%arc(iarc)%link(1)
              vtx2 = contr%arc(iarc)%link(2)
              idx_op1 = vertex(vtx1)%idx_op
              idx_op2 = vertex(vtx2)%idx_op
              found = .false.
              if (idx_op1.eq.idxcum) then
                ! no contraction between cumulants
                if (idx_op2.eq.idxcum)
     &             call quit(1,'select_mrcc_lag','contr. between cum.?')
                idx_op  = idx_op2
                svtx1 = contr%svertex(vtx1)
                svtx2 = contr%svertex(vtx2)
                found = .true.
              else if (idx_op2.eq.idxcum) then
                idx_op  = idx_op1
                svtx1 = contr%svertex(vtx2)
                svtx2 = contr%svertex(vtx1)
                found = .true.
              end if
              if (found) then
                icum = idxlist(svtx1,svtx_cum,ncum,1)
                if (cum_cls(icum).eq.0) then
                  if (idx_op.eq.idxham) then
                    cum_cls(icum) = 2
                  else if (idx_op.ne.idxtop) then
                    cum_cls(icum) = 1
                  end if
                else if (cum_cls(icum).eq.1) then
                  if (.not.pure_con.and.idx_op.eq.idxham)
     &                 cum_cls(icum) = 2
                else
                  if (pure_con.and.idx_op.ne.idxham
     &                .and.idx_op.ne.idxtop) cum_cls(icum) = 1
                end if

                if (cum_cnt(icum).eq.0) then
                  cum_cnt(icum) = svtx2
                else if (abs(cum_cnt(icum)).ne.svtx2
     &                   .and.cum_cnt(icum).ne.1000) then
                  if (no_tt) then
                    if (idx_op.eq.idxtop) then
                      if (vertex(idxlist(abs(cum_cnt(icum)),
     &                    contr%svertex,nvtx,1))%idx_op.eq.idxtop) then
                        cum_cnt(icum) = 1000
                      else
                        cum_cnt(icum) = -svtx2
                      end if
                    else
                      cum_cnt(icum) = -abs(cum_cnt(icum))
                    end if
                  else
                    cum_cnt(icum) = -1
                  end if
                end if
              end if
            end do            
          end if

          if (ntop.gt.0) then
            allocate(t_con(ntop))
            t_con = .false.
            itop = 0
            ivtx = 0
            do while(itop.lt.ntop)
              ivtx = ivtx + 1
              if (vertex(ivtx)%idx_op.ne.idxtop) cycle
              itop = itop + 1
              do iarc = 1, contr%narc
                if (t_con(itop)) exit
                vtx1 = contr%arc(iarc)%link(1)
                vtx2 = contr%arc(iarc)%link(2)
                idx_op1 = vertex(vtx1)%idx_op
                idx_op2 = vertex(vtx2)%idx_op
                found = .false.
                if (vtx1.eq.ivtx) then
                  ! no contraction between T op.s
                  if (idx_op2.eq.idxtop)
     &               call quit(1,'select_mrcc_lag','contr. between T?')
                  idx_op  = idx_op2
                  svtx1 = contr%svertex(vtx2)
                  found = .true.
                else if (vtx2.eq.ivtx) then
                  idx_op  = idx_op1
                  svtx1 = contr%svertex(vtx1)
                  found = .true.
                end if
                if (found) then
                  t_con(itop) = t_con(itop).or.idx_op.eq.idxham
                  if (idx_op.eq.idxcum) then
                    icum = idxlist(svtx1,svtx_cum,ncum,1)
                    t_con(itop) = t_con(itop).or.cum_cls(icum).eq.2
                  end if
                end if
              end do
            end do
          end if
c dbg
c          if (ncum.gt.0) print *,'cum_cls: ',cum_cls
c          if (ncum.gt.0) print *,'cum_cnt: ',cum_cnt
c          if (ntop.gt.0) print *,'t_con: ',t_con
c dbgend

          delete = .false.
          if (ncum.gt.0) then
            ! rule 1: no self-contractions
            ! (for no_tt=T: no more than one T in a contraction)
            delete = delete.or.any(cum_cnt.ge.0)
            ! rule 2: no pure T-T contractions
            delete = delete.or.any(cum_cls.eq.0)
            deallocate(svtx_cum,cum_cls,cum_cnt)
          end if
          if (ntop.gt.0) then
            ! rule 3: all T must be connected with H
            delete = delete.or.(any(.not.t_con).and..not.no_con)
            deallocate(t_con)
          end if

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
          call quit(1,'select_mrcc_lag','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
      
      
