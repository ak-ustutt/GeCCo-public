      subroutine select_mrcc_lag2(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*
*     matthias, Nov. 2010
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

      integer, parameter ::
     &     maxtt = 12, maxtop = 8

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     nlabels
      character(len=*), intent(in) ::
     &     labels(nlabels), mode

      logical ::
     &     delete, error, deldue2maxtt
      integer ::
     &     idxtop, idxham
      integer ::
     &     ii, idx_op, ivtx, nvtx, iarc, vtx1, vtx2,
     &     ntt, ntop, nham, idx_op1, idx_op2, itop,
     &     ntesting, itesting, maxcon_tt, ntt_save
      integer ::
     &     idxop(nlabels), bins(maxtt+1,maxtop+1), binsum(maxtop+1)
      integer, allocatable ::
     &     testing(:)
      logical, allocatable ::
     &     connected(:), bchpart(:)


      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_mrcc_lag2')
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
     &       call quit(1,'select_mrcc_lag2','need exactly 2 labels')
        call quit(1,'select_mrcc_lag2','Labels not on list!')
      end if

      call get_argument_value('method.MRCC','maxtt',
     &     ival=maxcon_tt)

      idxham  = idxop(1)
      idxtop  = idxop(2)

      bins = 0
      deldue2maxtt = .false.

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

          allocate(bchpart(nvtx),connected(nvtx),testing(nvtx))
          bchpart = .false.
          connected = .false.
          testing = 0
          ntesting = 0

          ! find out:
          ! - number of T operators
          ! - vertex number of Hamiltonian
          ntop  = 0
          nham  = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              bchpart(ivtx) = .true.
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              bchpart(ivtx) = .true.
              ntesting = ntesting + 1
              testing(ntesting) = ivtx
              connected(ivtx) = .true.
            end if
          end do
          if (nham.gt.1) call warn('select_mrcc_lag2',
     &       'More than one Hamiltonian? Are we completely crazy now?!')

          ! number of T-T contractions
          ntt = 0
          do iarc = 1, contr%narc
            vtx1 = contr%arc(iarc)%link(1)
            vtx2 = contr%arc(iarc)%link(2)
            idx_op1 = vertex(vtx1)%idx_op
            idx_op2 = vertex(vtx2)%idx_op
            if (idx_op1.eq.idxtop.and.idx_op2.eq.idxtop) ntt = ntt + 1
          end do            

          ! increment binning
          ntt_save = ntt
          ntt = min(ntt,maxtt)
          ntop = min(ntop,maxtop)
          bins(ntt+1,ntop+1) = bins(ntt+1,ntop+1) + 1

          ! are all Ts and H connected?
          itesting = 0
          if (nham.gt.0) then
            do while(itesting.lt.ntesting)
              itesting = itesting + 1
              ivtx = testing(itesting)
              do iarc = 1, contr%narc
                if (ivtx.eq.contr%arc(iarc)%link(1)) then
                  vtx1 = contr%arc(iarc)%link(1)
                  vtx2 = contr%arc(iarc)%link(2)
                else if (ivtx.eq.contr%arc(iarc)%link(2)) then
                  vtx1 = contr%arc(iarc)%link(2)
                  vtx2 = contr%arc(iarc)%link(1)
                else
                  cycle
                end if
                if (.not.bchpart(vtx2)) cycle
                if (connected(vtx2)) cycle
                connected(vtx2) = .true.
                ntesting = ntesting + 1
                testing(ntesting) = vtx2
              end do
              if (ntesting.eq.ntop+nham) exit ! all vtxs connected
            end do
          end if

          delete = ntesting.ne.ntop+nham
          if (delete) print *,'Deleting disconnected term!'

          deallocate(bchpart,connected,testing)

          ! delete if more T-T connections than requested
          if (maxcon_tt.ge.0) then
            delete = delete.or.ntt_save.gt.maxcon_tt
            deldue2maxtt = deldue2maxtt.or.ntt_save.gt.maxcon_tt
          end if

          if (delete) then
            ! Undo binning count
            bins(ntt+1,ntop+1) = bins(ntt+1,ntop+1) - 1
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
          call quit(1,'select_mrcc_lag2','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      do itop = 0,maxtop
        binsum(itop+1) = sum(bins(1:maxtt+1,itop+1))
      end do
      write(luout,'(x,76("-"))')
      write(luout,'(x,a)') 'Number of terms with n-fold commutators'
      write(luout,'(x,a)') '   n       0       1       2       3'//
     &                 '       4       5       6       7       8'
      write(luout,'(x,76("-"))')
      write(luout,'(5x,9i8)') binsum(1:maxtop+1)
      write(luout,'(x,76("-"))')
      write(luout,'(x,a)') 'By number of T-T contractions'
      do ii = 0, maxtt
        if (maxval(bins(ii+1,1:maxtop+1)).eq.0) cycle
        write(luout,'(x,i4,9i8)') ii, bins(ii+1,1:maxtop+1)
        if (deldue2maxtt.and.maxcon_tt.ge.0.and.ii.ge.maxcon_tt)
     &      write(luout,'(x,a,i4,a)') 'Truncated at ',maxcon_tt,
     &                              ' T-T connections'
      end do
      write(luout,'(x,76("-"))')

      return
      end
      
      
