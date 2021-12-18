      subroutine select_mrcc_f12(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     special selection routine for MRCC F12 methods
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
     &     delete, error, replace, red
      integer ::
     &     idxfeff, idxr12, idxham, idxtop, idxtbar, idxrsi, idxsr12
      integer ::
     &     trunc_e, trunc_t1, trunc_t2, ii, trunc_r12,
     &     idx_op, iblk_op, ivtx, nvtx, max2t1, max2r12,
     &     ntop, nt1, typ_tbar, ord_ham, nr12, nrdag,
     &     iblk_t1, iblknew, ham_vtx, nrsi, trunc_rsi,
     &     cgastp, agastp, nfeff, nrsidag
      integer ::
     &     idxop(nlabels), occ(ngastp,2), occ_ham(ngastp,2)


      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, iblk_occ, idxlist

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_mrcc_f12')
      endif

      ! get operator indices
      error = .false.
      select case(trim(mode))
      case('MRCC','mrcc')
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = error.or.nlabels.ne.5
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.5)
     &       call quit(1,'select_mrcc_f12','need exactly 5 labels')
        call quit(1,'select_mrcc_f12','Labels not on list!')
      end if

        trunc_e = 2
        trunc_r12 = 1
        trunc_t2 = 1
        max2r12 = 0

      idxtbar = idxop(1)
      idxham  = idxop(2)
      idxtop  = idxop(3)
      idxr12  = idxop(4)
      idxfeff = idxop(5)

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

          !nvtx: number of the vertex
          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! find out:
          ! - number of T operators
          ! - number of R+
          ! - number of R
          ! - number of T1
          ! - order and type of H
          ! - type of Tbar
          ntop  = 0
          typ_tbar = 0
          nrdag = 0
          nr12  = 0
          nt1 = 0
          nfeff = 0
          replace = .false.
          iblknew =0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              if (iblk_op.eq.iblk_t1 ) nt1  = nt1 + 1
              ntop = ntop + 1
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
              ham_vtx=ivtx
              ord_ham = 
     &             op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
              if (ord_ham.eq.0) then
                  iblknew = iblk_occ(op_info%op_arr(idxham)%
     &                      op%ihpvca_occ(1:ngastp,1:2,iblk_op),
     &                      .false.,
     &                      op_info%op_arr(idxfeff)%op,
     &                      op_info%op_arr(idxham)%op
     &                      %blk_version(vertex(ham_vtx)%iblk_op))
                if (iblknew.le.0) call quit(1,'select_mrcc_f12',
     &            'no matching block found in effective Fock operator')
                replace=.true.
              end if
            end if
          end do

          ! rules for deleting vertex
          delete=
     &    (nr12.gt.1.or.nrdag.ge.1.and.nr12.ge.1.and.ord_ham.ge.1)
          delete=delete.or.
     &    (nr12.ge.1.and.nrdag.ge.1.and.ntop.ge.1.and.ord_ham.eq.0)
          replace = replace .and. nr12.eq.1.and.nrdag.eq.1.and.ntop.eq.0

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
              write(lulog,*) 'nrdag,nr12,ntop,typ_tbar,ord_ham: ',
     &             nrdag,nr12,ntop,typ_tbar,ord_ham
            end if

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else if (replace) then
            ! Replace
            vertex(ham_vtx)%idx_op = idxfeff
            vertex(ham_vtx)%iblk_op = iblknew
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_mrcc_f12','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      case('MRCC_SI','mrcc_si','MRCC_SI_RED','mrcc_si_red')
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = error.or.nlabels.ne.7!6
      red = (mode.eq.'MRCC_SI_RED'.or.mode.eq.'mrcc_si_red')      

      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (nlabels.ne.6)
     &       call quit(1,'select_mrcc_f12','need exactly 6 labels')
        call quit(1,'select_mrcc_f12','Labels not on list!')
      end if

        trunc_e = 2
        trunc_r12 = 1
        trunc_rsi = 1
        trunc_t2 = 1
        max2r12 = 0

      idxtbar = idxop(1)
      idxham  = idxop(2)
      idxtop  = idxop(3)
      idxr12  = idxop(4)
      idxfeff = idxop(5)
      idxrsi  = idxop(6)
      idxsr12 = idxop(7)

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

          !nvtx: number of the vertex
          nvtx = form_pnt%contr%nvtx
          vertex => form_pnt%contr%vertex

          ! find out:
          ! - number of T operators
          ! - number of R+
          ! - number of R
          ! - number of Rsi
          ! - number of T1
          ! - order and type of H
          ! - type of Tbar
          ntop  = 0
          typ_tbar = 0
          nrdag = 0
          nr12  = 0
          nrsi  = 0
          nrsidag = 0
          nt1 = 0
          nfeff = 0
          replace = .false.
          iblknew =0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtop) then
              if (iblk_op.eq.iblk_t1 ) nt1  = nt1 + 1
              ntop = ntop + 1
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
            if (idx_op.eq.idxrsi .or. idx_op.eq.idxsr12) then
              if (vertex(ivtx)%dagger) then
                nrsidag = nrsidag+1
              else
                nrsi = nrsi+1
              end if
            end if

            if (idx_op.eq.idxham) then
              ham_vtx=ivtx
              ord_ham = 
     &             op_info%op_arr(idx_op)%op%ica_occ(1,iblk_op)-1
              if (ord_ham.eq.0) then
                  iblknew = iblk_occ(op_info%op_arr(idxham)%
     &                      op%ihpvca_occ(1:ngastp,1:2,iblk_op),
     &                      .false.,
     &                      op_info%op_arr(idxfeff)%op,
     &                      op_info%op_arr(idxham)%op
     &                      %blk_version(vertex(ham_vtx)%iblk_op))
                if (iblknew.le.0) call quit(1,'select_mrcc_f12',
     &            'no matching block found in effective Fock operator')
                replace=.true.
              end if
            end if
          end do

          ! rules for deleting vertex
          delete=
     &    (nr12+nrsi.gt.1.or.nrdag+nrsidag.ge.1
     &     .and.nr12+nrsi.ge.1.and.ord_ham.ge.1)
          delete=delete.or.
     &    (nr12+nrsi.ge.1.and.nrdag+nrsidag.ge.1
     &     .and.ntop.ge.1.and.ord_ham.eq.0)

          ! more strict delete rule for semi-internals:
          if (red) then
            delete=delete.or.
     &      (nrsidag.ge.1.and.ntop.gt.1).or.(nrsi.ge.1.and.ntop.gt.0)
          end if

          replace = replace .and. 
     &              nr12+nrsi.eq.1.and.nrdag+nrsidag.eq.1.and.ntop.eq.0

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000) then
              write(lulog,*) 'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
              write(lulog,*) 'nrdag,nr12,ntop,typ_tbar,ord_ham,
     &             nrsi,nrsidag: ',
     &             nrdag,nr12,ntop,typ_tbar,ord_ham,
     &             nrsi,nrsidag
            end if

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          else if (replace) then
            ! Replace
            vertex(ham_vtx)%idx_op = idxfeff
            vertex(ham_vtx)%iblk_op = iblknew
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_mrcc_f12','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo
      case default
        call quit(1,'select_mrcc_f12','mode?? "'//trim(mode)//'"')
      end select

      return
      end
      
      
