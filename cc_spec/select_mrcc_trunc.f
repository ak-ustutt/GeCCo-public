      subroutine select_mrcc_trunc(flist,labels,nlabels,mode,op_info)
*----------------------------------------------------------------------*
*     apply truncations to MRCC Lagrangian.
*     Each block of T is assigned a perturbation order by the user.
*     Five (predefined) parts of the Hamiltonian as well.
*     The five parts of H are:
*     1.) F0 (diagonal blocks of Fock op. + effective parts)
*     2.) purely active two-particle operator
*     3.) Difference between Fock operator and F0
*     4.) two-electron operators with up to two active lines
*     5.) two-electron operators with three active lines
*     matthias, oct. 2011
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
     &     delete, error, replace, count_l, prescreen, h1bar, 
     &     is_rept
      integer ::
     &     idxtop, idxham, norder, len, idxfeff, iorder, iblk, 
     &     ihampart, nact, iblknew, idx_op, nrank, cgastp, agastp,
     &     ivtx, ii, nvtx, nham, ham_vtx, iterm, idxl, t1ord,
     &     lrank, trank, ntop, simp, mode_rept
      integer ::
     &     idxop(nlabels), occ_ham(ngastp,2)

      integer, pointer ::
     &     torder(:), horder(:), lorder(:)

      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(operator), pointer ::
     &     op_ham, op_feff, op_top
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, idxlist, iblk_occ
      logical, external ::
     &     zero_ivec
 
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'select_mrcc_trunc')
        write(lulog,*) 'mode = ',trim(mode)
      endif

      count_l = mode(1:7).eq.'COUNT_L'
      prescreen = mode(9:17).eq.'PRESCREEN'

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = .not.count_l.and.nlabels.ne.3.or.count_l.and.nlabels.ne.4
      
      if (error) then
        write(lulog,*) 'Error for operator labels:'
        do ii = 1, nlabels
          if (idxop(ii).le.0) then
            write(lulog,'(a20," - ??")') trim(labels(ii))
          else
            write(lulog,'(a20," - OK")') trim(labels(ii))
          end if
        end do
        if (.not.count_l.and.nlabels.ne.3.or.count_l.and.nlabels.ne.4)
     &       call quit(1,'select_mrcc_trunc','wrong number of labels')
        call quit(1,'select_mrcc_trunc','Labels not on list!')
      end if

      idxham  = idxop(1)
      idxfeff  = idxop(2)
      idxtop = idxop(3)
      if (count_l) idxl = idxop(4)
      op_ham => op_info%op_arr(idxham)%op
      op_feff => op_info%op_arr(idxfeff)%op
      op_top => op_info%op_arr(idxtop)%op
      iterm = 0

      call get_argument_value('method.MRCC','trunc_order',
     &     ival=norder)
      if (is_argument_set('method.MRCC','trunc_top').gt.0) then
        call get_argument_dimension(len,'method.MRCC','trunc_top')
        allocate(torder(len),horder(5),lorder(len))
        call get_argument_value('method.MRCC','trunc_top',
     &                          iarr=torder)
        lorder(1:len) = torder(1:len)
      else
        len = op_top%n_occ_cls 
        allocate(torder(len),horder(5),lorder(len))
        torder(1) = -1
      end if
      call get_argument_value('method.MRCC','trunc_ham',
     &                        iarr=horder)
      if (horder(1).gt.horder(3)) call quit(1,'select_mrcc_trunc',
     &            'Order of F0 may not exceed order of Fdiff')
      call get_argument_value('method.MRCC','rept_ham',
     &                        ival=mode_rept)

      ! prescreen: assume minimum pert. order for all Hamiltonian blks
      if (prescreen) horder(1:5) = minval(horder(1:5))

      ! default: automatically assign perturbation orders to blocks of T
      if (torder(1).lt.0) then
        call get_argument_value('method.MRCC','T1ord',
     &       ival=t1ord)
        ! default order assigned to T1: 0 if sequential exponential,
        !                               1 if simultaneous exponential
        if (t1ord.lt.0) then
          call get_argument_value('method.MRCC','H1bar',
     &         lval=h1bar)
          t1ord = 1
          if (h1bar) t1ord = 0
        end if

        do iblk = 1, op_top%n_occ_cls
          nrank = op_top%ica_occ(1,iblk)
          if (nrank.eq.1) then
            torder(iblk) = t1ord
            lorder(iblk) = min(1,t1ord) ! Lambda_1 at most first order
          else
            torder(iblk) = nrank - 1
            lorder(iblk) = nrank - 1
          end if
        end do
      end if

      ! special simplifications requested?
      call get_argument_value('method.MRCC','simp',
     &     ival=simp)

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

          iterm = iterm + 1
c dbg
c          print *,'current term: ',iterm
c dbgend
          contr => form_pnt%contr
          nvtx = contr%nvtx
          vertex => contr%vertex

          ! find out:
          ! sum of perturbation orders of all operators
          ! operators except H and T are assigned order zero
          iorder  = 0
          nham = 0
          ihampart = 0
          replace = .false.
          ntop = 0
          lrank = 0
          trank = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              iblk = vertex(ivtx)%iblk_op
              if (iblk.gt.len) call quit(1,'select_mrcc_trunc',
     &           'Please define perturbation order for all blocks of T')
              iorder = iorder + torder(iblk)
              ntop = ntop + 1
              trank = op_top%ica_occ(1,iblk)
            else if (count_l.and.idx_op.eq.idxl) then
              iblk = vertex(ivtx)%iblk_op
              if (iblk.gt.len) call quit(1,'select_mrcc_trunc',
     &           'Please define perturbation order for all blocks of T')
              iorder = iorder + lorder(iblk)
              lrank = op_info%op_arr(idxl)%op%ica_occ(1,iblk)
            else if (idx_op.eq.idxham) then
              nham = nham+1
              ham_vtx = ivtx
              iblk = vertex(ivtx)%iblk_op
              nrank = op_ham%ica_occ(1,iblk)
              select case(nrank)
              case(0)
                ihampart = 0
                nham = nham - 1 !just ignore scalar part
              case(1)
                cgastp = idxlist(1,op_ham%ihpvca_occ(1:ngastp,1,iblk),
     &                           ngastp,1)
                agastp = idxlist(1,op_ham%ihpvca_occ(1:ngastp,2,iblk),
     &                           ngastp,1)
                if (cgastp.eq.agastp) then
                  ihampart = 1
                  replace = cgastp.ne.IVALE
                  if (replace) occ_ham(1:ngastp,1:2) =
     &                         op_ham%ihpvca_occ(1:ngastp,1:2,iblk)
                else
                  ihampart = 3
                end if
              case(2)
                nact = op_ham%ihpvca_occ(IVALE,1,iblk)
     &                +op_ham%ihpvca_occ(IVALE,2,iblk)
                is_rept = .false.
                if (mode_rept.eq.1) is_rept = 
     &                    zero_ivec(op_ham%ihpvca_occ(1:ngastp,1,iblk)
     &                             -op_ham%ihpvca_occ(1:ngastp,2,iblk),
     &                              ngastp)
                if (mode_rept.eq.2) is_rept =
     &                 (op_ham%ihpvca_occ(IHOLE,1,iblk).eq.2.and.
     &                  op_ham%ihpvca_occ(IHOLE,2,iblk).eq.2) .or.
     &                 (op_ham%ihpvca_occ(IPART,1,iblk).eq.2.and.
     &                  op_ham%ihpvca_occ(IPART,2,iblk).eq.2) .or.
     &                 (op_ham%ihpvca_occ(IVALE,1,iblk).eq.2.and.
     &                  op_ham%ihpvca_occ(IVALE,2,iblk).eq.2) 
                select case(nact)
                case(0,1,2)
                  ihampart = 4
                  if (mode_rept.gt.0.and.is_rept) ihampart=2
                case(3)
                  ihampart = 5
                case(4)
                  ihampart = 2
                case default
                  call quit(1,'select_mrcc_trunc','impossible!')
                end select
              case default
                call quit(1,'select_mrcc_trunc','nrank>2 for H?!')
              end select
              if (ihampart.gt.0) iorder = iorder + horder(ihampart)
            end if
          end do
          if (nham.gt.1) call warn('select_mrcc_trunc',
     &       'More than one Hamiltonian? Are we completely crazy now?!')

          ! delete if total order is too high
          delete = iorder.gt.norder
          ! do we have to replace H by effective Fock operator?
          if (.not.delete.and.replace) then
            replace = replace.and.iorder-horder(1)+horder(3).gt.norder
            replace = replace.and..not.mode_rept.eq.1 ! no replacement for REPT
          end if
          ! special simplifications requested (for (T))?
          if (.not.delete.and.simp.eq.2) then
            delete = lrank.eq.3.and.(ntop.ne.1.or.trank.eq.1)
     &               .or.lrank.eq.0.and.ntop.eq.1.and.trank.eq.3
          else if (.not.delete.and.simp.eq.1) then
            delete = lrank.eq.3.and.
     &               (ntop.eq.0.or.ntop.eq.1.and.trank.eq.1)
     &               .or.lrank.eq.0.and.ntop.eq.1.and.trank.eq.3
          end if

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Deleted formula item:'
              write(lulog,*) 'iterm, iorder:',iterm,iorder
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)

          else if (replace) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*) 'Replaced formula item:'
              write(lulog,*) 'iterm, iorder:',iterm,iorder
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif
            ! determine matching block of effective Fock op.
            iblknew = iblk_occ(occ_ham,.false.,op_feff,
     &                op_ham%blk_version(vertex(ham_vtx)%iblk_op))
            if (iblknew.le.0) call quit(1,'select_mrcc_trunc',
     &           'no matching block found in effective Fock operator')
            ! replace
            vertex(ham_vtx)%idx_op = idxfeff
            vertex(ham_vtx)%iblk_op = iblknew
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'select_mrcc_trunc','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      deallocate(torder,horder,lorder)

      return
      end
