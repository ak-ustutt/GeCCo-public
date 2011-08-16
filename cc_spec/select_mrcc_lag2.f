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
     &     delete, error, deldue2maxtt, check_fac, first, check, proj,
     &     alt_ansatz
      integer ::
     &     idxtop, idxham, idxc, idxl, nprojdel
      integer ::
     &     ii, idx_op, ivtx, nvtx, iarc, vtx1, vtx2,
     &     ntt, ntop, nham, idx_op1, idx_op2, itop,
     &     ntesting, itesting, maxcon_tt, ntt_save, nj, ioff,
     &     jvtx, kvtx, nhash, ham_vtx, kk, iterm, n_extra, idx,
     &     nopen(2), nclos(2), vtxc1, vtxc2, nvtxc, cnt(ngastp,2),
     &     nvtxl, vtxl, nfac_l, nfac_r
      integer ::
     &     idxop(nlabels), bins(maxtt+1,maxtop+1), binsum(maxtop+1)
      integer(8) ::
     &     hash_cur
      real(8) ::
     &     fac_tot, fac_cur, x_ansatz, fac_l, fac_r, fac_alt

      integer, allocatable ::
     &     testing(:)
      logical, allocatable ::
     &     connected(:), bchpart(:), topo_loc(:,:)
      integer, pointer ::
     &     svtx1(:),svtx2(:),ireo(:),iperm(:), extra_term(:), itmp(:)
      integer(8), pointer ::
     &     ivtx1(:),topo1(:,:),xlines1(:,:),
     &     ivtx2(:),topo2(:,:),xlines2(:,:), hash_list(:), i8tmp(:),
     &     extra_hash(:)
      real(8), pointer ::
     &     extra_fac(:), r8tmp(:)


      type(contraction), pointer ::
     &     contr
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

      integer, external ::
     &     idx_oplist2, njres_contr, idxlist, i8mltlist, ifac,
     &     i8common_entry, n_possible_reos
      integer(8), external ::
     &     topo_hash
      logical, external ::
     &     next_perm

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'select_mrcc_lag2')
        write(luout,*) 'mode = ',trim(mode)
      endif
      check = mode(1:5).eq.'CHECK'

      ! get operator indices
      error = .false.
      do ii = 1, nlabels
        idxop(ii) = idx_oplist2(trim(labels(ii)),op_info)
        error = error.or.idxop(ii).le.0
      end do
      error = nlabels.ne.2.and.nlabels.ne.4
      proj = nlabels.eq.4
      
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
     &       call quit(1,'select_mrcc_lag2','need 2 or 4 labels')
        call quit(1,'select_mrcc_lag2','Labels not on list!')
      end if

      call get_argument_value('method.MRCC','maxtt',
     &     ival=maxcon_tt)
      call get_argument_value('method.MRCC','x_ansatz',
     &     xval=x_ansatz)
      alt_ansatz = x_ansatz.eq.0d0.or.x_ansatz.eq.1d0

      idxham  = idxop(1)
      idxtop  = idxop(2)
      if (proj) then
        idxc    = idxop(3)
        idxl    = idxop(4)
      end if

      bins = 0
      n_extra = 0
      iterm = 0
      nprojdel = 0
      allocate(extra_term(1),extra_hash(1),extra_fac(1))

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
          nj = njres_contr(contr)
          vertex => contr%vertex

          allocate(bchpart(nvtx),connected(nvtx),testing(nvtx))
          bchpart = .false.
          connected = .false.
          deldue2maxtt = .false.
          testing = 0
          ntesting = 0

          ! find out:
          ! - number of T operators
          ! - vertex number of Hamiltonian
          ntop  = 0
          nham  = 0
          nvtxc = 0
          nvtxl = 0
          ham_vtx = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              bchpart(ivtx) = .true.
              if (nham.eq.0) ham_vtx = ham_vtx + 1
            end if
            if (idx_op.eq.idxham) then
              nham = nham+1
              bchpart(ivtx) = .true.
              ntesting = ntesting + 1
              testing(ntesting) = ivtx
              connected(ivtx) = .true.
              ham_vtx = ham_vtx + 1
            end if
            if (proj) then
              if (idx_op.eq.idxc) then
                nvtxc = nvtxc+1
                if (nvtxc.eq.1) vtxc1 = ivtx
                if (nvtxc.eq.2) vtxc2 = ivtx
              else if (idx_op.eq.idxl) then
                nvtxl = nvtxl+1
                vtxl = ivtx
              end if
            end if
          end do
          if (nham.gt.1) call warn('select_mrcc_lag2',
     &       'More than one Hamiltonian? Are we completely crazy now?!')
          if (proj.and.(nvtxc.gt.2.or.nvtxl.gt.1))
     &          call quit(1,'set_mrcc_lag2','max. C0^+, C0, and L')

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

          ! delete discon. term (only if Hamiltonian was found at all)
          delete = ntesting.ne.ntop+nham.and.nham.gt.0
c          delete = ntesting.ne.ntop+nham 
          if (delete.and.check.and.abs(contr%fac).ge.1d-12) 
     &       write(luout,'(x,a,i12)')
     &       'Deleting nonzero disconnected term with number: ',iterm

          ! delete if more T-T connections than requested
          if (maxcon_tt.ge.0) then
            delete = delete.or.ntt_save.gt.maxcon_tt
            deldue2maxtt = deldue2maxtt.or.ntt_save.gt.maxcon_tt
c            iterm = iterm - 1 ! do not count as term
          end if

          ! delete if term will be projected out anyways
          if (.not.delete.and.proj.and.nvtxc.eq.2.and.nvtxl.eq.1) then
            nopen = 0
            nclos = 0
            ! get number of open and closed lines to C0^+ and C0
            do iarc = 1, contr%narc
              if (contr%arc(iarc)%link(1).eq.vtxc1) then
                if (contr%arc(iarc)%link(2).eq.vtxl) then
                  nopen(1) = nopen(1)
     &                     + sum(contr%arc(iarc)%occ_cnt(1:ngastp,1:2))
                else if (contr%arc(iarc)%link(2).ne.vtxc2) then
                  nclos(1) = nclos(1)
     &                     + sum(contr%arc(iarc)%occ_cnt(1:ngastp,1:2))
                end if
              else if (contr%arc(iarc)%link(2).eq.vtxc2) then
                if (contr%arc(iarc)%link(1).eq.vtxl) then
                  nopen(2) = nopen(2)
     &                     + sum(contr%arc(iarc)%occ_cnt(1:ngastp,1:2))
                else if (contr%arc(iarc)%link(1).ne.vtxc1) then
                  nclos(2) = nclos(2)
     &                     + sum(contr%arc(iarc)%occ_cnt(1:ngastp,1:2))
                end if
              end if
            end do
            ! delete if more open line pairs than closed line pairs
            delete = minval(nopen(1:2)).gt.minval(nclos(1:2))
            if (delete) nprojdel = nprojdel + 1
          end if

          ! factor checking
          ! CAUTION: This procedure depends on a unique hash function!
          ! CAUTION: Will only work if open lines are at
          !          extreme ends of the diagram
          check_fac = mode(1:9).eq.'CHECK_FAC'.and.
     &                .not.delete.and.ntesting.gt.0
          if (check_fac) then
            fac_tot = 0d0
            fac_alt = 0d0
            nhash = 0
            allocate(ivtx1(nvtx),topo1(nvtx,nvtx),
     &         xlines1(nvtx,nj),svtx1(nvtx),ireo(nvtx),
     &         iperm(ntesting),ivtx2(nvtx),topo2(nvtx,nvtx),
     &         xlines2(nvtx,nj),hash_list(1))
            call pack_contr(svtx1,ivtx1,topo1,xlines1,contr,nj)
            if (x_ansatz.ne.0.5d0) call isort(testing,ntesting,1)
c dbg
c            write(luout,*) ' initial topology:'
c            call prt_contr_p(luout,svtx1,ivtx1,topo1,xlines1,nvtx,nj)
c dbgend
            do ivtx = 1, nvtx
              ireo(ivtx) = ivtx
              if (ivtx.le.ntesting) iperm(ivtx) = ivtx
            end do
            first = .true.
            perm_loop: do
              ivtx2 = ivtx1
              topo2 = topo1
              xlines2 = xlines1
              if (.not.first) then
                ! get ireo
                ioff = 0
                do ivtx = 1, nvtx
                  if (bchpart(ivtx)) then
                    ireo(ivtx) = iperm(ivtx-ioff)+ioff
                  else
                    ioff = ioff+1
                    ireo(ivtx) = ivtx
                  end if
                end do
                ! is this reo allowed by the contractions?
                ! order of contracted vertices must not be changed
                idx = 0
                do ivtx = 1, nvtx
                  if (bchpart(ivtx)) idx = idx + 1
                  do jvtx = ivtx+1, nvtx
                    if (ireo(ivtx).gt.ireo(jvtx)) then
                      if (topo1(ireo(ivtx),ireo(jvtx)).ne.0) then
                        ! change ireo to replace op at pos. ivtx
                        ioff = ntesting
                        do kvtx = idx+1, ntesting
                          do while(idxlist(ioff,
     &                              iperm(1:kvtx-1),kvtx-1,1).gt.0)
                            ioff = ioff - 1
                          end do
                          iperm(kvtx) = ioff
                        end do
                        if (.not.next_perm(iperm,ntesting))
     &                     exit perm_loop
c dbg
c                        write(luout,'(x,a,14i4)') 'skip to perm:',iperm
c dbgend
                        cycle perm_loop
                      end if
                    end if
                  end do
                end do
                ! reorder vtx, topo, xlines
                call reoi8mat(ivtx2,ireo,nvtx,1,1)
                call reoi8mat(topo2,ireo,nvtx,nvtx,3) !both rows and cols
                call reoi8mat(xlines2,ireo,nvtx,nj,1)
c dbg
c                write(luout,'(x,a,14i4)') 'reordering array: ',
c     &                                    ireo(1:nvtx)
c                write(luout,*) ' reordered topology:'
c                call prt_contr_p(luout,svtx1,ivtx2,topo2,xlines2,nvtx,
c     &                           nj)
c dbgend
              end if
              first = .false.
              ! get hash value
              hash_cur = topo_hash(ivtx2,topo2,xlines2,nvtx,nj)
c dbg
c              print *,'current hash value: ',hash_cur
c dbgend
              ! go ahead if we didn't have this term before
              if (nhash.eq.0.or.
     &            i8mltlist(hash_cur,hash_list,nhash,1).eq.0) then
                ! store hash value
                nhash = nhash + 1
                allocate(i8tmp(nhash))
                if (nhash.gt.1) i8tmp(1:nhash-1) = hash_list(1:nhash-1)
                i8tmp(nhash) = hash_cur
                deallocate(hash_list)
                hash_list => i8tmp
                i8tmp => null()
                kk = idxlist(ham_vtx,iperm,ntesting,1)-1
                fac_cur = 1d0/dble(ifac(kk)*ifac(ntop-kk))
                if (mod(kk,2).ne.0) fac_cur = -fac_cur
                fac_tot = fac_tot + fac_cur
                ! factor to be changed? Determine correct factor
                if (alt_ansatz) then
                  ! T's to the left of H
                  allocate(topo_loc(kk,kk))
                  do jvtx = 1, kk
                    do ivtx = 1, kk
                      topo_loc(ivtx,jvtx) = 
     &                      topo2(testing(ivtx),testing(jvtx)).ne.0
                    end do
                  end do
                  if (x_ansatz.eq.0d0) then ! {e^T}^-1 on left side
                    nfac_l = n_possible_reos(topo_loc,kk)
                    fac_l = 1d0/dble(nfac_l)
                  else !x_ansatz=1d0: {e^-T} on left side
                    if (any(topo_loc(1:kk,1:kk))) then
                      fac_l = 0d0
                    else
                      fac_l = 1d0/dble(ifac(kk))
                    end if
                  end if
                  if (mod(kk,2).ne.0) fac_l = -fac_l
                  deallocate(topo_loc)
                  ! T's to the right of H
                  allocate(topo_loc(ntop-kk,ntop-kk))
                  do jvtx = 1, ntop-kk
                    do ivtx = 1, ntop-kk
                      topo_loc(ivtx,jvtx) =
     &                      topo2(testing(ivtx+kk+1),
     &                            testing(jvtx+kk+1)).ne.0
                    end do
                  end do
                  if (x_ansatz.eq.1d0) then ! {e^-T}^-1 on right side
                    nfac_r = n_possible_reos(topo_loc,ntop-kk)
                    fac_r = 1d0/dble(nfac_r)
                  else !x_ansatz=0d0: {e^T} on right side
                    if (any(topo_loc(1:ntop-kk,1:ntop-kk))) then
                      fac_r = 0d0
                    else
                      fac_r = 1d0/dble(ifac(ntop-kk))
                    end if
                  end if
                  deallocate(topo_loc)
                  fac_alt = fac_alt + fac_l*fac_r
c dbg
c                  print *,'left side factor: ',fac_l
c                  print *,'right side factor:',fac_r
c dbgend
                end if
c dbg
c                print *,'current factor contribution: ',fac_cur
c dbgend
              end if
              if (.not.next_perm(iperm,ntesting)) exit perm_loop
            end do perm_loop
c dbg
c            print *,'calculated factor: ',fac_tot
c            print *,'existing factor  : ',contr%fac
c            print *,'factor changed to: ',fac_alt
c dbgend

            ! if factor does not match:
            if (abs(fac_tot-contr%fac).gt.1d-12) then
              ! if we already had an identical term:
              ! add the factor and delete this term
              ! else: correct the factor and create new "extra" term
              idx = i8common_entry(extra_hash,hash_list,n_extra,nhash)
              if (idx.gt.0) then
                extra_fac(idx) = extra_fac(idx) + contr%fac
                delete = .true.
                write(luout,'(x,2(a,i12),a,e10.3)') 'added to term ',
     &                 extra_term(idx),
     &                 ': term ',iterm,' with factor ',contr%fac
              else
                n_extra = n_extra + 1
                allocate(i8tmp(n_extra),itmp(n_extra),r8tmp(n_extra))
                if (n_extra.gt.1) then
                  i8tmp(1:n_extra-1) = extra_hash(1:n_extra-1)
                  itmp(1:n_extra-1) = extra_term(1:n_extra-1)
                  r8tmp(1:n_extra-1) = extra_fac(1:n_extra-1)
                end if
                i8tmp(n_extra) = hash_list(1) !arbitrary: first hash value
                itmp(n_extra) = iterm
                r8tmp(n_extra) = contr%fac-fac_tot ! factor difference
                deallocate(extra_hash,extra_term,extra_fac)
                extra_hash => i8tmp
                extra_term => itmp
                extra_fac  => r8tmp
                i8tmp => null()
                itmp  => null()
                r8tmp => null()
                contr%fac = fac_tot
                write(luout,'(x,a,i12,a,e10.3)') 'term ',iterm,
     &              ': new virtual term with factor ',extra_fac(n_extra)
              end if
            end if

            ! modify factor?
            if (.not.delete.and.alt_ansatz) then
              contr%fac = fac_alt
c              delete = (abs(contr%fac).lt.1d-12)
            else if (.not.delete.and.x_ansatz.ne.0.5d0.and.ntop.eq.2
     &          .and.ntt.eq.1.
     &          .and.abs(abs(contr%fac)-0.5d0).lt.1d-12) then
              if (contr%fac.gt.0d0.and.ham_vtx.eq.1) then
                contr%fac = x_ansatz
              else if (contr%fac.gt.0d0) then
                contr%fac = 1d0-x_ansatz
              else if (ham_vtx.eq.1) then
                contr%fac = x_ansatz-1d0
              else if (ham_vtx.eq.3) then
                contr%fac = -x_ansatz
              else if (topo1(testing(2),testing(3)).ne.0) then
                contr%fac = x_ansatz-1d0
              else if (topo1(testing(2),testing(1)).ne.0) then
                contr%fac = -x_ansatz
              else
                call quit(1,'select_mrcc_lag2','should not happen')
              end if
c              delete = (abs(contr%fac).lt.1d-12)
            end if
            delete = delete.or.abs(contr%fac).lt.1d-12

            deallocate(ivtx1,topo1,xlines1,svtx1,ireo,iperm,
     &                 ivtx2,topo2,xlines2,hash_list)
          end if

          deallocate(bchpart,connected,testing)

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

      if (proj) write(luout,'(x,a,i12)')
     &   'Number of terms deleted because they will be projected out:',
     &       nprojdel

      if (check) then
        do itop = 0,maxtop
          binsum(itop+1) = sum(bins(1:maxtt+1,itop+1))
        end do
        write(luout,'(x,76("-"))')
        write(luout,'(x,a)') 'Number of terms with n-fold commutators'
        write(luout,'(x,a)') '   n       0       1       2       3'//
     &                   '       4       5       6       7       8'
        write(luout,'(x,76("-"))')
        write(luout,'(5x,9i8)') binsum(1:maxtop+1)
        write(luout,'(x,76("-"))')
        write(luout,'(x,a)') 'By number of T-T contractions'
        do ii = 0, maxtt
          if (maxval(bins(ii+1,1:maxtop+1)).eq.0) cycle
          write(luout,'(x,i4,9i8)') ii, bins(ii+1,1:maxtop+1)
          if (deldue2maxtt.and.maxcon_tt.ge.0.and.ii.ge.maxcon_tt)
     &        write(luout,'(x,a,i4,a)') 'Truncated at ',maxcon_tt,
     &                                ' T-T connections'
        end do
        write(luout,'(x,76("-"))')
      end if

      do idx = 1, n_extra
        if (abs(extra_fac(idx)).gt.1d-12) then
          write(luout,'(x,a,i12,a,e10.3)')
     &         'OH NO! Virtual term belonging to term ',
     &         extra_term(idx),' has a factor ',extra_fac(idx)
          call warn('select_mrcc_lag2',
     &              'FATAL: We messed up the formula!')
        end if
      end do
      deallocate(extra_hash,extra_term,extra_fac)

      return
      end
      
*----------------------------------------------------------------------*
      integer function n_possible_reos(conmat,ndim)
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim
      logical, intent(in) ::
     &     conmat(ndim,ndim)

      integer ::
     &     ii, jj, kk, off, perm(ndim)
      logical ::
     &     first
      integer, external ::
     &     idxlist
      logical, external ::
     &     next_perm

      n_possible_reos = 1
      if (ndim.eq.1.and.conmat(1,1)) call quit(1,'n_possible_reos',
     &       'a single element can never be connected!')
      if (ndim.le.1) return

c dbg
c      write(luout,*) '-------------------------'
c      write(luout,*) 'n_possible_reos in action'
c      write(luout,*) 'input connectivity matrix:'
c      do ii = 1, ndim
c        write(luout,*) conmat(ii,1:ndim)
c      end do
c dbgend

      do ii = 1, ndim
        perm(ii) = ii
      end do
      first = .true.
      n_possible_reos = 0
      perm_loop: do
        if (.not.first) then
          ! is this permutation allowed by the connectivity?
          ! order of contracted vertices must not be changed
          do ii = 1, ndim
            do jj = ii+1, ndim
              if (perm(ii).gt.perm(jj)) then
                if (conmat(perm(ii),perm(jj))) then
                  ! change perm to replace op at pos. ii
                  off = ndim
                  do kk = ii+1, ndim
                    do while(idxlist(off,perm(1:kk-1),kk-1,1).gt.0)
                      off = off - 1
                    end do
                    perm(kk) = off
                  end do
                  if (.not.next_perm(perm,ndim))
     &               exit perm_loop
c dbg
c                  write(luout,'(x,a,14i4)') 'skip to perm:',perm
c dbge
                  cycle perm_loop
                end if
              end if
            end do
          end do
        end if
        first = .false.
        n_possible_reos = n_possible_reos + 1
        if (.not.next_perm(perm,ndim)) exit perm_loop
      end do perm_loop
c dbg
c      write(luout,*) 'number of possible reorderings:',n_possible_reos
c dbgend

      return
      end
