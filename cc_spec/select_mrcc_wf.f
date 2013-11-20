      subroutine select_mrcc_wf(flist,op_label,op_info)
*----------------------------------------------------------------------*
*     given a (piece of a ) CC-type wave function, two things are done:
*     1.) If disconnected terms should be forbidden, delete them
*     2.) Modify prefactors if an ansatz differing from e^T is used
*
*     matthias, Aug. 2011
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
      character(len=*), intent(in) ::
     &     op_label

      logical ::
     &     delete, first, alt_ansatz
      integer ::
     &     idxtop, idxl
      integer ::
     &     idx_op, ivtx, nvtx, iarc, vtx1, vtx2,
     &     ntt, ntop, idx_op1, idx_op2,
     &     ntesting, itesting, nj, ioff,
     &     jvtx, kvtx, nhash, iterm, n_extra, idx,
     &     nfac_r, tred
      integer(8) ::
     &     hash_cur
      real(8) ::
     &     fac_tot, fac_cur, x_ansatz, fac_r, fac_alt

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
        call write_title(lulog,wst_dbg_subr,'select_mrcc_wf')
      endif

      ! get operator indices
      idxtop = idx_oplist2(trim(op_label),op_info)
      
      if (idxtop.le.0)
     &    call quit(1,'select_mrcc_wf','Label not on list!')

      call get_argument_value('method.MRCC','Tred_mode',
     &     ival=tred)
      call get_argument_value('method.MRCC','x_ansatz',
     &     xval=x_ansatz)
      alt_ansatz = x_ansatz.eq.0d0.or.x_ansatz.eq.1d0

      n_extra = 0
      iterm = 0
      allocate(extra_term(1),extra_hash(1),extra_fac(1))

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
          nj = njres_contr(contr)
          vertex => contr%vertex

          allocate(bchpart(nvtx),connected(nvtx),testing(nvtx))
          bchpart = .false.
          connected = .false.
          testing = 0
          ntesting = 0

          ! find out:
          ! - number of T operators
          ntop  = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              ntop = ntop+1
              bchpart(ivtx) = .true.
              if (ntop.eq.1) then
                ntesting = 1
                testing(1) = ivtx
                connected(ivtx) = .true.
              end if
            end if
          end do

          ! number of T-T contractions
          ntt = 0
          do iarc = 1, contr%narc
            vtx1 = contr%arc(iarc)%link(1)
            vtx2 = contr%arc(iarc)%link(2)
            idx_op1 = vertex(vtx1)%idx_op
            idx_op2 = vertex(vtx2)%idx_op
            if (idx_op1.eq.idxtop.and.idx_op2.eq.idxtop) ntt = ntt + 1
          end do            

          ! are all Ts and H connected?
          itesting = 0
          do while(itesting.lt.ntesting)
            itesting = itesting + 1
            ivtx = testing(itesting)
            do iarc = 1, contr%narc
              if (ivtx.eq.contr%arc(iarc)%link(1)) then
                vtx2 = contr%arc(iarc)%link(2)
              else if (ivtx.eq.contr%arc(iarc)%link(2)) then
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
            if (ntesting.eq.ntop) exit ! all vtxs connected
          end do

          ! delete discon. term (if requested)
          delete = ntesting.ne.ntop.and.tred.eq.2

          ! now we need "testing" to comprise all T operators
          ntesting = 0
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
            if (idx_op.eq.idxtop) then
              ntesting = ntesting + 1
              testing(ntesting) = ivtx
            end if
          end do

          ! factor checking
          ! CAUTION: This procedure depends on a unique hash function!
          ! CAUTION: Will only work if open lines are at
          !          extreme ends of the diagram
          if (.not.delete.and.ntesting.gt.0) then
            fac_tot = 0d0
            fac_alt = 0d0
            nhash = 0
            allocate(ivtx1(nvtx),topo1(nvtx,nvtx),
     &         xlines1(nvtx,nj),svtx1(nvtx),ireo(nvtx),
     &         iperm(ntesting),ivtx2(nvtx),topo2(nvtx,nvtx),
     &         xlines2(nvtx,nj),hash_list(1))
            call pack_contr(svtx1,ivtx1,topo1,xlines1,contr,nj)
c            if (x_ansatz.ne.0.5d0) call isort(testing,ntesting,1)
c dbg
c            write(lulog,*) ' initial topology:'
c            call prt_contr_p(lulog,svtx1,ivtx1,topo1,xlines1,nvtx,nj)
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
c                        write(lulog,'(x,a,14i4)') 'skip to perm:',iperm
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
c                write(lulog,'(x,a,14i4)') 'reordering array: ',
c     &                                    ireo(1:nvtx)
c                write(lulog,*) ' reordered topology:'
c                call prt_contr_p(lulog,svtx1,ivtx2,topo2,xlines2,nvtx,
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
                fac_cur = 1d0/dble(ifac(ntop))
                fac_tot = fac_tot + fac_cur
                ! factor to be changed? Determine correct factor
                if (alt_ansatz) then
                  allocate(topo_loc(ntop,ntop))
                  do jvtx = 1, ntop
                    do ivtx = 1, ntop
                      topo_loc(ivtx,jvtx) =
     &                      topo2(testing(ivtx),testing(jvtx)).ne.0
                    end do
                  end do
                  if (x_ansatz.eq.1d0) then ! {e^-T}^-1
                    nfac_r = n_possible_reos(topo_loc,ntop)
                    fac_r = 1d0/dble(nfac_r)
                  else !x_ansatz=0d0: {e^T}
                    if (any(topo_loc(1:ntop,1:ntop))) then
                      fac_r = 0d0
                    else
                      fac_r = 1d0/dble(ifac(ntop))
                    end if
                  end if
                  deallocate(topo_loc)
                  fac_alt = fac_alt + fac_r
c dbg
c                  print *,'right side factor:',fac_r
c dbgend
                end if
c dbg
c                print *,'current factor contribution: ',fac_cur
c dbgend
              end if
              if (.not.next_perm(iperm,ntesting)) exit perm_loop
            end do perm_loop
            ! all factors will at least have the same sign
            ! this allows us to adjust the factors to a global sign
            fac_tot = sign(fac_tot,contr%fac)
            fac_alt = sign(fac_alt,contr%fac)
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
                write(lulog,'(x,2(a,i12),a,e10.3)') 'added to term ',
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
                write(lulog,'(x,a,i12,a,e10.3)') 'term ',iterm,
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
                contr%fac = x_ansatz
            end if
            delete = delete.or.abs(contr%fac).lt.1d-12

            deallocate(ivtx1,topo1,xlines1,svtx1,ireo,iperm,
     &                 ivtx2,topo2,xlines2,hash_list)
          end if

          deallocate(bchpart,connected,testing)

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
          call quit(1,'select_mrcc_wf','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      do idx = 1, n_extra
        if (abs(extra_fac(idx)).gt.1d-12) then
          write(lulog,'(x,a,i12,a,e10.3)')
     &         'OH NO! Virtual term belonging to term ',
     &         extra_term(idx),' has a factor ',extra_fac(idx)
          call warn('select_mrcc_wf',
     &              'FATAL: We messed up the formula!')
        end if
      end do
      deallocate(extra_hash,extra_term,extra_fac)

      return
      end
