*----------------------------------------------------------------------*
      subroutine set_mrcc_intermediates(form_out,
     &     title,label_int,label_op,nop,max_n,fac_in,typ_str,op_info)
*----------------------------------------------------------------------*
*     set the formal definition of intermediates in MRCC
*     label_int: label of intermediate
*     label_op(1:nop): labels of the operators which contribute 
*              (1 -- T, 2 -- H, 3 -- Heff)
*     typ_str:  describes the operator
*         Heff, Geff
*
*     matthias, march 2011
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 00

      type(formula), intent(inout), target ::
     &     form_out

      integer, intent(in) ::
     &     nop, max_n
      character*(*), intent(in) ::
     &     title,
     &     label_int,
     &     label_op(nop),
     &     typ_str
      type(operator_info), intent(inout) ::
     &     op_info
      real(8), intent(in) ::
     &     fac_in

      integer ::
     &     idx, idx_intm, idx_t, idx_h, idx_heff, idx_op4, idx_op5,
     &     nvtx, nn, ii, jj, iblk_int, iblk_h, mini, tmax,
     &     iterm, nterm, icall, icall0, nsumcalls,
     &     iu(2), ho(2), hu(2), tto(2), navoid
      character ::
     &     name*(form_maxlen_label*2)
      type(formula_item), target ::
     &     flist, flist_scr
      type(formula_item), pointer ::
     &     flist_pnt, start_pnt1, start_pnt2, start_pnt3
      type(operator), pointer ::
     &     op_int, op_t, op_h, op_heff
      real(8) ::
     &     fac
      logical ::
     &     next, set_zero

      integer, allocatable ::
     &     idx_op(:), iblk_min(:), iblk_max(:),
     &     idx_op_vtx(:), dist(:), imnmx(:,:), perm(:), avoid(:)

      integer, external ::
     &     idx_oplist2, ifac
      logical, external ::
     &     next_dist, next_perm

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0,
     &     cpu1, sys1, wall1, cpu10, sys10, wall10

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'set_mrcc_intermediates')
        write(lulog,*) 'setting: ',trim(label_int)
        write(lulog,*) 'type : ',trim(typ_str)
        write(lulog,*) 'input ops: '
        do idx = 1, nop
          write(lulog,*) '"',trim(label_op(idx)),'"'
        end do
      end if

      if (nop.lt.1)
     &     call quit(1,'set_mrcc_intermediates',
     &                 'nop < 1 ?? too few! phew!')

      ! get indices of operators on label_list
      idx_intm = idx_oplist2(label_int,op_info)
      if (idx_intm.lt.0)
     &     call quit(1,'set_mrcc_intermediates',
     &     'label not on list (1): '//label_int)
      op_int => op_info%op_arr(idx_intm)%op

      idx_t = idx_oplist2(label_op(1),op_info)
      if (idx_t.lt.0)
     &     call quit(1,'set_mrcc_intermediates',
     &     'label not on list (2): '//label_op(1))
      op_t => op_info%op_arr(idx_t)%op

      if (trim(typ_str).ne.'D_INT') then
        idx_h = idx_oplist2(label_op(2),op_info)
        if (idx_h.lt.0)
     &       call quit(1,'set_mrcc_intermediates',
     &       'label not on list (3): '//label_op(2))
        op_h => op_info%op_arr(idx_h)%op
      end if

      if (trim(typ_str).eq.'Geff'.or.trim(typ_str).eq.'D_INT2') then
        idx_heff = idx_oplist2(label_op(3),op_info)
        if (idx_heff.lt.0)
     &       call quit(1,'set_mrcc_intermediates',
     &       'label not on list (4): '//label_op(3))
        op_heff => op_info%op_arr(idx_heff)%op
      end if
      if (trim(typ_str).eq.'D_FORM') then
        idx_heff = idx_oplist2(label_op(3),op_info)
        if (idx_heff.lt.0.and.max_n.le.3)
     &       call quit(1,'set_mrcc_intermediates',
     &       'label not on list (4): '//label_op(3))
      end if

      if (trim(typ_str).eq.'D_FORM'.and.max_n.le.3) then
        idx_op4 = idx_oplist2(label_op(4),op_info)
        if (idx_op4.lt.0)
     &       call quit(1,'set_mrcc_intermediates',
     &       'label not on list (5): '//label_op(4))
      end if
      if (trim(typ_str).eq.'D_FORM'.and.max_n.le.2) then
        idx_op5 = idx_oplist2(label_op(5),op_info)
        if (idx_op5.lt.0)
     &       call quit(1,'set_mrcc_intermediates',
     &       'label not on list (6): '//label_op(5))
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'definition of intermediate:'
        call print_op_occ(lulog,op_int)
      end if

      tmax = op_t%n_occ_cls

      ! set up scratch formula
      call init_formula(flist_scr)
      flist_pnt => flist_scr
      call new_formula_item(flist_pnt,
     &     command_set_target_init,idx_intm)
      start_pnt1 => flist_pnt
      flist_pnt => flist_pnt%next

      ! First part of intermediate
      select case(trim(typ_str))
      case('Heff','Geff') ! (He^T)_c
        do iblk_int = 1, op_int%n_occ_cls
          iu(1) = op_int%ihpvca_occ(1,2,iblk_int)
          iu(2) = op_int%ihpvca_occ(2,1,iblk_int)
          start_pnt1 => flist_pnt
          do iblk_h = 1, op_h%n_occ_cls
            ho(1) = op_h%ihpvca_occ(1,2,iblk_h)
            ho(2) = op_h%ihpvca_occ(2,1,iblk_h)
            ! upper inactive lines of H must be connected to Int
            if (ho(1).gt.iu(1).or.ho(2).gt.iu(2)) cycle
            hu(1) = op_h%ihpvca_occ(1,1,iblk_h)
            hu(2) = op_h%ihpvca_occ(2,2,iblk_h)
            do nn = 0, max_n
              fac = 1d0/dble(ifac(nn))
              ! each T has at least one upper inactive line
              if (iu(1)+iu(2)-ho(1)-ho(2)+hu(1)+hu(2).lt.nn) cycle
              nvtx = nn+3

              allocate(idx_op(nvtx),iblk_min(nvtx),iblk_max(nvtx),
     &                 idx_op_vtx(nvtx),dist(nn),imnmx(2,nn),perm(nn))

              idx_op_vtx(1) = idx_intm
              idx_op_vtx(nvtx) = idx_intm
              do ii = 1, nvtx-1
                idx_op(ii) = ii
              end do
              idx_op(nvtx) = 1

              if (nn.gt.0) then
                dist = 1
                dist(1) = 0
                imnmx(1,1:nn) = 1
                imnmx(2,1:nn) = tmax
              end if

              set_zero = nn.eq.0
              do while(next_dist(dist,nn,imnmx,1).or.set_zero)
               set_zero = .false.
               tto = 0
               do ii = 1, nn
                 tto(1) = tto(1)
     &               + op_t%ihpvca_occ(1,2,dist(ii))
                 tto(2) = tto(2)
     &               + op_t%ihpvca_occ(2,1,dist(ii))
               end do
               if (iu(1)+hu(1).ne.ho(1)+tto(1).or.
     &             iu(2)+hu(2).ne.ho(2)+tto(2)) then
                 mini = 0
                 do ii = nn-1,1,-1
                   if (mini.eq.0.and.dist(ii).eq.tmax) mini=dist(ii+1)
                   if (mini.gt.0) imnmx(1,ii) = mini+1
                 end do
                 cycle
               end if
               start_pnt2 => flist_pnt
               nsumcalls = 0
               ! must be in increasing order for next_perm
               do ii = 1, nn
                 perm(ii) = dist(nn+1-ii)
               end do

               next = .true.
               do while(next)

                !all T's to the right of Hamiltonian: kk = 0
                idx_op_vtx(2:nvtx-1) = idx_t
                idx_op_vtx(2) = idx_h
                iblk_min = 1 
                iblk_max = -1 
                iblk_min(1) = iblk_int
                iblk_max(1) = iblk_int
                iblk_min(nvtx) = iblk_int
                iblk_max(nvtx) = iblk_int
                iblk_min(2) = iblk_h
                iblk_max(2) = iblk_h
                jj = 0
                start_pnt3 => flist_pnt
                do ii = 3, nvtx-1
                  jj = jj + 1
                  iblk_min(ii) = perm(jj)
                  iblk_max(ii) = perm(jj)
                end do

c dbg
c                write(lulog,'(a,3i6)') 'iblk_int,nh,n:',
c     &                                  iblk_int,iblk_h,nn
c                write(lulog,'(a,12i4)') 'iblk_min: ',iblk_min
c                write(lulog,'(a,12i4)') 'iblk_max: ',iblk_max
c dbgend
                call atim_csw(cpu10,sys10,wall10)
                icall = icall + 1

                call expand_op_product2(flist_pnt,idx_intm,
     &               fac,nvtx,nvtx-1,
     &               idx_op_vtx,
     &               idx_op,
     &               iblk_min,iblk_max,
     &               0,0,
     &               0,0,
     &               0,0,
     &               .true.,op_info)

                call atim_csw(cpu1,sys1,wall1)
                iterm = 0
                do while(flist_pnt%command.ne.command_end_of_formula)
                  iterm = iterm + 1
                  flist_pnt => flist_pnt%next
                end do
                nterm = nterm + iterm
                if (iterm.eq.0) icall0 = icall0 + 1
c dbg
c                print *,'# of generated terms: ',iterm
c                call prtim(lulog,'time in expand_op_product2 ',
c       &                   cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                nsumcalls = nsumcalls + 1
                call atim_csw(cpu10,sys10,wall10)
                call sum_terms(start_pnt3,op_info)
                call atim_csw(cpu1,sys1,wall1)
c dbg
c              call prtim(lulog,'time in sum_terms (2) ',
c       &                 cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                next = next_perm(perm,nn)
               end do
               mini = 0
               do ii = nn-1,1,-1
                 if (mini.eq.0.and.dist(ii).eq.imnmx(2,1))
     &              mini = dist(ii+1)
                 if (mini.gt.0) imnmx(1,ii) = mini+1
               end do

               if (nsumcalls.ne.1) then
                 call atim_csw(cpu10,sys10,wall10)
                 call sum_terms(start_pnt2,op_info)
                 call atim_csw(cpu1,sys1,wall1)
c dbg
c                 call prtim(lulog,'time in sum_terms ',
c       &                    cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
               end if
              end do
              deallocate(idx_op_vtx,idx_op,iblk_min,iblk_max,dist,imnmx,
     &                   perm)
            end do
          end do
          ! discard disconnected terms
          call select_mrcc_lag2(start_pnt1,
     &         (/label_op(2),label_op(1)/),2,'---',op_info)

          ! if required: add -(e^T Heff)
          if (trim(typ_str).eq.'Geff') then
          start_pnt1 => flist_pnt
          do iblk_h = 1, op_heff%n_occ_cls
            ho(1) = op_heff%ihpvca_occ(1,2,iblk_h)
            ho(2) = op_heff%ihpvca_occ(2,1,iblk_h)
            ! upper inactive lines of H must be connected to Int
            if (ho(1).gt.iu(1).or.ho(2).gt.iu(2)) cycle
            hu(1) = op_heff%ihpvca_occ(1,1,iblk_h)
            hu(2) = op_heff%ihpvca_occ(2,2,iblk_h)
            do nn = 1, max_n ! there should be no scalar term here
              fac = -1d0/dble(ifac(nn))
              ! each T has at least one upper inactive line
              if (iu(1)+iu(2)-ho(1)-ho(2)+hu(1)+hu(2).lt.nn) cycle
              nvtx = nn+3

              allocate(idx_op(nvtx),iblk_min(nvtx),iblk_max(nvtx),
     &                 idx_op_vtx(nvtx),dist(nn),imnmx(2,nn),perm(nn))

              idx_op_vtx(1) = idx_intm
              idx_op_vtx(nvtx) = idx_intm
              do ii = 1, nvtx-1
                idx_op(ii) = ii
              end do
              idx_op(nvtx) = 1
              if (nn.gt.0) then
                dist = 1
                dist(1) = 0
                imnmx(1,1:nn) = 1
                imnmx(2,1:nn) = tmax
              end if

              set_zero = nn.eq.0
              do while(next_dist(dist,nn,imnmx,1).or.set_zero)
               set_zero = .false.
               tto = 0
               do ii = 1, nn
                 tto(1) = tto(1)
     &               + op_t%ihpvca_occ(1,2,dist(ii))
                 tto(2) = tto(2)
     &               + op_t%ihpvca_occ(2,1,dist(ii))
               end do
               if (iu(1).ne.ho(1)+tto(1).or.
     &             iu(2).ne.ho(2)+tto(2)) then
                 mini = 0
                 do ii = nn-1,1,-1
                   if (mini.eq.0.and.dist(ii).eq.tmax) mini=dist(ii+1)
                   if (mini.gt.0) imnmx(1,ii) = mini+1
                 end do
                 cycle
               end if
               start_pnt2 => flist_pnt
               nsumcalls = 0
               ! must be in increasing order for next_perm
               do ii = 1, nn
                 perm(ii) = dist(nn+1-ii)
               end do

               next = .true.
               do while(next)

                ! all T's to the left of Hamiltonian: kk = nn
                idx_op_vtx(2:nvtx-1) = idx_t
                idx_op_vtx(2+nn) = idx_heff
                iblk_min = 1 
                iblk_max = -1 
                iblk_min(1) = iblk_int
                iblk_max(1) = iblk_int
                iblk_min(nvtx) = iblk_int
                iblk_max(nvtx) = iblk_int
                iblk_min(2+nn) = iblk_h
                iblk_max(2+nn) = iblk_h
                jj = 0
                start_pnt3 => flist_pnt
                do ii = 2, nn+1
                  jj = jj + 1
                  iblk_min(ii) = perm(jj)
                  iblk_max(ii) = perm(jj)
                end do

c dbg
c                write(lulog,'(a,3i6)') 'iblk_int,nh,n:',
c     &                                  iblk_int,iblk_h,nn
c                write(lulog,'(a,12i4)') 'iblk_min: ',iblk_min
c                write(lulog,'(a,12i4)') 'iblk_max: ',iblk_max
c dbgend
                call atim_csw(cpu10,sys10,wall10)
                icall = icall + 1

                call expand_op_product2(flist_pnt,idx_intm,
     &               fac,nvtx,nvtx-1,
     &               idx_op_vtx,
     &               idx_op,
     &               iblk_min,iblk_max,
     &               0,0,
     &               0,0,
     &               0,0,
     &               .true.,op_info)

                call atim_csw(cpu1,sys1,wall1)
                iterm = 0
                do while(flist_pnt%command.ne.command_end_of_formula)
                  iterm = iterm + 1
                  flist_pnt => flist_pnt%next
                end do
                nterm = nterm + iterm
                if (iterm.eq.0) icall0 = icall0 + 1
c dbg
c                print *,'# of generated terms: ',iterm
c                call prtim(lulog,'time in expand_op_product2 ',
c       &                   cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                nsumcalls = nsumcalls + 1
                call atim_csw(cpu10,sys10,wall10)
                call sum_terms(start_pnt3,op_info)
                call atim_csw(cpu1,sys1,wall1)
c dbg
c              call prtim(lulog,'time in sum_terms (2) ',
c       &                 cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                next = next_perm(perm,nn)
               end do
               mini = 0
               do ii = nn-1,1,-1
                 if (mini.eq.0.and.dist(ii).eq.imnmx(2,1))
     &              mini = dist(ii+1)
                 if (mini.gt.0) imnmx(1,ii) = mini+1
               end do

               if (nsumcalls.ne.1) then
                 call atim_csw(cpu10,sys10,wall10)
                 call sum_terms(start_pnt2,op_info)
                 call atim_csw(cpu1,sys1,wall1)
c dbg
c                 call prtim(lulog,'time in sum_terms ',
c       &                    cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
               end if
              end do

              deallocate(idx_op_vtx,idx_op,iblk_min,iblk_max,dist,imnmx,
     &                   perm)
            end do
          end do
          ! discard disconnected terms
          call select_mrcc_lag2(start_pnt1,
     &         (/label_op(3),label_op(1)/),2,'---',op_info)
          end if
        end do

      case('Tred')
        ! T(n)^red= [O(T^2)]^|| ...
        iterm = 0
        do nn = 2, max_n
          fac = 1d0/dble(ifac(nn))
          nvtx = nn+4
          navoid = nn+3

          allocate(idx_op(nvtx),idx_op_vtx(nvtx),avoid(2*navoid))

          idx_op_vtx(1:nvtx) = idx_intm
          idx_op_vtx(2:nvtx-1) = idx_h
          idx_op_vtx(3:nvtx-2) = idx_t
          do ii = 1, nvtx-2
            idx_op(ii) = ii
          end do
          idx_op(nvtx-1) = 2
          idx_op(nvtx) = 1
          avoid(1:6) = (/1,nvtx-1,2,nvtx-1,2,nvtx/)
          do idx = 1, nn
            avoid(2*(2+idx)+1) = idx+2
            avoid(2*(3+idx)) = nvtx
          end do

          call expand_op_product2(flist_pnt,idx_intm,
     &         fac,nvtx,nvtx-2,
     &         idx_op_vtx,
     &         idx_op,
     &         -1,-1,
     &         0,0,
     &         avoid,navoid,
     &         0,0,
     &         .true.,op_info)

          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
          end do
          deallocate(idx_op_vtx,idx_op,avoid)
        end do
        ! no external lines from T
        call select_line(start_pnt1,idx_intm,idx_t,1,3,'no_ext')

        ! ... -O(T^2)
        do nn = 2, max_n
          fac = -1d0/dble(ifac(nn))
          nvtx = nn+2
          allocate(idx_op(nvtx),idx_op_vtx(nvtx))
          idx_op_vtx(1:nvtx) = idx_intm
          idx_op_vtx(2:nvtx-1) = idx_t
          do ii = 1, nvtx-1
            idx_op(ii) = ii
          end do
          idx_op(nvtx) = 1

          call expand_op_product2(flist_pnt,idx_intm,
     &         fac,nvtx,nvtx-1,
     &         idx_op_vtx,
     &         idx_op,
     &         -1,-1,
     &         0,0,
     &         0,0,
     &         0,0,
     &         .true.,op_info)

          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          deallocate(idx_op_vtx,idx_op)
        end do
        if (iterm.gt.0) call sum_terms(start_pnt1,op_info)
        ! delete disconnected terms or modify prefactors if requested
        call select_mrcc_wf(start_pnt1,label_op(1),op_info)

      case('D_INT','D_INT2')
        ! gamma_2 (or gamma_4)
        call expand_op_product2(flist_pnt,idx_intm,
     &       1d0,6,2,
     &       (/idx_intm,idx_t,idx_intm,idx_intm,idx_t,idx_intm/),
     &       (/1,2,1,1,2,1/),
     &       -1,-1,
     &       0,0,
     &       (/2,5/),1,
     &       0,0,
     &       .true.,op_info)

        do while(flist_pnt%command.ne.command_end_of_formula)
          flist_pnt => flist_pnt%next
          iterm = iterm + 1
        end do

        ! + fac_in * (gamma_1)^2
        if (trim(typ_str).eq.'D_INT') then
          call expand_op_product2(flist_pnt,idx_intm,
     &         fac_in,8,3,
     &         (/idx_intm,idx_t,idx_t,idx_intm,idx_intm,
     &           idx_t,idx_t,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         -1,-1,
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
        else
          call expand_op_product2(flist_pnt,idx_intm,
     &         fac_in,8,3,
     &         (/idx_intm,idx_h,idx_heff,idx_intm,idx_intm,
     &           idx_heff,idx_h,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         -1,-1,
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
        end if

        do while(flist_pnt%command.ne.command_end_of_formula)
          flist_pnt => flist_pnt%next
          iterm = iterm + 1
        end do

      case('D_FORM') ! define cumulant approximation for RDMs

        select case(max_n) !maximum cumulant rank
        case(2)
          ! gamma_3: D_INT08 * gamma_1
          call expand_op_product2(flist_pnt,idx_intm,
     &         1d0,8,3,
     &         (/idx_intm,idx_t,idx_intm,idx_intm,idx_intm,
     &           idx_intm,idx_t,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/3,1,1/),(/3,1,1/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! gamma_4: 0.5* D_INT09 * D_INT10
          call expand_op_product2(flist_pnt,idx_intm,
     &         0.5d0,8,3,
     &         (/idx_intm,idx_h,idx_heff,idx_intm,idx_intm,
     &           idx_heff,idx_h,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/4,1,1/),(/4,1,1/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! gamma_5: 0.5 * D_INT11 * D_INT12 * gamma_1
          call expand_op_product2(flist_pnt,idx_intm,
     &         0.5d0,10,4,
     &         (/idx_intm,idx_op4,idx_op5,idx_intm,idx_intm,idx_intm,
     &           idx_intm,idx_op5,idx_op4,idx_intm/),
     &         (/1,2,3,4,1,1,4,3,2,1/),
     &         (/5,1,1,1/),(/5,1,1,1/),
     &         0,0,
     &         (/2,7,2,8,2,9,3,7,3,8,3,9,4,7,4,8,4,9/),9,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
        case(3)
          ! gamma_4: gamma_3 * gamma_1
          call expand_op_product2(flist_pnt,idx_intm,
     &         1d0,8,3,
     &         (/idx_intm,idx_intm,idx_intm,idx_intm,idx_intm,
     &           idx_intm,idx_intm,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/4,3,1/),(/4,3,1/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! + 0.5 * D_INT04 * D_INT05
          call expand_op_product2(flist_pnt,idx_intm,
     &         0.5d0,8,3,
     &         (/idx_intm,idx_t,idx_h,idx_intm,idx_intm,
     &           idx_h,idx_t,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/4,1,1/),(/4,1,1/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! gamma_5: gamma_3 * gamma_2
          call expand_op_product2(flist_pnt,idx_intm,
     &         1d0,8,3,
     &         (/idx_intm,idx_intm,idx_intm,idx_intm,idx_intm,
     &           idx_intm,idx_intm,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/5,3,2/),(/5,3,2/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! -0.5 * D_INT06 * D_INT07 * gamma_1
          call expand_op_product2(flist_pnt,idx_intm,
     &         -0.5d0,10,4,
     &         (/idx_intm,idx_heff,idx_op4,idx_intm,idx_intm,idx_intm,
     &           idx_intm,idx_op4,idx_heff,idx_intm/),
     &         (/1,2,3,4,1,1,4,3,2,1/),
     &         (/5,1,1,1/),(/5,1,1,1/),
     &         0,0,
     &         (/2,7,2,8,2,9,3,7,3,8,3,9,4,7,4,8,4,9/),9,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
        case(4)
          ! gamma_5: gamma_4 * gamma_1 or D_INT13 * gamma_1
          if (idx_heff.ge.0) then
            call expand_op_product2(flist_pnt,idx_intm,
     &           1d0,8,3,
     &           (/idx_intm,idx_intm,idx_intm,idx_intm,idx_intm,
     &             idx_intm,idx_intm,idx_intm/),
     &           (/1,2,3,1,1,3,2,1/),
     &           (/5,4,1/),(/5,4,1/),
     &           0,0,
     &           (/2,6,2,7,3,6,3,7/),4,
     &           0,0,
     &           .true.,op_info)
          else
            call expand_op_product2(flist_pnt,idx_intm,
     &           1d0,8,3,
     &           (/idx_intm,idx_h,idx_intm,idx_intm,idx_intm,
     &             idx_intm,idx_h,idx_intm/),
     &           (/1,2,3,1,1,3,2,1/),
     &           (/5,1,1/),(/5,1,1/),
     &           0,0,
     &           (/2,6,2,7,3,6,3,7/),4,
     &           0,0,
     &           .true.,op_info)
          end if
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          ! + gamma_3 * D_INT01
          call expand_op_product2(flist_pnt,idx_intm,
     &         1d0,8,3,
     &         (/idx_intm,idx_intm,idx_t,idx_intm,idx_intm,
     &           idx_t,idx_intm,idx_intm/),
     &         (/1,2,3,1,1,3,2,1/),
     &         (/5,3,1/),(/5,3,1/),
     &         0,0,
     &         (/2,6,2,7,3,6,3,7/),4,
     &         0,0,
     &         .true.,op_info)
          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
            iterm = iterm + 1
          end do
          if (idx_heff.ge.0) then
            ! - D_INT02 * D_INT03 * gamma_1
            call expand_op_product2(flist_pnt,idx_intm,
     &           -1d0,10,4,
     &           (/idx_intm,idx_h,idx_heff,idx_intm,idx_intm,idx_intm,
     &             idx_intm,idx_heff,idx_h,idx_intm/),
     &           (/1,2,3,4,1,1,4,3,2,1/),
     &           (/5,1,1,1/),(/5,1,1,1/),
     &           0,0,
     &           (/2,7,2,8,2,9,3,7,3,8,3,9,4,7,4,8,4,9/),9,
     &           0,0,
     &           .true.,op_info)
            do while(flist_pnt%command.ne.command_end_of_formula)
              flist_pnt => flist_pnt%next
              iterm = iterm + 1
            end do
          end if

        case default
          write(lulog,*) 'maxcum: ', max_n
          call quit(1,'set_mrcc_intermediates','Not yet available')
        end select

      case default
        call quit(1,'set_mrcc_intermediates',
     &       'unknown type: '//trim(typ_str))
      end select

      ! prepare file:
      form_out%comment = trim(title)
      write(name,'(a,".fml")') trim(form_out%label)
      call file_init(form_out%fhand,name,ftyp_sq_unf,0)
      ! and write list
      call write_form_list(form_out%fhand,flist_scr,form_out%comment)

      if (ntest.ge.100) then
        write(lulog,*) 'final formula: ',trim(op_int%name)
        call print_form_list(lulog,flist_scr,op_info)
      end if

      call dealloc_formula_list(flist_scr)

      return
      end
