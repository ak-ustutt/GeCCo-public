*----------------------------------------------------------------------*
      subroutine set_mrcc_lagrangian(form_lag,
     &     title,label,nlabels,max_n,max_n_en,approx,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*     lagrangian for internally contracted MRCC
*
*     matthias, summer 2010
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 10

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_lag
      integer, intent(in) ::
     &     max_n, max_n_en, nlabels
      character(*), intent(in) ::
     &     label(nlabels), title, approx

      type(operator_info), intent(inout) ::
     &     op_info

      type(orbinf), intent(inout) ::
     &     orb_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     flist_pnt, start_pnt, start_pnt2

      integer ::
     &     ilabel, idx, idxham,idxtbar, idxt, idxen, idxd,
     &     nvtx, min_n, nn, ii, jj, iblk, mini, maxc, off, nl, nham,
     &     kk, iblk_l, iblk_ham, iterm, icall, icall0, nterm,
     &     lu(2), ho(2), hu(2), tto(2), tto_l(2), idxddag, ioff, idx_h,
     &     max_n0, nsumcalls, G_level

      logical ::
     &     next, set_zero, set_scalar, esym, sym

      real(8) ::
     &     fac

      integer, allocatable ::
     &     idx_op(:), iblk_min(:), iblk_max(:), avoid(:,:),
     &     idx_op_vtx(:), dist(:), imnmx(:,:), perm(:)
 
      integer, external::
     &     idx_oplist2, ifac
      logical, external ::
     &     next_dist, next_perm

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0,
     &     cpu1, sys1, wall1, cpu10, sys10, wall10

      off = 0

      if (ntest.eq.100) then
        call write_title(luout,wst_dbg_subr,
     &       'output from set_mrcc_lagrangian')
        write(luout,*) ' max_n = ',max_n
        write(luout,*) ' approx = ',trim(approx)
      end if

      min_n = 0
      ioff = 0
      if (approx(1:4).eq.'NORM') then
        nham = 0
        min_n = 1
      else
        nham = 1
      end if
      if (approx(1:5).eq.'FIX_N') min_n = max_n
      if (approx(1:4).eq.'HBAR') then
        min_n = max_n
        ioff = 1
      end if
      if (approx(1:4).eq.'EMAX') then
        read(approx(5:5),'(i1)') max_n0
      else
        max_n0 = max_n_en
      end if
      set_scalar = .not.approx(1:6).eq.'NOSCAL'
      sym = approx(1:3).eq.'SYM'
      esym = approx(1:4).eq.'ESYM'.or.sym
      if (.not.set_scalar) min_n = 1

      call get_argument_value('method.MRCC','G_level',
     &     ival=G_level)
      if (G_level.lt.0) G_level = max(max_n,max_n0) ! no approximation

      call atim_csw(cpu0,sys0,wall0)
      nterm = 0
      icall = 0
      icall0 = 0

      ! transform labels into indices in this way
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_mrcc_lagrangian',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idxen = idx
        if (ilabel.eq.2) idxtbar = idx
        if (ilabel.eq.3) idxham = idx
        if (ilabel.eq.4) idxt = idx
        if (ilabel.eq.5) idxd = idx
      end do
      idxddag = idxd
      if (ioff.eq.0) idxddag = -idxd

      maxc = op_info%op_arr(idxt)%op%n_occ_cls-off

      ! initialize formula
      call init_formula(flist)
      flist_pnt => flist
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxen)
      flist_pnt => flist_pnt%next

      fac = 1d0

      if (nham.eq.1.and.set_scalar) then
      ! <0| C0^+ e^(-T) H e^T C0 |0>
      start_pnt => flist_pnt
      if (ioff.eq.0) then
        idx_h = idxham
      else
        idx_h = idxen
      end if
      do iblk_ham = 1, op_info%op_arr(idx_h)%op%n_occ_cls
      ho(1) = op_info%op_arr(idx_h)%op%ihpvca_occ(1,2,iblk_ham)
      ho(2) = op_info%op_arr(idx_h)%op%ihpvca_occ(2,1,iblk_ham)
      ! no upper inactive lines of H
      if (ho(1)+ho(2).gt.0.and.ioff.eq.0.and..not.esym) cycle
      hu(1) = op_info%op_arr(idx_h)%op%ihpvca_occ(1,1,iblk_ham)
      hu(2) = op_info%op_arr(idx_h)%op%ihpvca_occ(2,2,iblk_ham)
      do nn = min_n, max_n0
       ! each T has at least one upper inactive line
       if (hu(1)+hu(2).lt.nn.and.ioff.eq.0.and..not.esym) cycle
       if (ioff.eq.1.and.ho(1)+ho(2)-hu(1)-hu(2)+4.lt.nn
     &     .and..not.esym) cycle

       allocate(idx_op(nn+2+nham),iblk_min(nn+2+nham),
     &          iblk_max(nn+2+nham),idx_op_vtx(nn+2+nham))

       idx_op_vtx(1) = idxddag
       idx_op_vtx(nn+2+nham) = idxd
       do ii = 1, nn+2+nham
         idx_op(ii) = ii+1-ioff
       end do
       if (ioff.eq.1) idx_op(nn+2+nham) = 1
       do kk = 0, min(G_level,nn) !nn ! 0 !no inactive lines on the left
        idx_op_vtx(2:nn+1+nham) = idxt
        idx_op_vtx(2+kk:1+kk+nham) = idxham
        iblk_min = 1
        iblk_max = -1
        if (ioff.eq.0) then
          iblk_max(1) = op_info%op_arr(idxd)%op%n_occ_cls
          iblk_min(2+kk:1+kk+nham) = iblk_ham
          iblk_max(2+kk:1+kk+nham) = iblk_ham
        else
          iblk_min(1) = iblk_ham
          iblk_max(1) = iblk_ham
        end if
        fac = 1d0/dble(ifac(kk)*ifac(nn-kk))
        if (mod(kk,2).ne.0) fac = -fac
        if (esym.and.nn.gt.0) fac = 0.5d0*fac

c dbg
c        write(luout,'(a,4i6)') 'nl,nh,n,k:',0,iblk_ham,nn,kk
c        write(luout,'(a,10i4)') 'iblk_min: ',iblk_min
c        write(luout,'(a,10i4)') 'iblk_max: ',iblk_max
c dbgend

        call atim_csw(cpu10,sys10,wall10)
        icall = icall + 1

        call expand_op_product2(flist_pnt,idxen,
     &       fac,nn+2+nham,nn+2+nham-ioff,
     &       idx_op_vtx,
     &       idx_op,
     &       iblk_min,iblk_max,
     &       0,0,
     &       0,0,
     &       0,0,
     &       .true.,op_info)

         call atim_csw(cpu1,sys1,wall1)

        iterm = 0
        do while(flist_pnt%command.ne.command_end_of_formula)
          iterm = iterm + 1
          flist_pnt => flist_pnt%next
        end do
        nterm = nterm + iterm
        if (iterm.eq.0) icall0 = icall0 + 1
c dbg
c        print *,'# of generated terms: ',iterm
c        call prtim(luout,'time in expand_op_product2 ',
c     &             cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
        if (esym.and.nn.gt.0) then
          idx_op_vtx(2:nn+1+nham) = -idxt
          idx_op_vtx(2+kk:1+kk+nham) = idxham
          fac = 0.5d0/dble(ifac(kk)*ifac(nn-kk))
          if (mod(nn-kk,2).ne.0) fac = -fac
c dbg
c          write(luout,'(a,4i6)') 'nl,nh,n,k:',0,iblk_ham,nn,kk
c          write(luout,'(a,10i4)') 'iblk_min: ',iblk_min
c          write(luout,'(a,10i4)') 'iblk_max: ',iblk_max
c dbgend  
          call atim_csw(cpu10,sys10,wall10)
          icall = icall + 1

          call expand_op_product2(flist_pnt,idxen,
     &         fac,nn+2+nham,nn+2+nham-ioff,
     &         idx_op_vtx,
     &         idx_op,
     &         iblk_min,iblk_max,
     &         0,0,
     &         0,0,
     &         0,0,
     &         .true.,op_info)

           call atim_csw(cpu1,sys1,wall1)

          iterm = 0
          do while(flist_pnt%command.ne.command_end_of_formula)
            iterm = iterm + 1
            flist_pnt => flist_pnt%next
          end do
          nterm = nterm + iterm
          if (iterm.eq.0) icall0 = icall0 + 1
c dbg
c          print *,'# of generated terms: ',iterm
c          call prtim(luout,'time in expand_op_product2 ',
c     &               cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend  
        end if

       end do
       deallocate(idx_op_vtx,idx_op,iblk_min,iblk_max)
      end do
      end do
      ! sum terms
      call atim_csw(cpu10,sys10,wall10)
      call sum_terms(start_pnt,op_info)
      call atim_csw(cpu1,sys1,wall1)
c dbg
c      call prtim(luout,'time in sum_terms ',
c     &           cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend

      do while(flist_pnt%command.ne.command_end_of_formula)
        flist_pnt => flist_pnt%next
      end do

      end if

      if (ioff.eq.0.and.(.not.esym.or.sym)) then
      ! <0| C0^+ L e^(-T) H e^T C0 |0>
      do nl = 1, op_info%op_arr(idxtbar)%op%n_occ_cls
        lu(1) = op_info%op_arr(idxtbar)%op%ihpvca_occ(1,1,nl)
        lu(2) = op_info%op_arr(idxtbar)%op%ihpvca_occ(2,2,nl)
        do iblk_ham = 1, op_info%op_arr(idxham)%op%n_occ_cls
          ho(1) = op_info%op_arr(idxham)%op%ihpvca_occ(1,2,iblk_ham)
          ho(2) = op_info%op_arr(idxham)%op%ihpvca_occ(2,1,iblk_ham)
          ! upper inactive lines of H must be connected to L
          if (ho(1).gt.lu(1).or.ho(2).gt.lu(2).and..not.sym) cycle
          hu(1) = op_info%op_arr(idxham)%op%ihpvca_occ(1,1,iblk_ham)
          hu(2) = op_info%op_arr(idxham)%op%ihpvca_occ(2,2,iblk_ham)
          do nn = min_n, max_n
            ! each T has at least one upper inactive line
            if (lu(1)+lu(2)-ho(1)-ho(2)+hu(1)+hu(2).lt.nn
     &          .and..not.sym) cycle
c            if (lu-ho+hu.gt.4*nn) cycle
            nvtx = nn+3+nham

            allocate(idx_op(nvtx),iblk_min(nvtx),iblk_max(nvtx),
     &               idx_op_vtx(nvtx),dist(nn),imnmx(2,nn),perm(nn))

            idx_op_vtx(1) = idxddag
            idx_op_vtx(nvtx) = idxd
            idx_op_vtx(2) = idxtbar
            do ii = 1, nvtx
              idx_op(ii) = ii+1
            end do

            if (nn.gt.0) then
              dist = 1
              dist(1) = 0
              imnmx(1,1:nn) = 1
              imnmx(2,1:nn) = maxc
            end if

            set_zero = nn.eq.0
            do while(next_dist(dist,nn,imnmx,1).or.set_zero)
             set_zero = .false.
             tto = 0
             do ii = 1, nn
               tto(1) = tto(1)
     &             + op_info%op_arr(idxt)%op%ihpvca_occ(1,2,dist(ii))
               tto(2) = tto(2)
     &             + op_info%op_arr(idxt)%op%ihpvca_occ(2,1,dist(ii))
             end do
             if ((lu(1)+hu(1).ne.ho(1)+tto(1).or.
     &           lu(2)+hu(2).ne.ho(2)+tto(2)).and..not.sym) then
               mini = 0
               do ii = nn-1,1,-1
                 if (mini.eq.0.and.dist(ii).eq.maxc) mini=dist(ii+1)
                 if (mini.gt.0) imnmx(1,ii) = mini+1
               end do
               cycle
             end if
             start_pnt => flist_pnt
             nsumcalls = 0
             ! must be in increasing order for next_perm
             do ii = 1, nn
               perm(ii) = dist(nn+1-ii)
             end do

             next = .true.
             do while(next)

              do kk = 0, min(G_level,nn) !nn
                idx_op_vtx(2+nham:nvtx-1) = idxt
                if (nham.eq.1) idx_op_vtx(3+kk) = idxham
                iblk_min = 1 
                iblk_max = -1 
                iblk_min(2) = nl 
                iblk_max(2) = nl 
                iblk_max(1) = op_info%op_arr(idxd)%op%n_occ_cls
                jj = 0
                tto_l = 0
                start_pnt2 => flist_pnt
                do ii = 3, nvtx-1
                 if (ii.eq.3+kk) then
                  iblk_min(ii) = iblk_ham
                  iblk_max(ii) = iblk_ham
                 else
                  jj = jj + 1
                  iblk_min(ii) = perm(jj)
                  iblk_max(ii) = perm(jj)
                  if (ii.lt.3+kk) then
                   tto_l(1) = tto_l(1)
     &               + op_info%op_arr(idxt)%op%ihpvca_occ(1,2,perm(jj))
                   tto_l(2) = tto_l(2)
     &               + op_info%op_arr(idxt)%op%ihpvca_occ(2,1,perm(jj))
                  end if
                 end if
                end do
                ! inactive lines of T left from H must connect to L
                if ((lu(1).lt.tto_l(1).or.lu(2).lt.tto_l(2).or.
     &              hu(1).gt.tto(1)-tto_l(1).or.
     &              hu(2).gt.tto(2)-tto_l(2)).and..not.sym) cycle
                fac = 1d0/dble(ifac(kk)*ifac(nn-kk))
                if (mod(kk,2).ne.0) fac = -fac
                if (sym) fac = 0.5d0*fac

c dbg
c                write(luout,'(a,4i6)') 'nl,nh,n,k:',nl,iblk_ham,nn,kk
c                write(luout,'(a,12i4)') 'iblk_min: ',iblk_min
c                write(luout,'(a,12i4)') 'iblk_max: ',iblk_max
c dbgend
                call atim_csw(cpu10,sys10,wall10)
                icall = icall + 1

                call expand_op_product2(flist_pnt,idxen,
     &               fac,nvtx,nvtx,
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
c                call prtim(luout,'time in expand_op_product2 ',
c     &                     cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                nsumcalls = nsumcalls + 1
                call atim_csw(cpu10,sys10,wall10)
                call sum_terms(start_pnt2,op_info)
                call atim_csw(cpu1,sys1,wall1)
c dbg
c            call prtim(luout,'time in sum_terms (2) ',
c     &                 cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
                if (sym) then
                 idx_op_vtx(1+nham:nvtx-2) = -idxt
                 if (nham.eq.1) idx_op_vtx(2+kk) = idxham
                 idx_op_vtx(nvtx-1) = -idxtbar
                 iblk_min = 1 
                 iblk_max = -1 
                 iblk_min(nvtx-1) = nl 
                 iblk_max(nvtx-1) = nl 
                 iblk_max(1) = op_info%op_arr(idxd)%op%n_occ_cls
                 jj = 0
                 tto_l = 0
                 do ii = 2, nvtx-2
                  if (ii.eq.2+kk) then
                   iblk_min(ii) = iblk_ham
                   iblk_max(ii) = iblk_ham
                  else
                   jj = jj + 1
                   iblk_min(ii) = perm(jj)
                   iblk_max(ii) = perm(jj)
                   if (ii.lt.2+kk) then
                    tto_l(1) = tto_l(1)
     &                + op_info%op_arr(idxt)%op%ihpvca_occ(1,2,perm(jj))
                    tto_l(2) = tto_l(2)
     &                + op_info%op_arr(idxt)%op%ihpvca_occ(2,1,perm(jj))
                   end if
                  end if
                 end do
                 ! inactive lines of T left from H must connect to L
                 if ((lu(1).lt.tto_l(1).or.lu(2).lt.tto_l(2).or.
     &               hu(1).gt.tto(1)-tto_l(1).or.
     &               hu(2).gt.tto(2)-tto_l(2)).and..not.sym) cycle
                 fac = 0.5d0/dble(ifac(kk)*ifac(nn-kk))
                 if (mod(nn-kk,2).ne.0) fac = -fac

c dbg
c                 write(luout,'(a,4i6)') 'nl,nh,n,k:',nl,iblk_ham,nn,kk
c                 write(luout,'(a,12i4)') 'iblk_min: ',iblk_min
c                 write(luout,'(a,12i4)') 'iblk_max: ',iblk_max
c dbgend
                 call atim_csw(cpu10,sys10,wall10)
                 icall = icall + 1

                 call expand_op_product2(flist_pnt,idxen,
     &                fac,nvtx,nvtx,
     &                idx_op_vtx,
     &                idx_op,
     &                iblk_min,iblk_max,
     &                0,0,
     &                0,0,
     &                0,0,
     &                .true.,op_info)

                 call atim_csw(cpu1,sys1,wall1)
                 iterm = 0
                 do while(flist_pnt%command.ne.command_end_of_formula)
                   iterm = iterm + 1
                   flist_pnt => flist_pnt%next
                 end do
                 nterm = nterm + iterm
                 if (iterm.eq.0) icall0 = icall0 + 1
c dbg
c                 print *,'# of generated terms: ',iterm
c                 call prtim(luout,'time in expand_op_product2 ',
c     &                      cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend

                 idx_op_vtx(2) = idxtbar
                end if
              end do
              next = next_perm(perm,nn)
             end do
             mini = 0
             do ii = nn-1,1,-1
              if (mini.eq.0.and.dist(ii).eq.imnmx(2,1)) mini=dist(ii+1)
              if (mini.gt.0) imnmx(1,ii) = mini+1
             end do

             if (nsumcalls.ne.1) then
               call atim_csw(cpu10,sys10,wall10)
               call sum_terms(start_pnt,op_info)
               call atim_csw(cpu1,sys1,wall1)
c dbg
c               call prtim(luout,'time in sum_terms ',
c     &                    cpu1-cpu10,sys1-sys10,wall1-wall10)
c dbgend
             end if

            end do
            deallocate(idx_op_vtx,idx_op,iblk_min,iblk_max,dist,imnmx,
     &                 perm)
          end do
        end do
      end do
      end if

      ! delete terms with factor zero (disconnected terms)
      if (nterm.gt.0) call del_zero_terms(flist,op_info,1d-12)

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist,op_info)
      end if

      ! assign comment
      form_lag%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_lag%label)
      call file_init(form_lag%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_lag%fhand,flist,
     &     form_lag%comment)

      call dealloc_formula_list(flist)

      if (ntest.ge.10) then
        write(luout,*) 'total # of initially generated terms: ',nterm
        write(luout,'(i6,a,i6,a)') icall0,' out of ',icall,
     &              ' calls (of exp_op_pr2) lead to zero terms'
      end if

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for formula generation',
     &           cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
