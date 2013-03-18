*----------------------------------------------------------------------*
      subroutine leqc_core(iter,
     &       task,iroute,xrsnrm,
     &       use_s,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,me_met,me_scr,
     &       me_special,nspecial,
c     &       ffopt,fftrv,ffmvp,ffrhs,ffdia,
     &       nincore,lenbuf,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,
     &       fspc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     core driver for LEQ solver
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
c      include 'def_filinf.h'
      include 'mdef_operator_info.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_dependency_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(inout) ::
     &     task
      integer, intent(inout) ::
     &     iter
      integer, intent(in) ::
     &     iroute, nincore, lenbuf, nspecial, nspcfrm
      logical, intent(in) ::
     &     use_s(*)

      type(me_list_array), intent(inout) ::
     &     me_opt(*), me_trv(*), me_dia(*),
     &     me_mvp(*), me_rhs(*), me_scr(*),
     &     me_special(nspecial)
      type(me_list_array), intent(inout) ::
     &     me_met(*)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffrhs(*), ffdia(*)

      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend
      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xrsnrm(opti_info%nroot,opti_info%nopt)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

* local
      logical ::
     &     zero_vec(opti_stat%ndim_vsbsp), init, conv, trafo, getnewrec
      integer ::
     &     idx, jdx, kdx, iroot, irhs,  nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job,
     &     ndim_save, ndel, iopt, lenscr, ifree, restart_mode,
     &     nselect, irec, ioff_s, ioff, nsec, jopt, isec, stsec, ndsec
      real(8) ::
     &     cond, xdum, xnrm
      real(8), pointer ::
     &     gred(:), vred(:), mred(:), sred(:),
     &     xmat1(:), xmat2(:), xvec(:), xret(:), signsec(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &     iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:),
     &     nwfpar(:), nwfpsec(:), idstsec(:), nsec_arr(:),
     &     ipiv(:), iconv(:), idxroot(:), idxselect(:)
      type(file_array), pointer ::
     &     ffrsbsp(:), ffvsbsp(:), ffssbsp(:), ffscr(:), ffmet(:)
      type(filinf), pointer ::
     &     ffspc
      type(filinf), target ::
     &     fdum

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'leqc_core entered')

      nopt = opti_info%nopt
      nroot = opti_info%nroot
      mxsub = opti_stat%mxdim_sbsp
      mred => opti_stat%sbspmat(1:)
      gred => opti_stat%sbspmat(mxsub**2+1:)
      vred => opti_stat%sbspmat(2*mxsub**2+1:)
      sred => opti_stat%sbspmat(3*mxsub**2+1:)
      ndim_rsbsp => opti_stat%ndim_rsbsp
      ndim_vsbsp => opti_stat%ndim_vsbsp
      ndim_ssbsp => opti_stat%ndim_ssbsp
      iord_rsbsp => opti_stat%iord_rsbsp
      iord_vsbsp => opti_stat%iord_vsbsp
      iord_ssbsp => opti_stat%iord_ssbsp
      ffscr => opti_stat%ffscr
      ! take care: there is currently some mess with 
      ! the definition of the below variables 
      ! (types filinf and file_array)
      ffrsbsp => opti_stat%ffrsbsp
      ffvsbsp => opti_stat%ffvsbsp
      ffssbsp => opti_stat%ffssbsp
      nwfpar => opti_info%nwfpar

      allocate(ffmet(nopt))

c      if (nopt.gt.1)
c     &     call quit(1,'leqc_core','not yet adapted for nopt>1')

      ! check for previously converged roots
      ifree = mem_alloc_int(iconv,nopt*nroot,'LEQ_conv')
      idx = 0
      do iopt = 1, nopt
        do iroot = 1, nroot
          idx = idx+1
          iconv(idx) = 0
          if (iter.gt.1.and.xrsnrm(iroot,iopt)
     &              .lt.opti_info%thrgrd(iopt))
     &         iconv(idx) = 1
        end do
      end do

      ifree = mem_alloc_int(idxroot,nopt*nroot,'LEQ_idxroot')

c      iopt = 1  ! preliminary
      if (ndim_vsbsp.ne.ndim_rsbsp)
     &     call quit(1,'leqc_core','subspace dimensions differ?')
      nred = ndim_vsbsp
      do iopt = 1, nopt
        init = iopt.eq.1
        nsec = opti_info%nsec(iopt)
        ioff = sum(opti_info%nsec(1:iopt))-nsec
        nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
        idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
        signsec => opti_info%signsec(ioff+1:ioff+nsec)
        ! update reduced space:
        ! ffvsbsp and ffrsbsp point to ff_trv(iopt)%fhand ...
        if (.not.use_s(iopt).and.nopt.eq.1) then
          call optc_update_redsp3
     &       (mred,gred,nred,nroot,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
     &       iord_rsbsp,ffrsbsp(iopt)%fhand,me_rhs(iopt)%mel%fhand,
     &       nsec,nwfpsec,idstsec,signsec,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        else if (.not.use_s(iopt)) then
          call optc_update_redsp4
     &         (mred,sred,gred,nred,nroot,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,init,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
     &       iord_rsbsp,ffrsbsp(iopt)%fhand,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,me_rhs(iopt)%mel%fhand,
     &       nsec,nwfpsec,idstsec,signsec,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        else
          call optc_update_redsp4
     &         (mred,sred,gred,nred,nroot,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,init,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
     &       iord_rsbsp,ffrsbsp(iopt)%fhand,
     &       iord_ssbsp,ffssbsp(iopt)%fhand,me_rhs(iopt)%mel%fhand,
c     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
c     &       iord_rsbsp,ffrsbsp(iopt)%fhand,
c     &       iord_ssbsp,ffssbsp(iopt)%fhand,fdum,
     &       nsec,nwfpsec,idstsec,signsec,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        end if
      end do

      zero_vec(1:opti_stat%ndim_vsbsp) = .false.

      ! ------------------------
      !    solve reduced LEQ
      ! ------------------------ 
      
      ! allocate some scratch 
      ! (automatically deallocated after leaving leq_control() )
      lenmat = nred*nred
      ifree = mem_alloc_real(xmat1,lenmat,'LEQ_mat')
      ifree = mem_alloc_real(xvec,nred,'LEQ_vec')
      ifree = mem_alloc_int (ipiv,nred,'LEQ_piv')

      ! get a copy of the subspace matrix
      ! and apply shift (incl. metric if applicable)
      if (opti_info%shift.eq.0d0.or..not.any(use_s(1:nopt))) then
        kdx = 0
        do idx = 1, nred
          do jdx = 1, nred
            kdx = kdx+1
            xmat1(kdx) = mred((idx-1)*mxsub+jdx)
            ! shift matrix
            if (idx.eq.jdx) xmat1(kdx) = xmat1(kdx) + opti_info%shift
          end do
        end do
      else
        kdx = 0
        do idx = 1, nred
          do jdx = 1, nred
            kdx = kdx+1
            xmat1(kdx) = mred((idx-1)*mxsub+jdx)
     &           + opti_info%shift * sred((idx-1)*mxsub+jdx)
          end do
        end do
      end if

      ! condition number and pivot vector
      call dgeco(xmat1,nred,nred,ipiv,cond,xvec)

      if (ntest.ge.10) write(luout,*)'dimension ,condition number: ',
     &         nred,cond

      if (ntest.ge.50) then
        write(luout,*) 'Factorized matrix from dgeco:'
        call wrtmat2(xmat1,nred,nred,nred,nred)
        write(luout,*) 'Pivot array:'
        call iwrtma(ipiv,1,nred,1,nred)
      end if

      if (cond.lt.1d-20)
     &     call quit(1,'leqc_core',
     &     'bad condition number, something must be wrong!')

      do iroot = 1, nroot
        idx = (iroot-1)*mxsub+1
        xvec(1:nred) = -gred(idx:idx-1+nred)

        job = 0
        call dgesl(xmat1,nred,nred,ipiv,xvec,job)

        if (ntest.ge.50) then
          write(luout,*) 'subspace solution for root # ',iroot
          call wrtmat2(xvec,nred,1,nred,1)
        end if

        vred(idx:idx-1+nred) = xvec(1:nred)
      end do

      irecscr = 1
      do iroot = 1, nroot
        ! assemble residual in full space
        do iopt = 1, nopt
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
            ffspc => me_special(2)%mel%fhand
            trafo = .true.
          else
            ffspc => ffscr(iopt)%fhand
            trafo = .false.
          end if

          if (nincore.ge.2) then
            call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,
     &                       nwfpar(iopt))
          else
            call da_sccpvec(me_rhs(iopt)%mel%fhand,iroot,
     &                      ffspc,iroot,
     &                      1d0,nwfpar(iopt),xbuf1,lenbuf)
          end if

          idx = (iroot-1)*mxsub + 1
          call optc_expand_vec(vred(idx),ndim_rsbsp,
     &                                        xrsnrm(iroot,iopt),.true.,
     &         ffspc,irecscr,1d0,ffrsbsp(iopt)%fhand,
     &         iord_rsbsp,
     &         nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          ! care for shift (incl. metric, if applicable)
          if (opti_info%shift.ne.0d0 .and. .not.use_s(iopt))
     &      call optc_expand_vec(
     &           opti_info%shift*vred(idx:idx-1+ndim_vsbsp),ndim_vsbsp,
     &                     xrsnrm(iroot,iopt),.true.,
     &           ffspc,irecscr,1d0,ffvsbsp(iopt)%fhand,
     &           iord_rsbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          if (opti_info%shift.ne.0d0 .and. use_s(iopt))
     &      call optc_expand_vec(
     &           opti_info%shift*vred(idx:idx-1+ndim_vsbsp),ndim_ssbsp,
     &                     xrsnrm(iroot,iopt),.true.,
     &           ffscr(iopt)%fhand,irecscr,1d0,ffssbsp(iopt)%fhand,
     &           iord_rsbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

C ?????
c          ! fix the signs (if needed)
c          if (.not.use_s(iopt).and.trafo)
c     &         call optc_fix_signs(xrsnrm(iroot,iopt),xvec,
c     &                      ffvsbsp(iopt)%fhand,ndim_vsbsp,iord_vsbsp,
c     &                      ffspc(iopt)%fhand,irecscr,
c     &                      opti_info,iopt,
c     &                      nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          ! transform residual (if requested)
          if (trafo) then
            call optc_traf(me_scr(iopt)%mel,irecscr,xrsnrm(iroot,iopt),
     &                     me_special(2)%mel,irecscr,
     &                     fspc(1),'B',me_special,nspecial,
     &                     nwfpar(iopt),xbuf1,
     &                     orb_info,op_info,str_info,strmap_info)

          end if

        end do

        ! not yet converged? increase record counter
        conv = .true.
        do iopt = 1, nopt
          conv = conv.and.xrsnrm(iroot,iopt).lt.opti_info%thrgrd(iopt)
        end do
        if (.not.conv) then
          idxroot(irecscr) = iroot
          irecscr = irecscr+1 
        end if

      end do
      

      ! number of new directions
      nnew = irecscr-1
      ! if this is the last iteration: no new directions, but 
      ! assemble all vectors, converged or not
      if (iter.ge.opti_info%maxmacit) nnew = 0

      if (nnew.gt.0) then
        nsec_arr => opti_info%nsec(1:nopt)
        nsec = sum(nsec_arr)
        nwfpsec => opti_info%nwfpsec(1:nsec)
        idstsec => opti_info%idstsec(1:nsec)
        signsec => opti_info%signsec(1:nsec)!2(1:nsec)
c        signsec => opti_info%signsec2(1:nsec)

        ! reduced space exhausted?
        if (nred+nnew.gt.mxsub) then
          restart_mode = 0
          if (restart_mode.eq.0) then
            ! complete internal restart
            ! assemble orth. subspace exactly spanning the nroot 
            ! currently best solution vectors
            call optc_minspace(
     &           iord_vsbsp,ffvsbsp,
     &           iord_rsbsp,ffrsbsp,
     &           iord_ssbsp,ffssbsp,use_s,
     &           vred,gred,mred,sred,nred,nroot,nroot,mxsub,nopt,
     &           ffscr(1)%fhand,nnew,  ! only scratch
     &           nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
            ndim_vsbsp = nred
            ndim_rsbsp = nred
            ndim_ssbsp = nred
          else
            call quit(1,'leqc_core','baustelle')
          end if

        end if

        ! divide new directions by preconditioner
        stsec = 1
        ndsec = 0
        do iopt = 1, nopt
          if (iopt.gt.1) stsec = stsec + nsec_arr(iopt-1)
          ndsec = ndsec + nsec_arr(iopt)

          select case(opti_info%typ_prc(iopt))
          case(optinf_prc_file,optinf_prc_traf)
            if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
              ffspc => me_special(2)%mel%fhand
              trafo = .true.
            else
              ffspc => ffscr(iopt)%fhand
              trafo = .false.
            end if
            if (nincore.ge.2) then
              call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,
     &                         nwfpar(iopt))
              do iroot = 1, nnew
                call vec_from_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                         nwfpar(iopt))
                ! scale residual for numerical stability:
c                xnrm = dnrm2(nwfpar,xbuf1,1)
                xnrm = 0d0
                do jopt = 1, nopt
                  xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
                end do
                xnrm = sqrt(xnrm)
                ! account for sign changes if necessary
                do isec = stsec, ndsec
                  call diavc(xbuf1(idstsec(isec)),xbuf1(idstsec(isec)),
     &                       signsec(isec)/xnrm,xbuf2(idstsec(isec)),
     &                       opti_info%shift,nwfpsec(isec))
                end do
c                call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,
c     &                     opti_info%shift,nwfpar(iopt))
                call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                         nwfpar(iopt))
              end do
            else
              ! security trap
              do isec = stsec, ndsec
                if (idstsec(isec).ne.1) 
     &               call quit(1,'leqc_core',
     &               'not this route (only incore)')
              end do

              do iroot = 1, nnew
c                ! request (nroot-iroot+1)th-last root 
c                irec = ioptc_get_sbsp_rec(-nroot+iroot-1,
c     &               iord_vsbsp,ndim_vsbsp,mxsbsp)
                xnrm = 0d0
                do jopt = 1, nopt
                  xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
                end do
                xnrm = sqrt(xnrm)
                call da_diavec(ffscr(iopt)%fhand,iroot,1,0d0,
     &                         ffscr(iopt)%fhand,iroot,1,1d0/xnrm,
     &                          me_dia(iopt)%mel%fhand,1,1,
     &                          opti_info%shift,-1d0,
     &                          nwfpar(iopt),xbuf1,xbuf2,lenbuf)
              end do
            end if
          case(optinf_prc_blocked)
            if (nincore.lt.3)
     &           call quit(1,'leqc_core',
     &           'I need at least 3 incore vectors (prc_blocked)')
            do iroot = 1, nnew
              call vec_from_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                         nwfpar(iopt))
              xnrm = 0d0
              do jopt = 1, nopt
                xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
              end do
              xnrm = sqrt(xnrm)
              call dscal(nwfpar(iopt),1d0/xnrm,xbuf1,1)
              call optc_prc_special2(me_mvp(iopt)%mel,me_special,
     &                                                        nspecial,
     &                           me_opt(iopt)%mel%op%name,
     &                           opti_info%shift,
     &                           nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                           orb_info,op_info,str_info,strmap_info)
              call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,nwfpar(iopt))
            end do
          case(optinf_prc_mixed)
            if (nincore.lt.3)
     &           call quit(1,'leqc_core',
     &           'I need at least 3 incore vectors (prc_mixed)')
            do iroot = 1, nnew
              call vec_from_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                         nwfpar(iopt))
              call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,
     &                         nwfpar(iopt))
              xnrm = 0d0
              do jopt = 1, nopt
                xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
              end do
              xnrm = sqrt(xnrm)
              call dscal(nwfpar(iopt),1d0/xnrm,xbuf1,1)
              call optc_prc_mixed(me_mvp(iopt)%mel,me_special,
     &                                                        nspecial,
     &                           me_opt(iopt)%mel%op%name,
     &                           opti_info%shift,
     &                           nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                           orb_info,op_info,str_info,strmap_info)
              call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,nwfpar(iopt))
            end do
          case default
            call quit(1,'leqc_core','unknown preconditioner type')
          end select

          ! if requested, transform new subspace vectors
          if (trafo) then
            do iroot = 1, nnew
              call optc_traf(me_scr(iopt)%mel,iroot,xdum,
     &                     me_special(2)%mel,iroot,
     &                     fspc(1),'F',me_special,nspecial,
     &                     nwfpar(iopt),xbuf1,
     &                     orb_info,op_info,str_info,strmap_info)
            end do
          end if

        end do ! iopt

        do iroot = 1, nnew
          getnewrec = .true.
          do iopt = 1, nopt
            if (use_s(iopt)) then
              ! assign op. with list containing the scratch trial vector
              call assign_me_list(me_scr(iopt)%mel%label,
     &                            me_opt(iopt)%mel%op%name,op_info)

              ! calculate metric * scratch trial vector
              allocate(xret(depend%ntargets),idxselect(depend%ntargets))
              nselect = 0
              call select_formula_target(idxselect,nselect,
     &                    me_met(iopt)%mel%label,depend,op_info)
              if (getnewrec) then
                irec = ioptc_get_sbsp_rec(0,iord_ssbsp,ndim_ssbsp,mxsub)
                if (iroot.eq.1) ioff_s = irec-1
                getnewrec = .false.
              end if
              call switch_mel_record(me_met(iopt)%mel,irec)
              call switch_mel_record(me_scr(iopt)%mel,iroot)
              call frm_sched(xret,flist,depend,idxselect,nselect,
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)
              me_met(iopt)%mel%fhand%last_mod(irec) = -1
              deallocate(xret,idxselect)

              ! reassign op. with list containing trial vector
              call assign_me_list(me_trv(iopt)%mel%label,
     &                            me_opt(iopt)%mel%op%name,op_info)
              ffmet(iopt)%fhand => me_met(iopt)%mel%fhand
            else
              ffmet(iopt)%fhand => fdum
            end if
          end do
        end do

        ! orthogonalize new directions to existing subspace
        ! and add linear independent ones to subspace
        nsec_arr => opti_info%nsec(1:nopt)
        nsec = sum(nsec_arr)
        nwfpsec => opti_info%nwfpsec(1:nsec)
        idstsec => opti_info%idstsec(1:nsec)
        signsec => opti_info%signsec(1:nsec)
        call optc_orthvec(nadd,.false.,
     &                 ffssbsp,iord_ssbsp,sred,
     &                 ffvsbsp,
     &                 iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                 use_s,ioff_s,ffmet,ffscr,nnew,nopt,
     &                 nsec_arr,nwfpsec,idstsec,signsec,
     &                 nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)

        ! set nadd
        if (nadd.eq.0)
     &       call quit(0,'leqc_core',
     &       'solver in problems: only linear dependent '//
     &       'new directions?')
        opti_stat%nadd = nadd

        ! |Mv> subspace organisation should be identical to |v> subsp.
        ndim_rsbsp = ndim_vsbsp
        iord_rsbsp = iord_vsbsp
        ! dto. for |Sv> subspace
        ndim_ssbsp = ndim_vsbsp
        iord_ssbsp = iord_vsbsp

      else
        ! if all converged: assemble vectors 
        do iopt = 1, nopt
          do iroot = 1, nroot

            idx = (iroot-1)*mxsub + 1
            call optc_expand_vec(vred(idx),ndim_vsbsp,xdum,.false.,
     &           me_opt(iopt)%mel%fhand,iroot,0d0,ffvsbsp(iopt)%fhand,
     &           iord_vsbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          end do
        end do

      end if

      deallocate(ffmet)

      return
      end

