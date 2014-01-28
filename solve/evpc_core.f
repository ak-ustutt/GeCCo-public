*----------------------------------------------------------------------*
      subroutine evpc_core(iter,
     &       task,iroute,xrsnrm,xeig,
     &       use_s,
     &       me_opt,me_trv,me_mvp,me_dia,me_met,me_scr,
     &       me_special,nspecial,
c     &       ffopt,fftrv,ffmvp,ffdia,
     &       nincore,lenbuf,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,
     &       fspc,nspcfrm,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     core driver for EVP solver
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

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_dia(*),
     &     me_mvp(*), me_special(*), me_scr(*)
      type(me_list_array), intent(inout) ::
     &     me_met(*), me_trv(*)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffdia(*)

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
     &     xrsnrm(opti_info%nroot,opti_info%nopt),
     &     xeig(opti_info%nroot,2)

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
     &     idx, jdx, kdx, iroot, nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job,
     &     ndim_save, ndel, iopt, jopt, lenscr, ioff,
     &     ifree, restart_mode, ierr, nselect, irec, ioff_s
      real(8) ::
     &     cond, xdum, xnrm, xshf
      real(8), pointer ::
     &     gred(:), vred(:), mred(:), sred(:), eigr(:), eigi(:),
     &     xmat1(:), xmat2(:), xmat3(:), xvec(:), xret(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &     iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:),
     &     nwfpar(:), idxselect(:),
     &     ipiv(:), iconv(:), idxroot(:)
      type(file_array), pointer ::
     &     ffrsbsp(:), ffvsbsp(:), ffssbsp(:), ffscr(:), ffmet(:)
      type(filinf) ::
     &     fdum
      type(filinf), pointer ::
     &     ffspc
      type(filinf), target ::
     &     fdum2

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2, ddot

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'evpc_core entered')

      zero_vec(1:opti_stat%ndim_vsbsp) = .false.
      nopt = opti_info%nopt
      nroot = opti_info%nroot
      mxsub = opti_stat%mxdim_sbsp
      mred => opti_stat%sbspmat(1:)
      vred => opti_stat%sbspmat(mxsub**2+1:)
      sred => opti_stat%sbspmat(2*mxsub**2+1:)
      ndim_rsbsp => opti_stat%ndim_rsbsp
      ndim_vsbsp => opti_stat%ndim_vsbsp
      ndim_ssbsp => opti_stat%ndim_ssbsp
      iord_rsbsp => opti_stat%iord_rsbsp
      iord_vsbsp => opti_stat%iord_vsbsp
      iord_ssbsp => opti_stat%iord_ssbsp
      ffscr => opti_stat%ffscr
      ffrsbsp => opti_stat%ffrsbsp
      ffvsbsp => opti_stat%ffvsbsp
      ffssbsp => opti_stat%ffssbsp
      nwfpar => opti_info%nwfpar

      allocate(ffmet(nopt))

c      if (nopt.gt.1)
c     &     call quit(1,'evpc_core','not yet adapted for nopt>1')

      ! check for previously converged roots
      ifree = mem_alloc_int(iconv,nopt*nroot,'EVP_conv')
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

      ifree = mem_alloc_int(idxroot,nopt*nroot,'EVP_idxroot')

c      iopt = 1  ! preliminary
      if (ndim_vsbsp.ne.ndim_rsbsp)
     &     call quit(1,'evpc_core','subspace dimensions differ?')
c dbg
c      call optc_check_redsp(
c     &     ndim_vsbsp,ndim_vsbsp,nopt,
c     &     iord_vsbsp,ffvsbsp,
c     &     nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
c dbg
      nred = ndim_vsbsp
      do iopt = 1, nopt
        init = iopt.eq.1
        ! update reduced space:
        ! ffvsbsp and ffrsbsp point to ff_trv(iopt)%fhand ...
        if (.not.use_s(iopt)) then
          if (nopt.ne.1) call quit(1,'evpc_core','not this route')
          call optc_update_redsp3
     &       (mred,xdum,nred,0,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
     &       iord_rsbsp,ffrsbsp(iopt)%fhand,fdum,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        else
          call optc_update_redsp4
     &       (mred,sred,xdum,nred,0,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,init,
     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
     &       iord_rsbsp,ffrsbsp(iopt)%fhand,
     &       iord_ssbsp,ffssbsp(iopt)%fhand,fdum,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        end if
      end do ! iopt
      ! ------------------------
      !    solve reduced EVP
      ! ------------------------ 
      
      ! allocate some scratch 
      ! (automatically deallocated after leaving leq_evp_control() )
      lenmat = nred*nred
      ifree = mem_alloc_real(xmat1,lenmat,'EVP_mat1')
      ifree = mem_alloc_real(xmat2,lenmat,'EVP_mat2')
      ifree = mem_alloc_real(xmat3,lenmat,'EVP_mat3')
      ifree = mem_alloc_real(eigr,nred,'EVP_eig_r')
      ifree = mem_alloc_real(eigi,nred,'EVP_eig_i')
      ifree = mem_alloc_real(xvec,nred,'EVP_vec')

      ! get a copy of the subspace matrix
      kdx = 0
      do idx = 1, nred
        do jdx = 1, nred
          kdx = kdx+1
          xmat1(kdx) = mred((idx-1)*mxsub+jdx)
        end do
      end do
      kdx = 0
      do idx = 1, nred
        do jdx = 1, nred
          kdx = kdx+1
          xmat3(kdx) = sred((idx-1)*mxsub+jdx)
        end do
      end do

      if (.not.use_s(1)) then
        call eigen_asym(nred,xmat1,eigr,eigi,xmat2,xmat3,ierr)
      else
        call eigen_asym_met(nred,xmat1,xmat3,eigr,eigi,xmat2,xmat3,ierr)
      end if

      ! copy back to vred with mxsub as leading dim
      kdx = 0
      do idx = 1, nred
        do jdx = 1, nred
          kdx = kdx+1
          vred((idx-1)*mxsub+jdx) = xmat2(kdx) 
        end do
      end do

      if (ntest.ge.50) then
        write(lulog,*) 'Eigenvalues in subspace:'
        do idx = 1, nred
          write(lulog,'(2x,i4,x,g20.12,x,g20.12)')
     &         idx, eigr(idx), eigi(idx)
        end do
        if (ntest.ge.100) then
          write(lulog,*) 'eigenvectors:'
          call wrtmat2(vred,nred,nred,mxsub,mxsub)
        end if
      end if

      xeig(1:nroot,1) = eigr(1:nroot)
      xeig(1:nroot,2) = eigi(1:nroot)

      irecscr = 1
      do iroot = 1, nroot
        ! assemble residual in full space
        
        do iopt = 1, nopt
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
            ffspc => me_special(1)%mel%fhand
            trafo = .true.
          else
            ffspc => ffscr(iopt)%fhand
            trafo = .false.
          end if
          ! M.v ....
          idx = (iroot-1)*mxsub + 1
          xvec(1:nred) = vred(idx:idx+nred-1)
          call optc_expand_vec(xvec,ndim_rsbsp,
     &                                     xrsnrm(iroot,iopt),.false.,
     &         ffspc,irecscr,0d0,ffrsbsp(iopt)%fhand,
     &         iord_rsbsp,
     &         nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          xvec(1:nred) = -eigr(iroot)*xvec(1:nred)
          if (.not.use_s(iopt)) then
          ! - eig * v
            call optc_expand_vec(xvec,ndim_vsbsp,
     &                                    xrsnrm(iroot,iopt),.not.trafo,
     &       ffspc,irecscr,1d0,ffvsbsp(iopt)%fhand,
     &       iord_vsbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
          else
            ! - eig * S * v
            call optc_expand_vec(xvec,ndim_vsbsp,
     &                                    xrsnrm(iroot,iopt),.not.trafo,
     &           ffspc,irecscr,1d0,ffssbsp(iopt)%fhand,
     &           iord_ssbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)
          end if

          ! if requested, transform residual
          if (trafo) then
            ! assign op. with list containing the scratch trial vector
            call assign_me_list(me_scr(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)
            ! use daggered transformation matrix if requested
            if (nspecial.eq.3)
     &         call assign_me_list(me_special(3)%mel%label,
     &                             me_special(3)%mel%op%name,op_info)

            ! calculate transformed residual
            allocate(xret(depend%ntargets),idxselect(depend%ntargets))
            nselect = 0
            call select_formula_target(idxselect,nselect,
     &                  me_trv(iopt)%mel%label,depend,op_info)
            call switch_mel_record(me_special(1)%mel,irecscr)
            call switch_mel_record(me_scr(iopt)%mel,irecscr)
            ! pretend that me_trv is not up to date
            call reset_file_rec(me_trv(iopt)%mel%fhand)
            call frm_sched(xret,flist,depend,idxselect,nselect,
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)
            ! residual norm (nselect should be 1):
            xrsnrm(iroot,iopt) = xret(idxselect(1))
            deallocate(xret,idxselect)
c dbg
c            print *,'residual vector before transformation:'
c            call vec_from_da(ffspc,
c     &        irecscr,xbuf1,nwfpar(iopt))
c            do idx = 1, nwfpar(iopt)
c              print *,idx,xbuf1(idx)
c            end do
c            print *,'transformed residual vector:'
c            call vec_from_da(me_scr(iopt)%mel%fhand,
c     &        me_scr(iopt)%mel%fhand%current_record,xbuf1,nwfpar(iopt))
c            do idx = 1, nwfpar(iopt)
c              print *,idx,xbuf1(idx)
c            end do
c dbgend

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
      if (nnew.gt.0) then

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
     &           vred,xdum,mred,sred,nred,nroot,0,mxsub,nopt,
     &           ffscr(1)%fhand,nnew,  ! only scratch
     &           nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
            ndim_vsbsp = nred
            ndim_rsbsp = nred
            ndim_ssbsp = nred
          else
            call quit(1,'evpc_core','baustelle')
          end if

        end if

        ! divide new directions by preconditioner
        do iopt = 1, nopt

          select case(opti_info%typ_prc(iopt))
          case(optinf_prc_file,optinf_prc_traf,optinf_prc_spinp,
     &         optinf_prc_prj)
            if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
              ffspc => me_special(1)%mel%fhand
              trafo = .true.
            else
              ffspc => ffscr(iopt)%fhand
              trafo = .false.
            end if
            if (nincore.ge.2) then
              call vec_from_da(
     &             me_dia(iopt)%mel%fhand,1,xbuf2,nwfpar(iopt))
              do iroot = 1, nnew
                call vec_from_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                           nwfpar(iopt))
                ! scale residual for numerical stability:
                xnrm = 0d0
                do jopt = 1, nopt
                  xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
                end do
                xnrm = sqrt(xnrm)
c                xnrm = 1d0
                xshf = -xeig(idxroot(iroot),1)
                call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,xshf,nwfpar(iopt))
                if (nopt.eq.1) then
                  xnrm = dnrm2(nwfpar(iopt),xbuf1,1)
                  call dscal(nwfpar(iopt),1d0/xnrm,xbuf1,1)
                end if
                call vec_to_da(ffspc,iroot,xbuf1,
     &                         nwfpar(iopt))
              end do
            else
              do iroot = 1, nnew
c            ! request (nroot-iroot+1)th-last root 
c            irec = ioptc_get_sbsp_rec(-nroot+iroot-1,
c     &         iord_vsbsp,ndim_vsbsp,mxsbsp)
                xnrm = 0d0
                do jopt = 1, nopt
                  xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
                end do
                xnrm = sqrt(xnrm)
                xshf = -xeig(idxroot(iroot),1)
                call da_diavec(ffspc,iroot,0d0,
     &                   ffscr(iopt)%fhand,iroot,
     &                   1d0/xnrm,me_dia(iopt)%mel%fhand,
     &                   1,xshf,-1d0,
     &                   nwfpar(iopt),xbuf1,xbuf2,lenbuf)
              end do
            end if

            ! project out spin contaminations or other components?
            if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &          opti_info%typ_prc(iopt).eq.optinf_prc_prj) then      
              ! assign op. with list containing the scratch trial vector
              call assign_me_list(me_scr(iopt)%mel%label,
     &                            me_opt(iopt)%mel%op%name,op_info)
              do iroot = 1, nnew
                call switch_mel_record(me_scr(iopt)%mel,iroot)
                if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp) then
                  call spin_project(me_scr(iopt)%mel,me_special(1)%mel,
     &                              fspc(1),opti_info%nwfpar(iopt),
     &                              xbuf1,xbuf2,.true.,xnrm,
     &                              opti_info,orb_info,
     &                              op_info,str_info,strmap_info)
                else
                  call evaluate2(fspc(1),.false.,.false.,
     &                           op_info,str_info,strmap_info,orb_info,
     &                           xnrm,.false.)
                end if
                if (xnrm.lt.1d-12) call warn('evpc_core',
     &               'Nothing left after projection!')
              end do
              ! reassign op. with list containing trial vector
              call assign_me_list(me_trv(iopt)%mel%label,
     &                            me_opt(iopt)%mel%op%name,op_info)
            end if
          case(optinf_prc_blocked)
            if (nincore.lt.3)
     &           call quit(1,'evpc_core',
     &           'I need at least 3 incore vectors (prc_special)')
            do iroot = 1, nnew
              xnrm = 0d0
              do jopt = 1, nopt
                xnrm = xnrm+xrsnrm(idxroot(iroot),jopt)**2
              end do
              xnrm = sqrt(xnrm)
c              xnrm = 1d0
              call vec_from_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                         nwfpar(iopt))
              call dscal(nwfpar(iopt),1d0/xnrm,xbuf1,1)
              xshf = -xeig(idxroot(iroot),1)
              call optc_prc_special2(me_mvp(iopt)%mel,me_special,
     &                                                        nspecial,
     &                           me_opt(iopt)%mel%op%name,xshf,
     &                           nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                           orb_info,op_info,str_info,strmap_info)
              call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,
     &                       nwfpar(iopt))
            end do
          case default
            call quit(1,'evpc_core','unknown preconditioner type')
          end select

          ! if requested, transform new subspace vectors
          if (trafo) then
            ! use non-daggered transformation matrix if requested
            if (nspecial.eq.3)
     &         call assign_me_list(me_special(2)%mel%label,
     &                             me_special(2)%mel%op%name,op_info)
            ! assign op. with original list 
            ! (to ensure proper spin symmetry if needed)
            call assign_me_list(me_opt(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)

            ! calculate transformed vector
            allocate(xret(depend%ntargets),idxselect(depend%ntargets))
            nselect = 0
            call select_formula_target(idxselect,nselect,
     &                  me_trv(iopt)%mel%label,depend,op_info)
            do iroot = 1, nnew
              call switch_mel_record(me_special(1)%mel,iroot)
              call switch_mel_record(me_opt(iopt)%mel,iroot)
              ! pretend that me_trv is not up to date
              call reset_file_rec(me_trv(iopt)%mel%fhand)
              call frm_sched(xret,flist,depend,idxselect,nselect,
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)
              ! in reality me_trv is still up to date:
              call touch_file_rec(me_trv(iopt)%mel%fhand)
              ! copy to scr list and reassign
              call switch_mel_record(me_scr(iopt)%mel,iroot)
              call list_copy(me_opt(iopt)%mel,me_scr(iopt)%mel,.false.)
            end do
            call assign_me_list(me_scr(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)
            deallocate(xret,idxselect)
c dbg
c            print *,'back-transformed trial vector:'
c            call vec_from_da(me_scr(iopt)%mel%fhand,
c     &        me_scr(iopt)%mel%fhand%current_record,xbuf1,nwfpar(iopt))
c            do idx = 1, nwfpar(iopt)
c              print *,idx,xbuf1(idx)
c            end do
c            ! symmetrize!
c            if (nspecial.eq.3) then
c              call sym_ab_blk(xbuf2(2),xbuf1(2),0.5d0,dble(1),
c     &                        me_scr(iopt)%mel,2,
c     &                        str_info,strmap_info,orb_info)
c              xbuf2(1) = 0d0
c              print *,'back-transformed trial vector:'
c              do idx = 1, nwfpar(iopt)
c                print *,idx,xbuf2(idx)
c              end do
c              call vec_to_da(me_scr(iopt)%mel%fhand,
c     &         me_scr(iopt)%mel%fhand%current_record,xbuf2,nwfpar(iopt))
c            end if
c dbgend 

            ! reassign op. with list containing trial vector
            call assign_me_list(me_trv(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)
          end if

        end do ! iopt

        do iroot = 1, nnew
          getnewrec = .true.
          do iopt = 1, nopt

            if (use_s(iopt)) then
              ! assign op. with list containing the scratch trial vector
              call assign_me_list(me_scr(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)

              ! calculate metric * scratch trial vector
              allocate(xret(depend%ntargets),idxselect(depend%ntargets))
              nselect = 0
              call select_formula_target(idxselect,nselect,
     &                  me_met(iopt)%mel%label,depend,op_info)
              if (getnewrec) then
                irec = ioptc_get_sbsp_rec(0,iord_ssbsp,ndim_ssbsp,mxsub)
                if (iroot.eq.1) ioff_s = irec-1
                getnewrec = .false.
              end if
              call switch_mel_record(me_met(iopt)%mel,irec)
              call switch_mel_record(me_scr(iopt)%mel,iroot)
              call frm_sched(xret,flist,depend,idxselect,nselect,
     &            .true.,.false.,op_info,str_info,strmap_info,orb_info)
              ! apply sign-fix (if needed)
              call optc_fix_signs2(me_met(iopt)%mel%fhand,irec,
     &                             opti_info,iopt,
     &                             opti_info%nwfpar(iopt),xbuf1)
              call reset_file_rec(me_met(iopt)%mel%fhand)

              ! reassign op. with list containing trial vector
              call assign_me_list(me_trv(iopt)%mel%label,
     &                          me_opt(iopt)%mel%op%name,op_info)
              ffmet(iopt)%fhand => me_met(iopt)%mel%fhand
              deallocate(xret,idxselect)
            else
              ffmet(iopt)%fhand => fdum2
            end if

          end do ! iopt
        end do ! iroot

        ! orthogonalize new directions to existing subspace
        ! and add linear independent ones to subspace
        call optc_orthvec(nadd,nopt.gt.1,
     &                  ffssbsp,iord_ssbsp,sred,
     &                  ffvsbsp,iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                  use_s,ioff_s,ffmet,ffscr,nnew,nopt,
     &                  nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)
c dbg
c            print *,'orthogonalized trial vector:'
c            print *,'ndim= ',ndim_vsbsp!-1
c            call vec_from_da(ffvsbsp(1)%fhand,
cc     &        ndim_vsbsp-1,xbuf1,nwfpar(1))
c     &        ndim_vsbsp,xbuf1,nwfpar(1))
c            do idx = 1, nwfpar(1)
c              print *,idx,xbuf1(idx)
c            end do
c dbgend

        ! set nadd
        if (nadd.eq.0)
     &       call quit(0,'evpc_core',
     &       'solver in problems: only linear dependent '//
     &       'new directions?')
        opti_stat%nadd = nadd

        ! |Mv> subspace organisation should be identical to |v> subsp.
        ndim_rsbsp = ndim_vsbsp
        iord_rsbsp = iord_vsbsp
        ! dto. for |Sv> subspace ...
        ndim_ssbsp = ndim_vsbsp
        iord_ssbsp = iord_vsbsp

      else
        ! if all converged: assemble vectors 

        do iopt = 1, nopt
          do iroot = 1, nroot

            idx = (iroot-1)*mxsub + 1
            call optc_expand_vec(vred(idx),ndim_vsbsp,xdum,.false.,
     &           me_opt(iopt)%mel%fhand,iroot,0d0,
     &                              ffvsbsp(iopt)%fhand,iord_vsbsp,
     &           nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2)

          end do
        end do

      end if
    
      deallocate(ffmet)

      return
      end

