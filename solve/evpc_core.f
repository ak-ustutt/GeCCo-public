*----------------------------------------------------------------------*
      subroutine evpc_core(iter,
     &       task,iroute,xrsnrm,xeig,
     &       use_s,
     &       me_opt,me_trv,me_mvp,me_dia,me_met,me_scr,me_ext,
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
     &     ntest = 000

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
     &     me_mvp(*), me_special(*), me_scr(*),me_ext(*)
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
     &     idx, jdx, kdx, iroot, nred, nadd, nnew, irecscr,idxdbg,
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
     &     ipiv(:), idxroot(:)
      type(file_array), pointer ::
     &     ffrsbsp(:), ffvsbsp(:), ffssbsp(:), ffscr(:), ffmet(:),
     &     ffext(:)
      type(filinf) ::
     &     fdum
      type(filinf), pointer ::
     &     ffspc
      type(filinf), target ::
     &     fdum2

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2, ddot, xnormop

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'evpc_core entered')
      if (ntest.ge.100) write (lulog,*) "iteration ",iter
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
      ffext => opti_stat%ffext
      ffrsbsp => opti_stat%ffrsbsp
      ffvsbsp => opti_stat%ffvsbsp
      ffssbsp => opti_stat%ffssbsp
      nwfpar => opti_info%nwfpar

      allocate(ffmet(nopt))



      ifree = mem_alloc_int(idxroot,nopt*nroot,'EVP_idxroot')


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
!         if (nopt.ne.1) call quit(1,'evpc_core','not this route')
          call optc_update_redsp3
     &       (mred,xdum,nred,0,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,init,
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
      call rl8_compress_array(mred, nred, mxsub, xmat1)!mred -> xmat1
      call rl8_compress_array(sred, nred, mxsub, xmat3)!sred -> xmat3

      if (.not.use_s(1)) then !! @TODO results should be invariant with respect to exchange of optimized variables
        call eigen_asym(nred,xmat1,eigr,eigi,xmat2,xmat3,ierr)
      else
        call eigen_asym_met(nred,xmat1,xmat3,eigr,eigi,xmat2,xmat3,ierr)
      end if

      ! copy back to vred with mxsub as leading dim
      call rl8_expand_array(xmat2,vred,nred,mxsub)!xmat2 -> vred


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
          elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
             if (nspecial .lt. 4) call quit(1,'evpc_core',
     &            'TR0 with <4 me_lists')
             ffspc => me_special(4)%mel%fhand
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


            call trafo_forward_wrap(flist,depend,
     &            me_special,me_scr,me_trv,
     &            xrsnrm, nroot, 
     &            iroot, iopt, irecscr,
     &            op_info, str_info, strmap_info, orb_info)


!     set all single excitations to zero if requestes
!     after long deliberation, I decided to also include V,V so one can be sure to include
!     **all** singular excitations even if T1 changes
            if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
               call set_blks(me_scr(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
               xrsnrm(iroot,iopt)=xnormop(me_scr(iopt)%mel) 
            endif

c dbg
C            if (iopt .ne. 2) then
c            print *,'residual vector before transformation:'
c            call vec_from_da(ffspc,
c     &        irecscr,xbuf1,nwfpar(iopt))
c            do idx = 1, nwfpar(iopt)
c              print *,idx,xbuf1(idx)
c            end do
c            call print_list('transformed residual vector:',
c     &           me_scr(iopt)%mel,"LIST",
c     &           -1d0,0d0,
c     &           orb_info,str_info)
c            endif
c dbgend

          end if
        end do

        ! not yet converged? increase record counter
        conv = .true.
        do iopt = 1, nopt
          conv = conv.and.xrsnrm(iroot,iopt).lt.opti_info%thrgrd(iopt)
        end do
        if (.not.conv.and.iter.lt.opti_info%maxmacit) then
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
     &           ffext,0,  ! only scratch
     &           nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
            ndim_vsbsp = nred
            ndim_rsbsp = nred
            ndim_ssbsp = nred
          else
!> @TODO implement evpc_core with restart possibility.
            call quit(1,'evpc_core','baustelle')
          end if

        end if

        ! divide new directions by preconditioner
        do iopt = 1, nopt
          select case(opti_info%typ_prc(iopt))
          case(optinf_prc_file,optinf_prc_traf,optinf_prc_traf_spc
     &         ,optinf_prc_spinp,optinf_prc_prj,optinf_prc_spinrefp)
            if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
              ffspc => me_special(1)%mel%fhand
              trafo = .true.
            elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
               ffspc => me_special(4)%mel%fhand
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
c dbg
c                print *,"xnorm:",xnrm
c dbgend                
c                xnrm = 1d0
c dbg
c                print *,"precon not yet applied. "
c                do idxdbg = 1, nwfpar(iopt)
c                   print *,idxdbg,xbuf1(idxdbg)
c                end do
c dbgend
                xshf = -xeig(idxroot(iroot),1) 
                call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,xshf,nwfpar(iopt))
c dbg
c                print *,"debug shift:",xshf
c                print *,"precon applied. ",nwfpar(iopt),"entries"
c                do idxdbg = 1, nwfpar(iopt)
c                   print *,idxdbg,xbuf1(idxdbg),xbuf2(idxdbg)
c                end do
c dbgend
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
                ! decrease xshf in first iteration
                if (iter .eq. 1) xshf=0.8d0*xshf
                call da_diavec(ffspc,iroot,0d0,
     &                   ffscr(iopt)%fhand,iroot,
     &                   1d0/xnrm,me_dia(iopt)%mel%fhand,
     &                   1,xshf,-1d0,
     &                   nwfpar(iopt),xbuf1,xbuf2,lenbuf)
              end do
            end if

            ! project out spin contaminations or other components?
            if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &          opti_info%typ_prc(iopt).eq.optinf_prc_prj.or.
     &          opti_info%typ_prc(iopt).eq.optinf_prc_spinrefp) then
              ! assign op. with list containing the scratch trial vector
c dbg
c               print *, "assign ",me_scr(iopt)%mel%label," to ",
c     &              me_opt(iopt)%mel%op%name
c dbgend
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
                elseif (opti_info%typ_prc(iopt).eq.
     &                  optinf_prc_spinrefp)then
                  call spin_project(me_scr(iopt)%mel,me_special(1)%mel,
     &                              fspc(2),opti_info%nwfpar(iopt),
     &                              xbuf1,xbuf2,.true.,xnrm,
     &                              opti_info,orb_info,
     &                              op_info,str_info,strmap_info)

                  call evaluate2(fspc(1),.false.,.false.,
     &                           op_info,str_info,strmap_info,orb_info,
     &                           xnrm,.false.)

                else
c dbg
c                  print *,"iopt=",iopt
c                  call print_list('projected vector:',
c     &                 me_scr(iopt)%mel,"NORM",
c     &                 -1d0,0d0,
c     &                 orb_info,str_info)
c dbgend
                  call evaluate2(fspc(1),.false.,.false.,
     &                 op_info,str_info,strmap_info,orb_info,
     &                 xnrm,.false.)
c dbg
c                  call print_list('projected vector:',
c     &                 me_scr(iopt)%mel,"NORM",
c     &                 -1d0,0d0,
c     &                 orb_info,str_info)
c dbgend
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
            if (nspecial.ge.3)
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
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                  call  switch_mel_record(me_special(4)%mel,iroot)
               else
                  call switch_mel_record(me_special(1)%mel,iroot)
               endif

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

c dbg 
c              call print_list("scratch trial vec", me_scr(iopt)%mel,
c     &             "NORM",0,0,orb_info,str_info)
c dbgend
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
        ! if all converged or last iteration: assemble vectors 

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

      contains
*----------------------------------------------------------------------*
!>   compresses the the occupied parts of an array into a another array
!!
!!   assumes the input array is a transformed matrix where only the upper left
!!   quadrant is occupied
!!   nmat gives dimension of the quadrant lmat of the actual matrix
*----------------------------------------------------------------------*
      subroutine rl8_compress_array(arr_in, nmat, lmat, arr_out)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="rl8_copy_matrix_to_array"
      integer,parameter::
     &     ntest=00

      real(8),dimension(*),intent(in)::
     &     arr_in
      real(8),dimension(*),intent(out)::
     &     arr_out
      integer,intent(in)::
     &     lmat,                !> actual dimension of matrix
     &     nmat                 !> dimension of filled part
      
      integer::
     &     kdx, idx, jdx

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)

      kdx=0
      do idx = 1, nmat
         do jdx = 1, nmat
            kdx=kdx+1
            arr_out(kdx) = arr_in((idx-1)*lmat+jdx)
         end do
      end do
      return
      end subroutine
*----------------------------------------------------------------------*
!>  expands an array into a another array
!!
!!   assumes the output array is a transformed matrix where only the upper left
!!   quadrant is occupied
!!   nmat gives dimension of the quadrant lmat of the actual matrix
!!   note the changed order of arguments as nmat and lmat belong to the matrix
*----------------------------------------------------------------------*
      subroutine rl8_expand_array(arr_in, arr_out, nmat, lmat )
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="rl8_copy_matrix_to_array"
      integer,parameter::
     &     ntest=00

      real(8),dimension(*),intent(in)::
     &     arr_in
      real(8),dimension(*),intent(out)::
     &     arr_out
      integer,intent(in)::
     &     lmat,                !> actual dimension of matrix
     &     nmat                 !> dimension of filled part
      
      integer::
     &     kdx, idx, jdx

      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)

      kdx=0
      do idx = 1, nmat
         do jdx = 1, nmat
            kdx=kdx+1
            arr_out((idx-1)*lmat+jdx) = arr_in(kdx)
         end do
      end do
      return
      end subroutine


*----------------------------------------------------------------------*
!>    wrapper to encapsulate some stupid decisions
!!
!!
*----------------------------------------------------------------------*
      subroutine trafo_forward_wrap(flist,depend,
     &     me_special,me_scr,me_trv,
     &     xrsnrm, 
     &     nroot, iroot, iopt, irecscr,
     &     op_info, str_info, strmap_info, orb_info)
*----------------------------------------------------------------------*
      implicit none

      type(me_list_array), dimension(*)::
     &     me_special,me_scr,me_trv
      type(formula_item),intent(in)::
     &     flist
      type(dependency_info),intent(in)::
     &     depend
      integer, intent(in)::
     &     nroot,
     &     iroot, 
     &     iopt,
     &     irecscr

      real(8), Dimension(nroot,*), intent(inout)::
     &     xrsnrm

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info


      type(me_list),pointer::
     &     me_in,
     &     me_trf
      type(operator),pointer::
     &     op_in,
     &     op_trf
      real(8) ::
     &     xnrm

      if (nspecial.ge.3)then ! who thought it would be good 
         me_trf=> me_special(3)%mel
         op_trf=> me_special(3)%mel%op
      else
         me_trf=> null() ! leave it associated as it is now
         op_trf=> null() ! --
      end if

      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
         me_in => me_special(4)%mel
      else
         me_in => me_special(1)%mel
      endif


      call switch_mel_record(me_scr(iopt)%mel,irecscr)
      call  switch_mel_record(me_in,irecscr)

      call change_basis_old(flist, depend,
     &     me_in, me_in%op,
     &     me_scr(iopt)%mel, me_scr(iopt)%mel%op, xnrm,
     &     me_trf, op_trf,                         ! notice the difference : me_trf != me_trv
     &     me_trv(iopt)%mel,
     &     op_info, str_info, strmap_info, orb_info)

      xrsnrm(iroot,iopt) = xnrm
      return
      end subroutine





*----------------------------------------------------------------------*
!>    subroutine for the transformation into the orthogonal basis
!!
!!    uses the old convention where the transformation formula is a subset of the whole formula
!!    me_tgt and determines the me_list that was bound to the target operator as the dependencies where 
!!    evaluated
*----------------------------------------------------------------------*
      subroutine change_basis_old(flist, depend,
     &     me_in, op_in, 
     &     me_out, op_out, outnrm,
     &     me_trf, op_trf,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="change_basis_old"
      integer,parameter::
     &     ntest=100
      
      type(formula_item)::
     &     flist

      type(me_list)::
     &     me_in,              !> list to be transformed 
     &     me_out,             !> result list
     &     me_tgt             !>list with the transformation operator
      type(me_list),pointer::     
     &     me_trf               !> target list definition to specify which subformula should actually be evaluated (stupid design decision)
      type(operator)::
     &     op_in,
     &     op_out
      type(operator),pointer::
     &     op_trf
      real(8),intent(out)::
     &     outnrm               !norm of the output list
      type(dependency_info)::
     &     depend               !>dependency info for the formula 
      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

      integer,dimension(:),allocatable::
     &     idxselect
      real(8),dimension(:),allocatable::
     &     xret
      integer::
     &     nselect

      if (ntest.ge.100)then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) "out:",me_out%label,op_out%name
         end if 
      call assign_me_list(me_out%label,
     &     op_out%name, op_info)

      if(associated(me_trf).and. associated(op_trf)) then
         call assign_me_list(me_trf%label, op_trf%name, op_info)
      else if(associated(me_trf).or. associated(op_trf))then
         call quit(1,i_am,
     &        "please make sure that either both (transformation list"//
     &        "and operator) or neither is associated") 
      end if
      allocate(xret(depend%ntargets),idxselect(depend%ntargets))
      nselect=0
      call select_formula_target(idxselect,nselect,
     &     me_tgt%label,depend,op_info)
! pretend that me_trv is not up to date
      call reset_file_rec(me_tgt%fhand)
      call frm_sched(xret,flist,depend,idxselect, nselect,
     &     .true.,.false.,op_info,str_info,strmap_info,orb_info)
      outnrm=xret(idxselect(1))
      deallocate(xret,idxselect)
      return
      end subroutine
*----------------------------------------------------------------------*
      end

