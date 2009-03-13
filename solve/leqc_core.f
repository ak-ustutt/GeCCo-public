*----------------------------------------------------------------------*
      subroutine leqc_core(iter,
     &       task,iroute,xrsnrm,
     &       use_s,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,
     &       me_special,nspecial,
c     &       ffopt,fftrv,ffmvp,ffrhs,ffdia,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       opti_info,opti_stat)
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

      integer, parameter ::
     &     ntest = 00

      integer, intent(inout) ::
     &     task
      real(8), intent(inout) ::
     &     xrsnrm(*)
      integer, intent(inout) ::
     &     iter
      integer, intent(in) ::
     &     iroute, nincore, lenbuf, nspecial
      logical, intent(in) ::
     &     use_s(*)

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_trv(*), me_dia(*),
     &     me_mvp(*), me_rhs(*),
     &     me_special(nspecial)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffrhs(*), ffdia(*)
      type(filinf), intent(in) ::
     &     ffscr

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

* local
      logical ::
     &     zero_vec(opti_stat%ndim_vsbsp)
      integer ::
     &     idx, jdx, kdx, iroot, irhs,  nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job,
     &     ndim_save, ndel, iopt, lenscr, ifree, restart_mode
      real(8) ::
     &     cond, xdum, xnrm
      real(8), pointer ::
     &     gred(:), vred(:), mred(:), sred(:),
     &     xmat1(:), xmat2(:), xvec(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &     iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:),
     &     nwfpar(:),
     &     ipiv(:), iconv(:), idxroot(:)
      type(filinf), pointer ::
     &     ffrsbsp, ffvsbsp, ffssbsp
      type(filinf) ::
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
      ffrsbsp => opti_stat%ffrsbsp(1)%fhand
      ffvsbsp => opti_stat%ffvsbsp(1)%fhand
      ffssbsp => opti_stat%ffssbsp(1)%fhand
      nwfpar => opti_info%nwfpar

      if (nopt.gt.1)
     &     call quit(1,'leqc_core','not yet adapted for nopt>1')

      ! check for previously converged roots
      ifree = mem_alloc_int(iconv,nopt*nroot,'LEQ_conv')
      idx = 0
      do iopt = 1, nopt
        do iroot = 1, nroot
          idx = idx+1
          iconv(idx) = 0
          if (iter.gt.1.and.xrsnrm(idx).lt.opti_info%thrgrd(iopt))
     &         iconv(idx) = 1
        end do
      end do

      ifree = mem_alloc_int(idxroot,nopt*nroot,'LEQ_idxroot')

      iopt = 1  ! preliminary
      if (ndim_vsbsp.ne.ndim_rsbsp)
     &     call quit(1,'leqc_core','subspace dimensions differ?')
      nred = ndim_vsbsp
      ! update reduced space:
      ! ffvsbsp and ffrsbsp point to ff_trv(iopt)%fhand ...
        if (.not.use_s(iopt)) then
          if (nopt.ne.1) call quit(1,'evpc_core','not this route')
          call optc_update_redsp3
     &       (mred,gred,nred,nroot,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,
     &       iord_vsbsp,ffvsbsp,iord_rsbsp,ffrsbsp,
     &                              me_rhs(iopt)%mel%fhand,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
        else
          call optc_update_redsp4
     &         (mred,sred,gred,nred,nroot,mxsub,
     &       opti_stat%nadd,opti_stat%ndel, .true.,!init,
     &       iord_vsbsp,ffvsbsp,
     &       iord_rsbsp,ffrsbsp,
     &       iord_ssbsp,ffssbsp,me_rhs(iopt)%mel%fhand,
c     &       iord_vsbsp,ffvsbsp(iopt)%fhand,
c     &       iord_rsbsp,ffrsbsp(iopt)%fhand,
c     &       iord_ssbsp,ffssbsp(iopt)%fhand,fdum,
     &       nincore,nwfpar(iopt),lenbuf,xbuf1,xbuf2,xbuf3)
        end if

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
      if (opti_info%shift.eq.0d0.or..not.use_s(1)) then
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
        if (nincore.ge.2) then
          call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,nwfpar)
        else
          call da_sccpvec(me_rhs(iopt)%mel%fhand,iroot,
     &                    ffscr,iroot,
     &                    1d0,nwfpar,xbuf1,lenbuf)
        end if

        idx = (iroot-1)*mxsub + 1
        call optc_expand_vec(vred(idx),ndim_rsbsp,xrsnrm(iroot),.true.,
     &       ffscr,irecscr,1d0,ffrsbsp,iord_rsbsp,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        ! care for shift (incl. metric, if applicable)
        if (opti_info%shift.ne.0d0 .and. .not.use_s(iopt))
     &    call optc_expand_vec(
     &         opti_info%shift*vred(idx:idx-1+ndim_vsbsp),ndim_vsbsp,
     &                   xrsnrm(iroot),.true.,
     &         ffscr,irecscr,1d0,ffvsbsp,iord_rsbsp,
     &         nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        if (opti_info%shift.ne.0d0 .and. use_s(iopt))
     &    call optc_expand_vec(
     &         opti_info%shift*vred(idx:idx-1+ndim_vsbsp),ndim_ssbsp,
     &                   xrsnrm(iroot),.true.,
     &         ffscr,irecscr,1d0,ffssbsp,iord_rsbsp,
     &         nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        ! not yet converged? increase record counter
        if (xrsnrm(iroot).gt.opti_info%thrgrd(iopt)) then
          idxroot(irecscr) = iroot
          irecscr = irecscr+1 
        end if

      end do
      
      ! number of new directions
      nnew = irecscr-1
c dbg
c      print *,'nnew = ',nnew
c dbg

      if (nnew.gt.0) then

        ! reduced space exhausted?
        if (nred+nnew.gt.mxsub) then
          restart_mode = 0
          if (restart_mode.eq.0) then
            ! complete internal restart
            ! assemble orth. subspace exactly spanning the nroot 
            ! currently best solution vectors
            call optc_minspace(
     &           iord_vsbsp,opti_stat%ffvsbsp,
     &           iord_rsbsp,opti_stat%ffrsbsp,
     &           iord_ssbsp,opti_stat%ffssbsp,use_s,
     &           vred,gred,mred,sred,nred,nroot,nroot,mxsub,nopt,
     &           ffscr,nnew,
     &           nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
            ndim_vsbsp = nred
            ndim_rsbsp = nred
            ndim_ssbsp = nred
          else
            call quit(1,'leqc_core','baustelle')
          end if

        end if

        ! divide new directions by preconditioner
        if (nincore.ge.2) then
          call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,nwfpar)
          do iroot = 1, nnew
            call vec_from_da(ffscr,iroot,xbuf1,nwfpar)
            ! scale residual for numerical stability:
c            xnrm = dnrm2(nwfpar,xbuf1,1)
            xnrm = xrsnrm(idxroot(iroot))
            call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,
     &                 opti_info%shift,nwfpar)
            call vec_to_da(ffscr,iroot,xbuf1,nwfpar)
          end do
        else
          do iroot = 1, nnew
c            ! request (nroot-iroot+1)th-last root 
c            irec = ioptc_get_sbsp_rec(-nroot+iroot-1,
c     &           iord_vsbsp,ndim_vsbsp,mxsbsp)
            xnrm = xrsnrm(idxroot(iroot))
            call da_diavec(ffscr,iroot,0d0,
     &                     ffscr,iroot,1d0/xnrm,
     &                      me_dia(1)%mel%fhand,1,0d0,-1d0,
     &                      nwfpar,xbuf1,xbuf2,lenbuf)
          end do
        end if

        ! orthogonalize new directions to existing subspace
        ! and add linear independent ones to subspace
        call optc_orthvec(nadd,.false.,
     &                  opti_stat%ffvsbsp,
     &                      iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                  ffscr,nnew,nopt,
     &                  nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)

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

        do iroot = 1, nroot

          idx = (iroot-1)*mxsub + 1
          call optc_expand_vec(vred(idx),ndim_vsbsp,xdum,.false.,
     &         me_opt(1)%mel%fhand,iroot,0d0,ffvsbsp,iord_vsbsp,
     &         nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        end do

      end if

      return
      end

