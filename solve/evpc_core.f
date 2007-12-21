*----------------------------------------------------------------------*
      subroutine evpc_core(iter,
     &       task,iroute,xrsnrm,xeig,
     &       ffopt,fftrv,ffmvp,ffdia,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       opti_info,opti_stat)
*----------------------------------------------------------------------*
*     core driver for EVP solver
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_file_array.h'
      include 'def_optimize_info.h'
      include 'def_optimize_status.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(inout) ::
     &     task
      integer, intent(inout) ::
     &     iter
      integer, intent(in) ::
     &     iroute, nincore, lenbuf

      type(file_array), intent(in) ::
     &     ffopt(*), fftrv(*), ffmvp(*), ffdia(*)
      type(filinf), intent(in) ::
     &     ffscr

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xrsnrm(opti_info%nroot), xeig(opti_info%nroot,2)

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

* local
      logical ::
     &     zero_vec(opti_stat%ndim_vsbsp)
      integer ::
     &     idx, jdx, kdx, iroot, nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job,
     &     ndim_save, ndel, iopt, lenscr, ifree, restart_mode, ierr
      real(8) ::
     &     cond, xdum, xnrm, xshf
      real(8), pointer ::
     &     gred(:), vred(:), mred(:), eigr(:), eigi(:),
     &     xmat1(:), xmat2(:), xmat3(:), xvec(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, iord_rsbsp(:), iord_vsbsp(:),
     &     nwfpar(:),
     &     ipiv(:), iconv(:), idxroot(:)
      type(filinf), pointer ::
     &     ffrsbsp, ffvsbsp
      type(filinf) ::
     &     fdum

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'evpc_core entered')

      zero_vec(1:opti_stat%ndim_vsbsp) = .false.
      nopt = opti_info%nopt
      nroot = opti_info%nroot
      mxsub = opti_stat%mxdim_sbsp
      mred => opti_stat%sbspmat(1:)
      vred => opti_stat%sbspmat(mxsub**2+1:)
c      vred => opti_stat%sbspmat(2*mxsub**2+1:)
      ndim_rsbsp => opti_stat%ndim_rsbsp
      ndim_vsbsp => opti_stat%ndim_vsbsp
      iord_rsbsp => opti_stat%iord_rsbsp
      iord_vsbsp => opti_stat%iord_vsbsp
      ffrsbsp => opti_stat%ffrsbsp(1)%fhand
      ffvsbsp => opti_stat%ffvsbsp(1)%fhand
      nwfpar => opti_info%nwfpar

      if (nopt.gt.1)
     &     call quit(1,'evpc_core','not yet adapted for nopt>1')

      ! check for previously converged roots
      ifree = mem_alloc_int(iconv,nopt*nroot,'EVP_conv')
      idx = 0
      do iopt = 1, nopt
        do iroot = 1, nroot
          idx = idx+1
          iconv(idx) = 0
          if (iter.gt.1.and.xrsnrm(idx).lt.opti_info%thrgrd(iopt))
     &         iconv(idx) = 1
        end do
      end do

      ifree = mem_alloc_int(idxroot,nopt*nroot,'EVP_idxroot')

      iopt = 1  ! preliminary
      if (ndim_vsbsp.ne.ndim_rsbsp)
     &     call quit(1,'evpc_core','subspace dimensions differ?')
      nred = ndim_vsbsp
      ! update reduced space:
      ! ffvsbsp and ffrsbsp point to ff_trv(iopt)%fhand ...
      call optc_update_redsp3
     &       (mred,xdum,nred,0,mxsub,
     &       opti_stat%nadd,opti_stat%ndel,
     &       iord_vsbsp,ffvsbsp,iord_rsbsp,ffrsbsp,fdum,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)

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

      call eigen_asym(nred,xmat1,eigr,eigi,xmat2,xmat3,ierr)

      ! copy back to vred with mxsub as leading dim
      kdx = 0
      do idx = 1, nred
        do jdx = 1, nred
          kdx = kdx+1
          vred((idx-1)*mxsub+jdx) = xmat2(kdx) 
        end do
      end do

      if (ntest.ge.50) then
        write(luout,*) 'Eigenvalues in subspace:'
        do idx = 1, nred
          write(luout,'(2x,i4,x,g20.12,x,g20.12)')
     &         idx, eigr(idx), eigi(idx)
        end do
        if (ntest.ge.100) then
          write(luout,*) 'eigenvectors:'
          call wrtmat2(vred,nred,nred,mxsub,mxsub)
        end if
      end if

      xeig(1:nroot,1) = eigr(1:nroot)
      xeig(1:nroot,2) = eigi(1:nroot)

      irecscr = 1
      do iroot = 1, nroot
        ! assemble residual in full space
        
        idx = (iroot-1)*mxsub + 1
        xvec(1:nred) = vred(idx:idx+nred-1)
        call optc_expand_vec(xvec,ndim_rsbsp,xrsnrm(iroot),.false.,
     &       ffscr,irecscr,0d0,ffrsbsp,iord_rsbsp,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        xvec(1:nred) = -eigr(iroot)*xvec(1:nred)
        call optc_expand_vec(xvec,ndim_vsbsp,xrsnrm(iroot),.true.,
     &       ffscr,irecscr,1d0,ffvsbsp,iord_vsbsp,
     &       nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        ! not yet converged? increase record counter
        if (xrsnrm(iroot).gt.opti_info%thrgrd(iopt)) then
          idxroot(irecscr) = iroot
          irecscr = irecscr+1 
        end if

      end do
      
      ! number of new directions
      nnew = irecscr-1
c dbg
      print *,'nnew = ',nnew
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
     &           iord_vsbsp,ffvsbsp,iord_rsbsp,ffrsbsp,
     &           vred,xdum,mred,nred,nroot,0,mxsub,
     &           ffscr,nnew,
     &           nincore,nwfpar,lenbuf,xbuf1,xbuf2,xbuf3)
            ndim_vsbsp = nred
            ndim_rsbsp = nred
          else
            call quit(1,'evpc_core','baustelle')
          end if

        end if

        ! divide new directions by preconditioner
        if (nincore.ge.2) then
          call vec_from_da(ffdia(iopt)%fhand,1,xbuf2,nwfpar)
          do iroot = 1, nnew
            call vec_from_da(ffscr,iroot,xbuf1,nwfpar)
            ! scale residual for numerical stability:
c            xnrm = dnrm2(nwfpar,xbuf1,1)
            xnrm = xrsnrm(idxroot(iroot))
            xshf = -xeig(idxroot(iroot),1)
            call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,xshf,nwfpar)
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
     &                      ffdia,1,0d0,-1d0,
     &                      nwfpar,xbuf1,xbuf2,lenbuf)
          end do
        end if

        ! orthogonalize new directions to existing subspace
        ! and add linear independent ones to subspace
        call optc_orthvec(nadd,
     &                  ffvsbsp,iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                  ffscr,nnew,
     &                  nwfpar,nincore,xbuf1,xbuf2,xbuf3,lenbuf)

        ! set nadd
        if (nadd.eq.0)
     &       call quit(0,'evpc_core',
     &       'solver in problems: only linear dependent '//
     &       'new directions?')
        opti_stat%nadd = nadd

        ! |Mv> subspace organisation should be identical to |v> subsp.
        ndim_rsbsp = ndim_vsbsp
        iord_rsbsp = iord_vsbsp

      else
        ! if all converged: assemble vectors 

        do iroot = 1, nroot

          idx = (iroot-1)*mxsub + 1
          call optc_expand_vec(vred(idx),ndim_vsbsp,xdum,.false.,
     &         ffopt(1)%fhand,iroot,1d0,ffvsbsp,iord_vsbsp,
     &         nincore,nwfpar,lenbuf,xbuf1,xbuf2)

        end do

      end if

      return
      end

