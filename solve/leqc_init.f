*----------------------------------------------------------------------*
      subroutine leqc_init(xrsnrm,iroute,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,
     &       nincore,lenbuf,ffscr,
     &       xbuf1,xbuf2,xbuf3,
     &       opti_info,opti_stat)
*----------------------------------------------------------------------*
*     initialization of LEQ solver
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

      real(8), intent(inout) ::
     &     xrsnrm(*)
      integer, intent(in) ::
     &     iroute, nincore, lenbuf

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_dia(*), 
     &     me_trv(*), me_mvp(*), me_rhs(*)
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
     &     gred(:), vred(:), mred(:),
     &     xmat1(:), xmat2(:), xvec(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &     iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:),
     &     nwfpar(:),
     &     ipiv(:), iconv(:), idxroot(:)
      type(filinf), pointer ::
     &     ffvsbsp

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'leqc_init entered')

      nopt = opti_info%nopt
      nroot = opti_info%nroot
      mxsub = opti_stat%mxdim_sbsp
      mred => opti_stat%sbspmat(1:)
      gred => opti_stat%sbspmat(mxsub**2+1:)
      vred => opti_stat%sbspmat(2*mxsub**2+1:)
      ndim_rsbsp => opti_stat%ndim_rsbsp
      ndim_vsbsp => opti_stat%ndim_vsbsp
      ndim_ssbsp => opti_stat%ndim_ssbsp
      iord_rsbsp => opti_stat%iord_rsbsp
      iord_vsbsp => opti_stat%iord_vsbsp
      iord_ssbsp => opti_stat%iord_ssbsp
      ffvsbsp => opti_stat%ffvsbsp(1)%fhand
      nwfpar => opti_info%nwfpar

      if (nopt.gt.1)
     &     call quit(1,'leqc_init','not yet adapted for nopt>1')
      if (nroot.gt.mxsub)
     &     call quit(1,'leqc_init','?? nroot>mxsub ??')

      if (nincore.ge.2) then
        do iopt = 1, nopt
          ! read diagonal pre-conditioner
          call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,nwfpar)
          do iroot = 1, nroot
            ! divide rhs's by preconditioner
            call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,nwfpar)
            xnrm = dnrm2(nwfpar,xbuf1,1)
            xrsnrm(iroot) = xnrm
            call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,0d0,nwfpar)
            call vec_to_da(ffscr,iroot,xbuf1,nwfpar)
          end do
        end do
      else

        call quit(1,'leqc_init','incore<2: do the programming')

        ! ... something using:
            call da_diavec(ffscr,iroot,0d0,
     &                     ffscr,iroot,1d0/xnrm,
     &                      me_dia(iopt)%mel%fhand,1,0d0,-1d0,
     &                      nwfpar,xbuf1,xbuf2,lenbuf)

      end if

      if (ndim_vsbsp.ne.0)
     &     call quit(1,'leqc_init','ndim_vsbsp.ne.0 ???')
      ! orthogonalize initial subspace
      call optc_orthvec(nadd,
     &                  opti_stat%ffvsbsp,
     &                     iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                  ffscr,nroot,nopt,
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

      return
      end

