*----------------------------------------------------------------------*
      subroutine leqc_init(xrsnrm,iroute,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,me_met,me_scr,
     &       nincore,lenbuf,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,use_s,
     &       opti_info,opti_stat,
     &       orb_info,op_info,str_info,strmap_info)
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

      real(8), intent(inout) ::
     &     xrsnrm(*)
      integer, intent(in) ::
     &     iroute, nincore, lenbuf

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_dia(*), 
     &     me_trv(*), me_mvp(*), me_rhs(*), me_scr(*)
      type(me_list_array), intent(inout) ::
     &     me_met(*)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffrhs(*), ffdia(*)

      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend

      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      logical, intent(in) ::
     &     use_s(*)

* local
      logical ::
     &     zero_vec(opti_stat%ndim_vsbsp)
      integer ::
     &     idx, jdx, kdx, iroot, irhs,  nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job,
     &     ndim_save, ndel, iopt, lenscr, ifree, restart_mode, nselect
      real(8) ::
     &     cond, xdum, xnrm
      real(8), pointer ::
     &     gred(:), vred(:), mred(:),
     &     xmat1(:), xmat2(:), xvec(:), xret(:)
      integer, pointer ::
     &     ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &     iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:),
     &     nwfpar(:),
     &     ipiv(:), iconv(:), idxroot(:), idxselect(:)
      type(filinf), pointer ::
     &     ffvsbsp, ffscr, ffmet
      type(filinf), target ::
     &     fdum

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
      ffscr => opti_stat%ffscr(1)%fhand
      ffvsbsp => opti_stat%ffvsbsp(1)%fhand
      nwfpar => opti_info%nwfpar

      if (nopt.gt.1)
     &     call quit(1,'leqc_init','not yet adapted for nopt>1')
      if (nroot.gt.mxsub)
     &     call quit(1,'leqc_init','?? nroot>mxsub ??')

      if (nincore.ge.2) then
        do iopt = 1, nopt
          ! read diagonal pre-conditioner
          if (ntest.ge.100) then
            write(luout,*) 'current ME-list: ',
     &            trim(me_dia(iopt)%mel%label)
          end if
          call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,nwfpar(iopt))
          if (ntest.ge.100)
     &        write(luout,*) 'xbuf2 norm = ',
     &                       dnrm2(nwfpar(iopt),xbuf2,1)
          do iroot = 1, nroot
            ! divide rhs's by preconditioner
            call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,
     &                       nwfpar(iopt))
            if (ntest.ge.100)
     &           write(luout,*) 'xbuf1 norm = ',
     &                          dnrm2(nwfpar(iopt),xbuf1,1) 
            xnrm = dnrm2(nwfpar(iopt),xbuf1,1)
            xrsnrm(iroot) = xnrm
            ! %shift is for shifted LEQ
            call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,opti_info%shift,
     &                 nwfpar(iopt))
            if (ntest.ge.100)
     &           write(luout,*) 'xbuf1 after division: ' //
     &                          ' norm = ', dnrm2(nwfpar(iopt),xbuf1,1)
            call vec_to_da(ffscr,iroot,xbuf1,nwfpar(iopt))
          end do
        end do
      else

        call quit(1,'leqc_init','incore<2: do the programming')

        ! ... something using:
            call da_diavec(ffscr,iroot,0d0,
     &                     ffscr,iroot,1d0/xnrm,
     &                      me_dia(iopt)%mel%fhand,1,0d0,-1d0,
     &                      nwfpar(iopt),xbuf1,xbuf2,lenbuf)

      end if

      if (ndim_vsbsp.ne.0)
     &     call quit(1,'leqc_init','ndim_vsbsp.ne.0 ???')

      do iopt = 1, nopt
        if (use_s(iopt)) then
          ! assign op. with list containing the scratch trial vector
          call assign_me_list(me_scr(iopt)%mel%label,
     &                        me_opt(iopt)%mel%op%name,op_info)

          ! calculate metric * scratch trial vector
          allocate(xret(depend%ntargets),idxselect(depend%ntargets))
          nselect = 0
          call select_formula_target(idxselect,nselect,
     &                me_met(iopt)%mel%label,depend,op_info)
          do iroot = 1, nroot
            call switch_mel_record(me_met(iopt)%mel,iroot)
            call switch_mel_record(me_scr(iopt)%mel,iroot)
            call frm_sched(xret,flist,depend,idxselect,nselect,
     &                  op_info,str_info,strmap_info,orb_info)
            me_met(iopt)%mel%fhand%last_mod(iroot) = -1
          end do
          deallocate(xret,idxselect)

          ! reassign op. with list containing trial vector
          call assign_me_list(me_trv(iopt)%mel%label,
     &                        me_opt(iopt)%mel%op%name,op_info)
          ffmet => me_met(1)%mel%fhand
        else
          ffmet => fdum
        end if
      end do

      ! orthogonalize initial subspace
      call optc_orthvec(nadd,.false.,
     &                  opti_stat%ffvsbsp,
     &                  iord_vsbsp,ndim_vsbsp,mxsub,zero_vec,
     &                  use_s,0,ffmet,ffscr,nroot,nopt,
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

