*----------------------------------------------------------------------*
      subroutine leqc_init(xrsnrm,iroute,
     &       me_opt,me_trv,me_mvp,me_rhs,me_dia,me_met,me_scr,
     &       me_special,nspecial,
     &       nincore,lenbuf,
     &       xbuf1,xbuf2,xbuf3,
     &       flist,depend,use_s,
     &       fspc,nspcfrm,
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

      integer, intent(in) ::
     &     iroute, nincore, lenbuf, nspecial, nspcfrm

      type(me_list_array), intent(inout) ::
     &     me_opt(*), me_dia(*), 
     &     me_trv(*), me_mvp(*), me_rhs(*), me_scr(*),
     &     me_special(nspecial)
      type(me_list_array), intent(inout) ::
     &     me_met(*)
c      type(file_array), intent(in) ::
c     &     ffopt(*), fftrv(*), ffmvp(*), ffrhs(*), ffdia(*)

      type(formula_item), intent(inout) ::
     &     flist, fspc(nspcfrm)
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
     &     xrsnrm(opti_info%nroot,opti_info%nopt)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      logical, intent(in) ::
     &     use_s(*)

* local
      logical ::
     &     zero_vec(opti_stat%ndim_vsbsp), trafo
      integer ::
     &     idx, jdx, kdx, iroot, irhs,  nred, nadd, nnew, irecscr,
     &     imet, idamp, nopt, nroot, mxsub, lenmat, job, jopt,
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
      type(file_array), pointer ::
     &     ffvsbsp(:), ffscr(:), ffmet(:)
      type(filinf), target ::
     &     fdum

      integer, external ::
     &     ioptc_get_sbsp_rec
      real(8), external ::
     &     dnrm2

      trafo = .false.

      nopt = opti_info%nopt
      nroot = opti_info%nroot
      mxsub = opti_stat%mxdim_sbsp
c      mred => opti_stat%sbspmat(1:)
c      gred => opti_stat%sbspmat(mxsub**2+1:)
c      vred => opti_stat%sbspmat(2*mxsub**2+1:)
      ndim_rsbsp => opti_stat%ndim_rsbsp
      ndim_vsbsp => opti_stat%ndim_vsbsp
      ndim_ssbsp => opti_stat%ndim_ssbsp
      iord_rsbsp => opti_stat%iord_rsbsp
      iord_vsbsp => opti_stat%iord_vsbsp
      iord_ssbsp => opti_stat%iord_ssbsp
      ffscr => opti_stat%ffscr
      ffvsbsp => opti_stat%ffvsbsp
      nwfpar => opti_info%nwfpar
      zero_vec = .false.

      allocate(ffmet(nopt))

c      if (nopt.gt.1)
c     &     call quit(1,'leqc_init','not yet adapted for nopt>1')
      if (nroot.gt.mxsub)
     &     call quit(1,'leqc_init','?? nroot>mxsub ??')

      if (nincore.ge.2) then

        do iopt = 1, nopt

         !select case(opti_info%typ_prc(iopt))
         if (opti_info%typ_prc(iopt).eq.optinf_prc_file.or.
     &       opti_info%typ_prc(iopt).eq.optinf_prc_mixed.or.
     &       opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
         !case(optinf_prc_file,optinf_prc_mixed,optinf_prc_traf)

          trafo = opti_info%typ_prc(iopt).eq.optinf_prc_traf            

          ! read diagonal pre-conditioner
          if (ntest.ge.100) then
            write(lulog,*) 'current ME-list: ',
     &            trim(me_dia(iopt)%mel%label)
          end if
          call vec_from_da(me_dia(iopt)%mel%fhand,1,xbuf2,nwfpar(iopt))
          if (ntest.ge.100)
     &        write(lulog,*) 'xbuf2 norm = ',
     &                       dnrm2(nwfpar(iopt),xbuf2,1)
          do iroot = 1, nroot
            if (trafo) then
              ! the present version works only, if we do this for root
              ! number 1; this is, because we modify the norm below
              ! => we have to switch around the loops
              if (iopt.ne.1)
     &             call quit(1,'leqc_init','route with trafo: problem')
              call optc_traf(me_special(2)%mel,1,xrsnrm(iroot,iopt),
     &                    me_rhs(iopt)%mel,iroot,
     &                    fspc(1),'B',me_special,nspecial,
     &                    nwfpar(iopt),xbuf1,
     &                    orb_info,op_info,str_info,strmap_info)
              call vec_from_da(me_special(2)%mel%fhand,1,xbuf1,
     &                       nwfpar(iopt))
            else
              call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,
     &                       nwfpar(iopt))
            end if
            ! divide rhs's by preconditioner
            if (ntest.ge.100)
     &           write(lulog,*) 'xbuf1 norm = ',
     &                          dnrm2(nwfpar(iopt),xbuf1,1) 
c            xnrm = dnrm2(nwfpar(iopt),xbuf1,1)
c            xrsnrm(iroot,iopt) = xnrm
            xnrm = 0d0
            do jopt = 1, nopt
              xnrm = xnrm+xrsnrm(iroot,jopt)**2
            end do
            xnrm = sqrt(xnrm)
            ! %shift is for shifted LEQ
            call diavc(xbuf1,xbuf1,1d0/xnrm,xbuf2,opti_info%shift,
     &                 nwfpar(iopt))
            if (ntest.ge.100)
     &           write(lulog,*) 'xbuf1 after division: ' //
     &                          ' norm = ', dnrm2(nwfpar(iopt),xbuf1,1)
            if (trafo) then
              call vec_to_da(me_special(2)%mel%fhand,1,
     &                       xbuf1,nwfpar(iopt))
              call optc_traf(me_opt(iopt)%mel,iroot,xdum,
     &                    me_special(2)%mel,1,
     &                    fspc(1),'F',me_special,nspecial,
     &                    nwfpar(iopt),xbuf1,
     &                    orb_info,op_info,str_info,strmap_info)
              ! copy to scr list
              ! original list was used to ensure spin symmetry if needed
              call switch_mel_record(me_scr(iopt)%mel,iroot)
              call list_copy(me_opt(iopt)%mel,me_scr(iopt)%mel,.false.)
            else
              call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,nwfpar(iopt))
            end if
          end do ! iroot
         !case(optinf_prc_blocked)
         else if (opti_info%typ_prc(iopt).eq.optinf_prc_blocked) then
          if (nincore.lt.3)
     &         call quit(1,'leqc_init',
     &         'I need at least 3 incore vectors (prc_blocked)')
          do iroot = 1, nroot
            call vec_from_da(me_rhs(iopt)%mel%fhand,iroot,xbuf1,
     &                       nwfpar(iopt))
c            xnrm = dnrm2(nwfpar(iopt),xbuf1,1)
c            xrsnrm(iroot,iopt) = xnrm
            xnrm = 0d0
            do jopt = 1, nopt
              xnrm = xnrm+xrsnrm(iroot,jopt)**2
            end do
            xnrm = sqrt(xnrm)
            call dscal(nwfpar(iopt),1d0/xnrm,xbuf1,1)
            call optc_prc_special2(me_rhs(iopt)%mel,me_special,
     &                                                      nspecial,
     &                         me_opt(iopt)%mel%op%name,
     &                         opti_info%shift,
     &                         nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                         orb_info,op_info,str_info,strmap_info)
            call vec_to_da(ffscr(iopt)%fhand,iroot,xbuf1,nwfpar(iopt))
          end do
         !case default
         else
           call quit(1,'leqc_init','unknown preconditioner type')
         end if
         !end select
        end do
      else

        call quit(1,'leqc_init','incore<2: do the programming')

        ! ... something using:
            call da_diavec(ffscr(iopt)%fhand,iroot,1,0d0,
     &                     ffscr(iopt)%fhand,iroot,1,1d0/xnrm,
     &                      me_dia(iopt)%mel%fhand,1,1,0d0,-1d0,
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
     &             .true.,.false.,op_info,str_info,strmap_info,orb_info)
            me_met(iopt)%mel%fhand%last_mod(iroot) = -1
          end do
          deallocate(xret,idxselect)

          ! reassign op. with list containing trial vector
          call assign_me_list(me_trv(iopt)%mel%label,
     &                        me_opt(iopt)%mel%op%name,op_info)
          ffmet(iopt)%fhand => me_met(iopt)%mel%fhand
        else
          ffmet(iopt)%fhand => fdum
        end if
      end do

      ! orthogonalize initial subspace
      call optc_orthvec(nadd,.false.,
     &                  opti_stat%ffssbsp,iord_ssbsp,1d0, !1d0:dummy
     &                  ffvsbsp,
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

      deallocate(ffmet)

      return
      end

