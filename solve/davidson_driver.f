
!>    driver for the davidson algoritmus
!!
!!    Davidson
!!    on entry: Mv-product on me_mvp me

      subroutine davidson_driver(
     &     dvdsbsp,
     &     iter,  task, nnew,
     &     nroot, nopt,
     &     trafo, use_s,
     &     xrsnrm , xeig, reig,
     &     me_opt,me_dia,
     &     me_met,me_metort,
     &     me_scr,me_res,
     &     me_trv,me_mvp,me_vort,me_mvort,
     &     me_special, nspecial,
     &     xbuf1,xbuf2, xbuf3, nincore,lenbuf,
     &     flist,depend,
     &     fspc,nspcfrm,
     &     opti_info, opti_stat,
     &     orb_info, op_info, str_info,strmap_info
     &     )

      implicit none

      include 'stdunit.h'
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
      include 'def_davidson_subspace.h'

      integer, parameter ::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="davidson_driver"
      real(8),parameter::
     &     thrgrd_e=0.1d0

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer, intent(inout) ::
     &     iter, task, nnew
      logical, intent(in)::
     &     trafo(*), use_s(*)
      integer, intent(in) ::
     &     nincore, lenbuf, nspecial, nspcfrm,
     &     nroot, nopt

      type(me_list_array), intent(inout) :: !inout to be sure
     &     me_opt(*), me_dia(*),
     &     me_met(*), me_metort(*),
     &     me_scr(*),me_res(*),
     &     me_trv(*), me_mvp(*),me_vort(*),me_mvort(*),
     &     me_special(*)


      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat

      real(8), intent(inout) ::
     &     xrsnrm(nroot,nopt),
     &     xeig(nroot),reig(nroot)      ! xeig persistent copy of eigenvalues to detect root switching
                                        ! reig copy for outer routine to commutincate eigevalues (only nnew entries are defined)
      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)
      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ,intent(in)::
     &     depend
      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)


      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in)::
     &     strmap_info


*----------------------------------------------------------------------*
!     locals

      integer::
     &     iopt, jopt,
     &     iroot,
     &     maxvec,
     &     irecres,        !record pointer for contracting lists to be updated on me_res, number of new_lists
     &     maxiter              !aliases for opti_info fields

      integer,pointer::
     &     typ_prc(:),nwfpar(:)
      logical::
     &     conv, zero_nrm
      real(8)::
     &     xnrm                !temporary variable for some norms
! lists for the mv-product and vector in the orthogonal space
      integer,external::
     &     dvdsbsp_get_nnew_vvec
      real(8),external::
     &     xnormop
      real(8)::
     &     lrsnrm(nroot,nopt), ! local copy of xrsnrm
     &     leig(nroot) !local variable for eigenvalues





      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.ge.100) write (lulog,*) "iteration ",iter

      maxiter = opti_info%maxmacit
      typ_prc => opti_info%typ_prc
      nwfpar => opti_info%nwfpar
cc     the lists are alialized as if the following block had been executed ARNE 8.9.16
cc
c      do iopt = 1,nopt
c         if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
c            me_vort(iopt)%mel => me_special(1)%mel
c            me_mvort(iopt)%mel => me_scr(iopt)%mel
c         elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
c            if (nspecial .lt. 4) call quit(1,'evpc_core',
c     &           'TR0 with <4 me_lists')
c            me_vort(iopt)%mel => me_special(4)%mel
c            me_mvort(iopt)%mel => me_scr(iopt)%mel
c         else !no transformation needed
c            me_vort(iopt)%mel => me_trv(iopt)%mel
c            me_mvort(iopt)%mel => me_mvp(iopt)%mel
c         end if
c      end do

      irecres=0
      do iroot=1,nnew
         
         if (ntest.ge.100)then
            do iopt=1,nopt
               write(lulog,*) "root no.",iroot
              call print_list("entry trv",me_trv(iopt)%mel,
     &              "LIST",0d0,0d0,
     &              orb_info,str_info)
            end do
         end if
         do iopt=1,nopt
            call switch_mel_record(me_mvp(iopt)%mel,iroot)
            call switch_mel_record(me_mvort(iopt)%mel,iroot)
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            call switch_mel_record(me_vort(iopt)%mel,iroot)
            if (use_s(iopt))then
               call switch_mel_record(me_met(iopt)%mel,iroot)
               call switch_mel_record(me_metort(iopt)%mel,iroot)
            end if
            if (trafo(iopt)) then
               call transform_forward_wrap(flist,depend,
     &              me_special,me_mvp,me_mvort, !mvp-> mvort
     &              xrsnrm,
     &              nroot, iroot, iopt, iroot, nspecial,
     &              me_trv,
     &              op_info, str_info, strmap_info, orb_info, opti_info)
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                 call set_blks(me_Mvort(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
               endif
              
              if (use_s(iopt))then
                 call transform_forward_wrap(flist,depend,
     &                me_special,me_met,me_metort, !met-> metort
     &                xrsnrm, nroot,
     &                iroot, iopt, iroot, nspecial,
     &                me_trv,
     &                op_info,str_info,strmap_info, orb_info, opti_info)
                 if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                    call set_blks(me_metort(iopt)%mel,
     &                   "P,H|P,V|V,H|V,V",0d0)
                 end if
              else
                 me_metort(iopt)%mel=> null()
              end if
           else
              if (use_s(iopt))then
                 me_metort(iopt)%mel=> me_met(iopt)%mel
              else
                 me_metort(iopt)%mel=>null()
              end if
           end if
        end do                  !iopt
      end do
      


      call dvdsbsp_append_vecs(dvdsbsp,
     &     me_metort,
     &     me_mvort,
     &     me_vort,
     &     nnew,
     &     nopt,                !mvort & metort into dvdsbsp and update the vMv-matrix and vSv-matrix
     &     xbuf1, xbuf2, xbuf3, nincore, lenbuf)
      
      if ( nnew.eq.0)
     &     call quit(0,i_am,
     &     "only linear depended new directions generated")

      if (nroot.gt. dvdsbsp_get_nfree(dvdsbsp) )then
         call dvdsbsp_compress(dvdsbsp, nroot,
     &        nopt, me_scr,     !scr as scratch
     &        xbuf1, xbuf2, lenbuf ,nincore)
      end if



      call davidson_assemble_residuals(dvdsbsp,
     &     leig, nnew,
     &     me_scr, nopt, nroot, xrsnrm, !temporary assemble on scr
     &     xbuf1, xbuf2, nincore, lenbuf)


!................. check convergence
      irecres=0
      do iroot=1,nnew
         ! test energy criteria
         if (iter .gt.1      ! if iter was 1 old energy was 0 convergence check is stupid also
     &        .and. check_e_convergence(iter,leig(iroot) ,xeig(iroot),
     &        thrgrd_e) )then
            conv = .true. ! yet
         else if(iter.gt.1)then
            conv=.false.
            call warn_e_convergence(lulog,iroot,leig(iroot),xeig(iroot))
         else
            conv=.true.
         end if
         ! test gradient(residual) criteria
         conv= conv .and.
     &        check_r_convergence(xrsnrm,opti_info%thrgrd,iroot, nopt)

         if (.not.conv)then
            irecres=irecres+1
            do iopt=1,nopt
               call switch_mel_record(me_scr(iopt)%mel,iroot)
               call switch_mel_record(me_res(iopt)%mel,irecres) ! all updated vectors are one after the other
               call list_copy(me_scr(iopt)%mel,me_res(iopt)%mel, !scr ->res
     &              .false.)    !collecting vectors which will lead to new direction on me_res (me_residual)
            end do
            reig(irecres)=leig(iroot)
         end if
      end do
      xeig=leig

      nnew=irecres              ! 0 if all converged
      

!................. if all converged, assemble results.
      if (nnew .eq. 0
     &     .or. iter .ge. maxiter) then
         do iopt=1,nopt
            if(.not.trafo(iopt))me_vort(iopt)%mel => me_opt(iopt)%mel
         end do
         maxvec=min(nroot,mel_get_maxrec(me_opt(1)%mel))
         call davidson_assemble_results(dvdsbsp,
     &        xeig,
     &        me_vort, nopt, maxvec, lrsnrm, ! results on me_vort
     &        xbuf1, xbuf2, nincore, lenbuf)
         
         do iroot=1,maxvec
            do iopt=1,nopt
               call switch_mel_record(me_vort(iopt)%mel,iroot)
               call switch_mel_record(me_opt(iopt)%mel,iroot)
            end do
            if (ntest.ge.00)then
               write(lulog,*) "root:",iroot
               do iopt=1,nopt
                  call print_list("result in orth space",
     &                 me_vort(iopt)%mel,
     &                 "LIST",0d0,0d0,
     &                 orb_info,str_info)
               end do
            end if
           
!            call vec_normalize(me_vort, nopt, xbuf1, xbuf2, lenbuf)
            do iopt=1,nopt
               if (trafo(iopt) ) then
                  call transform_back_wrap(flist,depend,
     &                 me_special,me_vort,me_opt, !vort -> opt !
     &                 iroot, iopt, nspecial,
     &                 me_trv,
     &                 op_info, str_info, strmap_info,
     &                 orb_info, opti_info)


!     else me_vort => me_trv => me_opt
               end if
            end do
            if (ntest.ge.00)then
               write(lulog,*) "root:",iroot
               do iopt=1,nopt
                  call print_list("result in normal space",
     &                 me_opt(iopt)%mel,
     &                 "LIST",0d0,0d0,
     &                 orb_info,str_info)
               end do
            end if
         end do
         if(nnew.eq.0)then
            write(lulog,'(x,a,i5,a)')
     &           'CONVERGED IN ',iter,' ITERATIONS'
            if (luout.ne.lulog.and.iprlvl.ge.5)
     &           write(luout,'(x,a,i5,a)')
     &           'CONVERGED IN ',iter,' ITERATIONS'
         else
            write(lulog,'(x,a,i5,a)') "Stopping after",iter,"iterations"
            call warn('linear solver', 'NO CONVERGENCE OBTAINED')
         end if
         task=8
         return                 !END of method !!!!!!!!!!!!!!!!!!!
      else
         task=4       ! calculate new mv product when returning
      end if

      if (ntest.ge.100)then
         do iopt=1,nopt
            call print_list("preconditioner",me_dia(iopt)%mel,
     &           "LIST",0d0,0d0,
     &           orb_info,str_info)
         end do
      end if
      do iroot = 1, nnew
         do iopt=1,nopt
            call switch_mel_record(me_res(iopt)%mel,iroot)
            call switch_mel_record(me_vort(iopt)%mel,iroot)
         end do

         if (ntest.ge.100)then
            do iopt=1,nopt
               write(lulog,*) "root no.",iroot
               call print_list("residual",me_res(iopt)%mel,
     &              "LIST",0d0,0d0,
     &              orb_info,str_info)
            end do
         end if
         do iopt = 1,nopt
           
            xnrm=0.0
            do jopt = 1,nopt
               xnrm = xnrm+xrsnrm(iroot,jopt)**2
            end do

            xnrm=sqrt(xnrm)
            call apply_preconditioner(
     &           me_res(iopt)%mel, me_dia(iopt)%mel, me_vort(iopt)%mel,  !res -> vort
     &           nwfpar(iopt),
     &           opti_info%typ_prc(iopt),
     &           xeig(iroot), xnrm, nnew, iroot,
     &           me_opt(iopt)%mel%op, me_trv(iopt)%mel,me_scr(iopt)%mel,
     &           me_special, nspecial,
     &           fspc, nspcfrm,
     $           xbuf1, xbuf2, xbuf3, lenbuf, nincore,
     &           nopt.eq.1, op_info, str_info, strmap_info)
         end do
         if (ntest.ge.100)then
            do iopt=1,nopt
               write(lulog,*) "root no.",iroot
               call print_list("unprojected vort",me_vort(iopt)%mel,
     &              "LIST",0d0,0d0,
     &              orb_info,str_info)
            end do
         end if
         call vecsp_orthvec(dvdsbsp%vspace, me_vort, nopt,
     &        xbuf1, xbuf2, lenbuf)
         if (ntest.ge.100)then
            do iopt=1,nopt
               write(lulog,*) "root no.",iroot
               call print_list("projected trv",me_vort(iopt)%mel,
     &              "LIST",0d0,0d0,
     &              orb_info,str_info)
            end do
         end if

      end do                    !iroot

      do iopt=1,nopt
         do iroot=1,nnew
            call switch_mel_record(me_vort(iopt)%mel,iroot)
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            if (trafo(iopt) ) then
               call transform_back_wrap(flist,depend,
     &              me_special,me_vort,me_trv, !vort -> trv !new_trialvector created
     &              iroot, iopt,nspecial,
     &              me_trv,
     &              op_info, str_info, strmap_info,
     &              orb_info, opti_info)
!            else me_vort => me_trv
            end if
         end do
      end do                    !iopt
      
      if (ntest.ge.100)then
         do iopt=1,nopt
            write(lulog,*) "root no.",iroot
            call print_list("retransformed trv",me_trv(iopt)%mel,
     &           "LIST",0d0,0d0,
     &           orb_info,str_info)
         end do
      end if
      return
      contains



!#######################################################################
!   apply the preconditioner
!######################################################################
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine apply_preconditioner(
     &     me_resid, me_diag, me_result, lenme,
     &     typ_prc,
     &     xeig, xnrm, nnew, iroot,
     &     op_opt, me_opt, me_scr, ! opt needed for some assignments me_scr as scratch list, may or may not be used
     &     me_special, nspecial,!
     &     fspc, nfspc,         !head of special formulars
     $     xbuf1, xbuf2, xbuf3, lbuf, nincore,
     &     renormalize, op_info, str_info, strmap_info)
*----------------------------------------------------------------------*
      implicit none
      character(len=*),parameter::
     &     i_am="apply_preconditioner"
      integer,parameter::
     &     ntest=00


      type(me_list), intent(in), target::
     &     me_resid,
     &     me_diag,
     &     me_result,
     &     me_scr               !me_list that might be used


      integer,intent(in)::
     &     lenme,
     &     nnew,
     &     nspecial,
     &     nfspc,
     &     iroot,
     &     lbuf, nincore

      integer,intent(in)::
     &     typ_prc

      real(8),intent(in)::
     &     xeig
      real(8),intent(inout)::
     &     xbuf1(*), xbuf2(*), xbuf3(*)

      real(8),intent(inout)::
     &     xnrm                 !> norm of the associated vectors
      logical, intent(in)::
     &     renormalize
      type(formula_item), intent(in)::
     &     fspc(*)
      type(me_list_array)::
     &     me_special(*)

      type(operator),intent(inout)::
     &     op_opt

      type(me_list),intent(inout)::
     &     me_opt

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in)::
     &     strmap_info
      type(me_list), pointer::
     &     me_intm
      logical ::
     &     spin_prj
      real(8)::
     &     xshf
      real(8),external::
     &     dnrm2
      integer::
     &     idxdbg


      select case(typ_prc)
      case(optinf_prc_file,optinf_prc_traf,optinf_prc_traf_spc
     &     ,optinf_prc_spinp,optinf_prc_prj,optinf_prc_spinrefp)
      if (typ_prc.eq.optinf_prc_spinp.or.
     &     typ_prc.eq.optinf_prc_prj.or.
     &     typ_prc.eq.optinf_prc_spinrefp)then
         spin_prj=.true.
         me_intm => me_result
      else
         spin_prj=.false.
         me_intm => me_result
      end if

      if (nincore.ge.2) then
         call vec_from_da(
     &        me_diag%fhand,1, xbuf2,lenme)
         call vec_from_da(me_resid%fhand,iroot,xbuf1,
     &        lenme)


c     dbg
c     print *,"xnorm:",xnrm
c     dbgend
c     xnrm = 1d0
c     dbg
c     print *,"precon not yet applied. "
c     do idxdbg = 1, nwfpar(iopt)
c     print *,idxdbg,xbuf1(idxdbg)
c     end do
c     dbgend
! scale residual for numerical stability:
         xshf = -xeig
         call diavc(xbuf1,xbuf1,1d0/xnrm,
     &        xbuf2,xshf,lenme)
c     dbg
c     print *,"debug shift:",xshf
c     print *,"precon applied. ",nwfpar(iopt),"entries"
c     do idxdbg = 1, nwfpar(iopt)
c     print *,idxdbg,xbuf1(idxdbg),xbuf2(idxdbg)
c     end do
c     dbgend
         if (renormalize) then
            xnrm = dnrm2(lenme,xbuf1,1)
            call dscal(lenme,1d0/xnrm,xbuf1,1)
         end if
         call vec_to_da(me_result%fhand,iroot,xbuf1,
     &        lenme)
      else
c     ! request (nroot-iroot+1)th-last root
c     irec = ioptc_get_sbsp_rec(-nroot+iroot-1,
c     &         iord_vsbsp,ndim_vsbsp,mxsbsp)
         xshf = -xeig
! decrease xshf in first iteration
         if (iter .eq. 1) xshf=0.8d0*xshf
         call da_diavec(me_resid%fhand,iroot,0d0,
     &        me_intm%fhand,iroot,
     &        1d0/xnrm,me_diag%fhand,
     &        1,xshf,-1d0,
     &        lenme,xbuf1,xbuf2,lenbuf)
      end if

      if (spin_prj)then !intm != result
! assign op. with list containing the scratch trial vector
c     dbg
c     print *, "assign ",me_scr(iopt)%mel%label," to ",
c     &              me_opt(iopt)%mel%op%name
c     dbgend
         call assign_me_list(me_intm%label,
     &        op_opt%name,op_info)
         call switch_mel_record(me_intm,iroot)
         call reset_file_rec(me_intm%fhand)
         
         if (typ_prc.eq.optinf_prc_spinp) then
            call spin_project(me_intm,me_result,
     &              fspc(1),lenme,
     &              xbuf1,xbuf2,.true.,xnrm,
     &              opti_info,orb_info,
     &              op_info,str_info,strmap_info)
         elseif (typ_prc.eq.
     &              optinf_prc_spinrefp)then
               call spin_project(me_intm,me_result,
     &              fspc(2),lenme,
     &              xbuf1,xbuf2,.true.,xnrm,
     &              opti_info,orb_info,
     &              op_info,str_info,strmap_info)

               call evaluate2(fspc(1),.false.,.false.,
     &              op_info,str_info,strmap_info,orb_info,
     &              xnrm,.false.)

            else
c     dbg
c     print *,"iopt=",iopt
c     call print_list('projected vector:',
c     &                 me_scr(iopt)%mel,"NORM",
c     &                 -1d0,0d0,
c     &                 orb_info,str_info)
c     dbgend
               call evaluate2(fspc(1),.false.,.false.,
     &              op_info,str_info,strmap_info,orb_info,
     &              xnrm,.false.)
c     dbg
c     call print_list('projected vector:',
c     &                 me_scr(iopt)%mel,"NORM",
c     &                 -1d0,0d0,
c     &                 orb_info,str_info)
c     dbgend
            end if
            if (xnrm.lt.1d-12) call warn('evpc_core',
     &           'Nothing left after projection!')
            call assign_me_list(me_opt%label,
     &           op_opt%name,op_info)
      end if
      case(optinf_prc_blocked)
! not sure, what happens here
         if (nincore.lt.3)
     &        call quit(1,'evpc_core',
     &        'I need at least 3 incore vectors (prc_special)')
         call vec_from_da(me_resid%fhand,iroot,xbuf1,
     &        lenme)
         call dscal(lenme,1d0/xnrm,xbuf1,1)
         xshf = -xeig
         call optc_prc_special2(me_mvp(iopt)%mel,me_special,
     &        nspecial,
     &        op_opt%name,xshf,
     &        nincore,xbuf1,xbuf2,xbuf3,lbuf,
     &        orb_info,op_info,str_info,strmap_info)
         call vec_to_da(me_result%fhand,iroot,xbuf1,
     &        lenme)
      case default
         call quit(1,'evpc_core','unknown preconditioner type')
      end select
      return
      end subroutine

*----------------------------------------------------------------------*
!>    norms a vector of nlists me_lists
*----------------------------------------------------------------------*
c$$$      subroutine orthogonalize_records(
c$$$     &     me_lists, me_Svlists, me_Mvlists, nlists,
c$$$     &     zero_nrm, use_s,
c$$$     &     xbuf1, xbuf2, lbuf)
c$$$      implicit none
c$$$
c$$$      do irec=1,nrec
c$$$         do jrec=1,irec
c$$$            do ilist=1,nlists
c$$$               lenlist=me_lists(ilist)%mel%len_op
c$$$               ffme=>me_lists(ilist)%mel%fhand
c$$$               
c$$$               call vec_from_da(me_lists(ilist)%mel%fhand,irec,xbuf1,
c$$$     &              lenlist)
c$$$               call vec_from_da(me_clists(ilist)%mel%fhand,jrec,xbuf2,
c$$$     &              lenlist)
c$$$               overlapp=overlapp+ddot(lenlist,xbuf1,1,xbuf2,1)
c$$$            end do
c$$$            do ilist=1,nlists
c$$$               call vec_from_da(ffme,irec,xbuf1,lenlist)
c$$$               call vec_from_da(ffme,jrec,xbuf2,lenlist)
c$$$               do idx=1,lenlist
c$$$                  
c$$$               
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      end subroutine

      subroutine rescale_vecs(me_lists, me_clists, me_Mvlists, nlists,
     &     zero_nrm, use_s,
     &     xbuf1, xbuf2, lbuf)
      implicit none
      include 'par_scale_copy_modes.h'

      
      interface
         subroutine mel_scale_copy(
     &     me_inp, me_res,
     &     buf1,buf2, lbuf, nbuf,
     &     fac, nfac,
     &     mode,
     &        opti_info)
         import optimize_info,me_list
         integer,intent(in)::
     &        lbuf, nbuf, nfac,
     &        mode
         type(me_list), intent(inout)::
     &        me_inp, me_res
         type(optimize_info),intent(in)::
     &        opti_info
         real(8), intent(inout) ::
     &        buf1(:), buf2(:)
         real(8),intent(in)::
     &        fac(nfac)
         end subroutine
      end interface
      
      real(8),parameter::
     &     zero_norm2_thresh=1.0D-12
      integer, parameter::
     &     ntest = 30
      character(len=*),parameter::
     &     i_am="rescale_vecs"

      integer,intent(in)::
     &     nlists,
     &     lbuf
      type(me_list_array),intent(inout)::
     &     me_lists(*), me_clists(*),me_MVlists(*)
      real(8),intent(inout)::
     &     xbuf1(lbuf), xbuf2(lbuf)
      logical,intent(out)::
     &     zero_nrm
      logical,intent(in)::
     &     use_s(nlists)
      
      integer::
     &     ilist,
     &     ifree
      real(8)::
     &     xnrm2
      logical::
     &     use_metric
      type(filinf),pointer::
     &     ffme

      real(8),external::
     &     me_ddot
      type(optimize_info)::
     &     opti_info2           

      use_metric = .False.
      do ilist=1,nlists
         use_metric = use_metric .or. use_s(ilist)
      end do
      xnrm2=0
      do ilist=1,nlists
         if(use_s(ilist))then
            xnrm2=xnrm2+
     &           me_ddot(me_lists(ilist)%mel, me_clists(ilist)%mel,
     &           xbuf1, xbuf2, 2, lbuf)
         else
            xnrm2=xnrm2+
     &           me_ddot(me_lists(ilist)%mel, me_lists(ilist)%mel,
     &           xbuf1, xbuf2, 2, lbuf)
         end if
      end do

      print *,"norm**2 was", xnrm2
      zero_nrm = xnrm2.le.zero_norm2_thresh

      if(zero_nrm .and. ntest.gt.20)
     &     write(lulog,*) "detected vector with norm**2 of",xnrm2

      if(zero_nrm) return
      
         
      ifree = mem_alloc_int(opti_info2%nsec,1,'nsec')
      ifree = mem_alloc_int(opti_info2%nwfpsec,1,'nwfpsec')
      ifree = mem_alloc_int(opti_info2%idstsec,1,'idstsec')
      ifree = mem_alloc_real(opti_info2%signsec,1,'signsec')
      opti_info2%nsec(1) = 1
      opti_info2%idstsec(1) = 1
      opti_info2%signsec(1) = 1d0

      do ilist=1,nlists
         opti_info2%nwfpsec(1) = me_lists(ilist)%mel%len_op
         call mel_scale_copy(
     &        me_lists(ilist)%mel, me_lists(ilist)%mel,
     &        xbuf1,xbuf2, lbuf, 2,
     &        (/ 1/sqrt(xnrm2) /), 1,
     &        MOD_SCALE,
     &        opti_info2)
         if(use_metric)then
            call mel_scale_copy(
     &           me_clists(ilist)%mel, me_clists(ilist)%mel,
     &           xbuf1,xbuf2, lbuf, 2,
     &           (/ 1/sqrt(xnrm2) /), 1,
     &           MOD_SCALE,
     &           opti_info2)
         end if
          call mel_scale_copy(
     &        me_Mvlists(ilist)%mel, me_Mvlists(ilist)%mel,
     &        xbuf1,xbuf2, lbuf, 2,
     &        (/ 1/sqrt(xnrm2) /), 1,
     &        MOD_SCALE,
     &        opti_info2)
       end do
       end subroutine
      
      
*----------------------------------------------------------------------*
!>    checks if the energy change between iterations is small enough to count as converged
!!
!!    @param iter iteration counter
!!    @param old_eig old eigenvalue of current root
!!    @param new_eig new_eigenvalue of current root
!!    @param thrgrd_e threshold for energy change
*----------------------------------------------------------------------*
      pure function check_e_convergence(iter,old_eig,new_eig, thrgrd_e)
*----------------------------------------------------------------------*
      implicit none
      logical :: check_e_convergence
      integer, intent(in)::
     &     iter

      real(8),intent(in)::
     &     old_eig,new_eig,
     &     thrgrd_e
      check_e_convergence=.true.
      if (abs(old_eig-new_eig).gt.thrgrd_e)
     &     check_e_convergence=.false.
      return
      end function


*----------------------------------------------------------------------*
!>  warns that energy changed too much
!!
!!
*----------------------------------------------------------------------*
      subroutine warn_e_convergence(lu, iroot, old_eig, new_eig)
      implicit none
      !iprlvl from include in paren routine
      integer,intent(in)::
     &     lu,
     &     iroot
      real(8),intent(in)::
     &     old_eig,
     &     new_eig

      character(len=*),parameter::
     &     i_am="davidson_driver:check_convergence"

      call warn (i_am,"root order may have changed")
      if (iprlvl.ge.5)then
         write(lu,*) "root no.:",iroot
         write(lu,*) "old,new:" , old_eig,new_eig
      end if
      end subroutine



*----------------------------------------------------------------------*
!!   check if the norm of the residual is below threshold
!>
*----------------------------------------------------------------------*
      pure function check_r_convergence( xrsnrm, thrgrd_r, iroot, nopt)
      implicit none
      logical :: check_r_convergence
      integer,intent(in)::
     &     iroot, nopt
      real(8),dimension(:,:),intent(in)::
     &     xrsnrm
      real(8),dimension(nopt),intent(in)::
     &     thrgrd_r
      integer::
     &     iopt
      check_r_convergence=.true.
      do iopt=1, nopt
         check_r_convergence=  check_r_convergence.and.
     &        (xrsnrm(iroot,iopt) .le. thrgrd_r(iopt) )
      end do
      return
      end function
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      pure function mel_get_maxrec(mel)
      implicit none
      integer:: mel_get_maxrec
      type(me_list),intent(in)::mel
      mel_get_maxrec=mel%fhand%active_records(2)

      end function
      pure function dvdsbsp_get_nfree(dvdsbsp)

      integer::
     &     dvdsbsp_get_nfree

      type(davidson_subspace_t),intent(in)::
     &     dvdsbsp

      dvdsbsp_get_nfree=dvdsbsp%nmaxsub-dvdsbsp%ncursub
      end function
      end subroutine
