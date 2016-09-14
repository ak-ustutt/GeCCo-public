
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
     &     me_met,me_scr,me_res,
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
     &     ntest = 000
      character(len=*),parameter::
     &     i_am="davidson_driver"
      real(8),parameter::
     &     thrgrd_e=0.1
      
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
     &     me_met(*),me_scr(*),me_res(*),
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
     &     conv
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
      
      do iroot=1,nnew
         do iopt=1,nopt
            

            if (trafo(iopt)) then
               
               call transform_forward_wrap(flist,depend,
     &              me_special,me_mvp,me_mvort, !mvp-> mvort
     &              xrsnrm, nroot, 
     &              iroot, iopt, iroot,
     &              me_opt,
     &              op_info, str_info, strmap_info, orb_info, opti_info)
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                 call set_blks(me_mvort(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
               endif
            end if
         end do
         call dvdsbsp_update(dvdsbsp,
     &        me_vort,
     &        me_mvort,nopt,    !mvort into dvdsbsp and update the vMv-matrix
     &        xbuf1, xbuf2, xbuf3, nincore, lenbuf)
         
      end do
      call davidson_assemble_residuals(dvdsbsp,
     &     leig, nnew,
     &     me_scr, nopt, nroot, xrsnrm, !temporary assemble on scr
     &     xbuf1, xbuf2, nincore, lenbuf)
      
      irecres=0
      do iroot=1,nnew
         conv=.true.
         if (abs(leig(iroot)-xeig(iroot)).gt.thrgrd_e)then !energy changing?
            conv=.false.
            if (iprlvl.ge.5)then
               write(lulog,*) "root order may have changed:"
               write(lulog,*) "root no.:",iroot
               write(lulog,*) "old,new:" , xeig(iroot),leig(iroot)
            end if
         end if 
         do iopt=1,nopt
            conv= conv .and. (xrsnrm(iroot,iopt) .le.
     &           opti_info%thrgrd(iopt) )
         end do
         if (.not.conv)then 
            irecres=irecres+1
            do iopt=1,nopt
               call switch_mel_record(me_res(iopt)%mel,irecres) ! all updated vectors are one after the other 
               call list_copy(me_scr(iopt)%mel,me_res(iopt)%mel, !scr ->res
     &              .false.)    !collecting vectors which will lead to new direction on me_res (me_residual)
            end do
            reig(irecres)=leig(iroot)
         end if 
      end do
      xeig=leig

      nnew=irecres
      

      if (nnew .eq. 0
     &     .or. iter .ge. maxiter) then ! check for end condition
         do iopt=1,nopt
            if(.not.trafo(iopt))me_vort(iopt)%mel => me_opt(iopt)%mel
         end do
         call davidson_assemble_results(dvdsbsp,
     &        xeig,
     &        me_vort, nopt, nroot, lrsnrm, ! results on me_vort
     &        xbuf1, xbuf2, nincore, lenbuf)
         maxvec=min(nroot,mel_get_maxrec(me_vort(1)%mel)) 
         do iroot=1,maxvec
            do iopt=1,nopt
               call switch_mel_record(me_vort(iopt)%mel,iroot)
            end do 
            call vec_normalize(me_vort, nopt, xbuf1, xbuf2, lenbuf)
            do iopt=1,nopt
               if (trafo(iopt) ) then
                  call transform_back_wrap(flist,depend,
     &                 me_special,me_vort,me_opt, !vort -> opt !
     &                 iroot, iopt,
     &                 me_opt,
     &                 op_info, str_info, strmap_info, 
     &                 orb_info, opti_info)
                  
                  
!     else me_vort => me_trv => me_opt
               end if
            end do 
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
      
      do iroot = 1, nnew 
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
         call dvdsbsp_append_vvec(dvdsbsp, !ortho scratchvec to subspace
     &        me_vort,nopt, 
     &        xbuf1, xbuf2, nincore, lenbuf)
      end do                    !iroot
      nnew=dvdsbsp_get_nnew_vvec(dvdsbsp)
      if ( nnew.eq.0)
     &     call quit(0,i_am,
     &     "only linear depended new directions generated")
      
      do iopt=1,nopt
         do iroot=1,nnew
            if (trafo(iopt) ) then
               call transform_back_wrap(flist,depend,
     &              me_special,me_vort,me_trv, !scr -> trv !new_trialvector created
     &              iroot, iopt,
     &              me_opt,
     &              op_info, str_info, strmap_info, 
     &              orb_info, opti_info)
!            else me_vort => me_trv
            end if
         end do
      end do                    !iopt
      return
      contains

!######################################################################
! subroutines for the transformation 
!######################################################################
*----------------------------------------------------------------------*
!>    wrapper for forward transformation to encapsulate some stupid decisions
!!
!!
*----------------------------------------------------------------------*
      subroutine transform_forward_wrap(flist,depend,
     &     me_special,me_in,me_out,
     &     xrsnrm, 
     &     nroot, iroot, iopt, irecscr,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info, opti_info)
*----------------------------------------------------------------------*
      implicit none

      type(me_list_array), dimension(*)::
     &     me_special,me_in,me_out, me_tgt
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
      type(optimize_info)::
     &     opti_info


      type(me_list),pointer::
     &     me_trf
      type(operator),pointer::
     &     op_in,
     &     op_trf
      real(8) ::
     &     xnrm

      if (nspecial.ge.3)then 
         me_trf=> me_special(3)%mel
         op_trf=> me_special(3)%mel%op
      else
         me_trf=> null() ! leave it associated as it is now
         op_trf=> null() ! --
      end if

      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
         op_in => me_special(4)%mel%op
      else
         op_in => me_special(1)%mel%op
      endif



      call switch_mel_record(me_out(iopt)%mel,irecscr)
      call  switch_mel_record(me_in(iopt)%mel,irecscr)
      
c dbg     
c      call print_list('residual vector before transformation:',
c     &     me_in,"LIST",
c     &     -1d0,0d0,
c     &     orb_info,str_info)
c dbg end 

      call change_basis_old(flist, depend,
     &     me_in(iopt)%mel, op_in,
     &     me_out(iopt)%mel, me_out(iopt)%mel%op, xnrm,
     &     me_trf, op_trf,                         
     &     me_tgt(iopt)%mel,
     &     op_info, str_info, strmap_info, orb_info)
c dbg
c            call print_list('transformed residual vector:',
c     &           me_scr(iopt)%mel,"LIST",
c     &           -1d0,0d0,
c     &           orb_info,str_info)
c dbgend
      xrsnrm(iroot,iopt) = xnrm
      return
      end subroutine
*----------------------------------------------------------------------*
!>    wrapper for back transformateion to encapsulate some stupid decisions
!!
!!
*----------------------------------------------------------------------*
      subroutine transform_back_wrap(flist,depend,
     &     me_special, me_in,me_out, 
     &     iroot, iopt,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info, opti_info)
*----------------------------------------------------------------------*
      implicit none

      type(me_list_array), dimension(*)::
     &     me_special,me_in, me_out, me_tgt
      type(formula_item),intent(in)::
     &     flist
      type(dependency_info),intent(in)::
     &     depend
      integer, intent(in)::
     &     iroot, 
     &     iopt


      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(optimize_info)::
     &     opti_info


      type(me_list),pointer::
     &     me_trf
      type(operator),pointer::
     &     op_in,
     &     op_trf
      real(8) ::
     &     xnrm

      if (nspecial.ge.3)then ! who thought it would be a good idea to determine the algorithm by the number of arguments? 
         me_trf=> me_special(2)%mel
         op_trf=> me_special(2)%mel%op
      else
         me_trf=> null() ! leave it associated as it is now
         op_trf=> null() ! --
      end if

      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
         op_in => me_special(4)%mel%op
      else
         op_in => me_special(1)%mel%op
      endif


      call  switch_mel_record(me_in(iopt)%mel,iroot)

      call switch_mel_record(me_out(iopt)%mel,iroot)
      
c dbg     
c      call print_list('trial vector before back transformation:',
c     &     me_in,"LIST",
c     &     -1d0,0d0,
c     &     orb_info,str_info)
c dbg end 





      call change_basis_old(flist, depend,
     &     me_in(iopt)%mel, op_in,
     &     me_out(iopt)%mel, me_out(iopt)%mel%op, xnrm,
     &     me_trf, op_trf,                         ! 
     &     me_tgt(iopt)%mel,
     &     op_info, str_info, strmap_info, orb_info)

c dbg     
c      call print_list('trial vector after back transformation:',
c     &     me_opt(iopt)%mel,"LIST",
c     &     -1d0,0d0,
c     &     orb_info,str_info)
c dbg end 

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
     &     ntest=00
      
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
         write(lulog,*) "out:",me_out%label,"bound to:",op_out%name
         write(lulog,*) " in:",me_in%label,"bound to:",op_in%name
         write(lulog,*) "tgt:",me_tgt%label,"bound to:",me_tgt%op%name
         if (associated(me_trf))
     &        write(lulog,*) "trf:",me_trf%label
         if (associated(op_trf))
     &        write(lulog,*) "trf bound:",op_trf%name
      end if 
      call assign_me_list(me_out%label,
     &     op_out%name, op_info)
      call assign_me_list(me_in%label,
     &     op_in%name,op_info)
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
! pretend that me_tgt is not up to date
      call reset_file_rec(me_tgt%fhand)
      call frm_sched(xret,flist,depend,idxselect, nselect,
     &     .true.,.false.,op_info,str_info,strmap_info,orb_info)
      ! actually it stays up to date
      call touch_file_rec(me_tgt%fhand)
      outnrm=xret(idxselect(1))
      deallocate(xret,idxselect)
      return

      end subroutine

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
         me_intm => me_scr
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
         call switch_mel_record(me_resid,iroot)
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
! not sure, what happens her
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
      subroutine vec_normalize(me_lists, nlists, xbuf1, xbuf2, lbuf)
      implicit none
      integer, parameter::
     &     ntest = 100
      character(len=*),parameter::
     &     i_am="vec_normalize"

      integer,intent(in)::
     &     nlists,
     &     lbuf
      type(me_list_array),intent(in)::
     &     me_lists(*)
      real(8),intent(inout)::
     &     xbuf1(*), xbuf2(*)

      integer::
     &     ilist,
     &     lenlist,
     &     irec,
     &     ii
      real(8)::
     &     xnrm2,xnrm
   
      type(filinf),pointer::
     &     ffme
      
      xnrm2=0
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf1,lenlist)
         do ii=1,lenlist
            xnrm2=xnrm2+xbuf1(ii)**2
         end do
      end do
      xnrm=sqrt(xnrm2)
      print *,"old norm was",xnrm
      
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf1,lenlist)
         xbuf1(1:lenlist)=xbuf1(1:lenlist)/xnrm
         call vec_to_da(ffme,irec,xbuf1,lenlist)
      end do



      end subroutine
*-------------------------------------------------------------*
!>   orthogonalizes all nroot vectors on me_lists()nlists with respect to each other 
!!   and orthogonalizes them
*-------------------------------------------------------------*
      subroutine vecs_orthonorm(me_lists,nlists,nroot,
     &     xbuf1,xbuf2,lbuf)
      implicit none
      integer,parameter::
     &     ntest=00
      type(me_list_array),intent(in)::
     &     me_lists(*)

      integer,intent(in)::
     &     nlists, nroot,
     &     lbuf
      
      real(8),intent(inout)::
     &     xbuf1(*), xbuf2(*)
      
      integer::
     &     iroot,ilist, jroot,kroot,
     &     lenlist
      real(8)::
     &     overlapp(nroot,nroot),
     &     orth_norm2(nroot)     ! squared norm of the orthogonalized vectors
      type(filinf),pointer::
     &     ffme1,ffme2
      real(8)::
     &     tmp_orth_norm2,dbg_nrm
      real(8),external::
     &     ddot
      overlapp(1:nroot,1:nroot)=0.0
      do ilist=1,nlists
         ffme1 => me_lists(ilist)%mel%fhand
         ffme2 => me_lists(ilist)%mel%fhand
         lenlist= me_lists(ilist)%mel%len_op
         do iroot=1,nroot
            call vec_from_da(ffme1,iroot,xbuf1,lenlist)
            do jroot=1,iroot
               call vec_from_da(ffme2,jroot,xbuf2,lenlist)
               overlapp(iroot,jroot) = overlapp(iroot,jroot) +
     &              ddot(lenlist, xbuf1,1, xbuf2 ,1)
            end do
         end do
      end do
     
      if (ntest.ge.100) then
         write(lulog,*) 'trafo-matrix:'
         call wrtmat2(overlapp,nroot,nroot,nroot,nroot)
      end if

!|i'>=|i>-sum_j(|j> <j|i>)
! assuming <i|j>=<j|i>
!<i'|i'>=<i|i>+sum_j(<j|i>(-2<j|i>+sum_k( <i|k><k|j> ) ) )
     
      orth_norm2(1:nroot)=0.0
      do iroot=1, nroot
         do jroot=1,iroot-1
            do kroot=1,jroot    !k <=j or no values in matrix but <k|j> +<j|k>=2<j|k>
               tmp_orth_norm2=tmp_orth_norm2+
     &              2*overlapp(iroot,kroot)*overlapp(jroot,kroot)
            end do
            orth_norm2(iroot)=(tmp_orth_norm2
     &           -2*overlapp(iroot,jroot)
     &           -overlapp(iroot,kroot)*overlapp(jroot,jroot)! added k=j term twice, remove once
     &           )*overlapp(iroot,jroot)
         end do
         orth_norm2(iroot)=orth_norm2(iroot)+overlapp(iroot,iroot)
      end do
      dbg_nrm=0d0
      do iroot=1,nroot
         do ilist=1,nlists
            ffme1 => me_lists(ilist)%mel%fhand
            ffme2 => me_lists(ilist)%mel%fhand
            lenlist= me_lists(ilist)%mel%len_op
            if(ntest.ge.100)
     &           write(lulog,*)"normalizing",me_lists(ilist)%mel%label
            call vec_from_da(ffme1,iroot,xbuf1,lenlist)
            do jroot=iroot+1,nroot
               call vec_from_da(ffme2,jroot,xbuf2,lenlist)
               xbuf2(1:lenlist)=xbuf2(1:lenlist)-
     &              overlapp(iroot,jroot)*xbuf1(1:lenlist)
               call vec_to_da(ffme2,jroot,xbuf2,lenlist)
            end do
            xbuf1(1:lenlist)=xbuf1(1:lenlist)/sqrt(orth_norm2(iroot)) !norm
            call vec_to_da(ffme1,jroot,xbuf1,lenlist)
            if(ntest.ge.100)
     &           dbg_nrm=dbg_nrm+ddot(lenlist,xbuf1,1,xbuf1,1)
         end do
         if(ntest.ge.100)then
            write(lulog, *) "norm of vector: ", iroot,"is",dbg_nrm
            do ilist=1,lenlist
               write(lulog, *) ilist,xbuf1(ilist)
            end do
           end if 
      end do                    !iroot
      return
      end subroutine
*----------------------------------------------------------------------*
      pure function mel_get_maxrec(mel)
      integer:: mel_get_maxrec
      type(me_list),intent(in)::mel
      mel_get_maxrec=mel%fhand%active_records(2)

      end function
      end subroutine
