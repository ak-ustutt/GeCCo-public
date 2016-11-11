
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
            call switch_mel_record(me_mvp(iopt)%mel,iroot)
            call switch_mel_record(me_mvort(iopt)%mel,iroot)
            if (use_s(iopt))then
               call switch_mel_record(me_met(iopt)%mel,iroot)
               call switch_mel_record(me_metort(iopt)%mel,iroot)
            end if
            if (trafo(iopt)) then
               
               call transform_forward_wrap(flist,depend,
     &              me_special,me_mvp,me_mvort, !mvp-> mvort
     &              xrsnrm, 
     &              nroot, iroot, iopt, iroot, nspecial,
     &              me_opt,
     &              op_info, str_info, strmap_info, orb_info, opti_info)
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                 call set_blks(me_mvort(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
              endif
              ! Not yet sure how to treat the metric
              if (use_s(iopt))then
                 call transform_forward_wrap(flist,depend,
     &                me_special,me_met,me_metort, !met-> metort
     &                xrsnrm, nroot, 
     &                iroot, iopt, iroot, nspecial,
     &                me_opt,
     &                op_info,str_info,strmap_info, orb_info, opti_info)
                 if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                    call set_blks(me_mvort(iopt)%mel,
     &                   "P,H|P,V|V,H|V,V",0d0)
                 end if
              else
                 me_metort(iopt)%mel=> null()
              end if
           else
              if (use_s(iopt))then
                 me_metort(iopt)%mel=> me_met(iopt)%mel
              else
                 me_metort(iopt)%mel=> null()
              end if
           end if

         end do !iopt
         call dvdsbsp_update(dvdsbsp,
     &        me_metort,
     &        me_mvort,nopt,    !mvort into dvdsbsp and update the vMv-matrix
     &        xbuf1, xbuf2, xbuf3, nincore, lenbuf)
         
      end do

      if (nroot.gt.dvdsbsp_get_nfree(dvdsbsp) )then
         call dvdsbsp_compress(dvdsbsp, nroot,
     &        nopt, me_scr,     !scr as scratch
     &        xbuf1, xbuf2, lenbuf ,nincore)
      end if

      
      call davidson_assemble_residuals(dvdsbsp,
     &     leig, nnew,
     &     me_scr, nopt, nroot, xrsnrm, !temporary assemble on scr
     &     xbuf1, xbuf2, nincore, lenbuf)
      
      irecres=0
      do iroot=1,nnew
         conv=.true.
         if (iter.gt.1
     &        .and. abs(leig(iroot)-xeig(iroot)).gt.thrgrd_e)then !energy changing?
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
               call switch_mel_record(me_scr(iopt)%mel,iroot)
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
               call switch_mel_record(me_opt(iopt)%mel,iroot)
            end do 
            call vec_normalize(me_vort, nopt, xbuf1, xbuf2, lenbuf)
            do iopt=1,nopt
               if (trafo(iopt) ) then
                  call transform_back_wrap(flist,depend,
     &                 me_special,me_vort,me_opt, !vort -> opt !
     &                 iroot, iopt, nspecial,
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
         do iopt=1,nopt
            call switch_mel_record(me_res(iopt)%mel,iroot)
            call switch_mel_record(me_vort(iopt)%mel,iroot)
         end do

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
         if (ntest.gt.100)then
            do iopt=1,nopt
               write(lulog,*) "root no.",iroot
               call print_list("unprojected trv",me_trv(iopt)%mel,
     &              "LIST",0d0,0d0,
     &              orb_info,str_info)
            end do
         end if
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
            call switch_mel_record(me_vort(iopt)%mel,iroot)
            call switch_mel_record(me_trv(iopt)%mel,iroot)
            if (trafo(iopt) ) then
               call transform_back_wrap(flist,depend,
     &              me_special,me_vort,me_trv, !vort -> trv !new_trialvector created
     &              iroot, iopt,nspecial,
     &              me_opt,
     &              op_info, str_info, strmap_info, 
     &              orb_info, opti_info)
!            else me_vort => me_trv
            end if
         end do
      end do                    !iopt
      return
      contains

c$$$!######################################################################
c$$$! subroutines for the transformation 
c$$$!######################################################################
c$$$*----------------------------------------------------------------------*
c$$$!>    wrapper for forward transformation to encapsulate some stupid decisions
c$$$!!
c$$$!!
c$$$*----------------------------------------------------------------------*
c$$$      subroutine transform_forward_wrap(flist,depend,
c$$$     &     me_special,me_in,me_out,
c$$$     &     xrsnrm, 
c$$$     &     nroot, iroot, iopt, irecscr,
c$$$     &     me_tgt,
c$$$     &     op_info, str_info, strmap_info, orb_info, opti_info)
c$$$*----------------------------------------------------------------------*
c$$$      implicit none
c$$$
c$$$      type(me_list_array), dimension(*)::
c$$$     &     me_special,me_in,me_out, me_tgt
c$$$      type(formula_item),intent(in)::
c$$$     &     flist
c$$$      type(dependency_info),intent(in)::
c$$$     &     depend
c$$$      integer, intent(in)::
c$$$     &     nroot,
c$$$     &     iroot, 
c$$$     &     iopt,
c$$$     &     irecscr
c$$$
c$$$      real(8), Dimension(nroot,*), intent(inout)::
c$$$     &     xrsnrm
c$$$
c$$$      type(orbinf), intent(in) ::
c$$$     &     orb_info
c$$$      type(operator_info), intent(inout) ::
c$$$     &     op_info
c$$$      type(strinf), intent(in) ::
c$$$     &     str_info
c$$$      type(strmapinf) ::
c$$$     &     strmap_info
c$$$      type(optimize_info)::
c$$$     &     opti_info
c$$$
c$$$
c$$$      type(me_list),pointer::
c$$$     &     me_trf
c$$$      type(operator),pointer::
c$$$     &     op_in,
c$$$     &     op_trf
c$$$      real(8) ::
c$$$     &     xnrm
c$$$
c$$$      if (nspecial.ge.3)then 
c$$$         me_trf=> me_special(3)%mel
c$$$         op_trf=> me_special(3)%mel%op
c$$$      else
c$$$         me_trf=> null() ! leave it associated as it is now
c$$$         op_trf=> null() ! --
c$$$      end if
c$$$
c$$$      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
c$$$         op_in => me_special(4)%mel%op
c$$$      else
c$$$         op_in => me_special(1)%mel%op
c$$$      endif
c$$$
c$$$
c$$$
c$$$      call switch_mel_record(me_out(iopt)%mel,irecscr)
c$$$      call  switch_mel_record(me_in(iopt)%mel,irecscr)
c$$$      
c$$$c dbg     
c$$$c      call print_list('residual vector before transformation:',
c$$$c     &     me_in,"LIST",
c$$$c     &     -1d0,0d0,
c$$$c     &     orb_info,str_info)
c$$$c dbg end 
c$$$
c$$$      call change_basis_old(flist, depend,
c$$$     &     me_in(iopt)%mel, op_in,
c$$$     &     me_out(iopt)%mel, me_out(iopt)%mel%op, xnrm,
c$$$     &     me_trf, op_trf,                         
c$$$     &     me_tgt(iopt)%mel,
c$$$     &     op_info, str_info, strmap_info, orb_info)
c$$$c dbg
c$$$c            call print_list('transformed residual vector:',
c$$$c     &           me_scr(iopt)%mel,"LIST",
c$$$c     &           -1d0,0d0,
c$$$c     &           orb_info,str_info)
c$$$c dbgend
c$$$      xrsnrm(iroot,iopt) = xnrm
c$$$      return
c$$$      end subroutine
c$$$*----------------------------------------------------------------------*
c$$$!>    wrapper for back transformateion to encapsulate some stupid decisions
c$$$!!
c$$$!!
c$$$*----------------------------------------------------------------------*
c$$$      subroutine transform_back_wrap(flist,depend,
c$$$     &     me_special, me_in,me_out, 
c$$$     &     iroot, iopt,
c$$$     &     me_tgt,
c$$$     &     op_info, str_info, strmap_info, orb_info, opti_info)
c$$$*----------------------------------------------------------------------*
c$$$      implicit none
c$$$
c$$$      type(me_list_array), dimension(*)::
c$$$     &     me_special,me_in, me_out, me_tgt
c$$$      type(formula_item),intent(in)::
c$$$     &     flist
c$$$      type(dependency_info),intent(in)::
c$$$     &     depend
c$$$      integer, intent(in)::
c$$$     &     iroot, 
c$$$     &     iopt
c$$$
c$$$
c$$$      type(orbinf), intent(in) ::
c$$$     &     orb_info
c$$$      type(operator_info), intent(inout) ::
c$$$     &     op_info
c$$$      type(strinf), intent(in) ::
c$$$     &     str_info
c$$$      type(strmapinf) ::
c$$$     &     strmap_info
c$$$      type(optimize_info)::
c$$$     &     opti_info
c$$$
c$$$
c$$$      type(me_list),pointer::
c$$$     &     me_trf
c$$$      type(operator),pointer::
c$$$     &     op_in,
c$$$     &     op_trf
c$$$      real(8) ::
c$$$     &     xnrm
c$$$
c$$$      if (nspecial.ge.3)then ! who thought it would be a good idea to determine the algorithm by the number of arguments? 
c$$$         me_trf=> me_special(2)%mel
c$$$         op_trf=> me_special(2)%mel%op
c$$$      else
c$$$         me_trf=> null() ! leave it associated as it is now
c$$$         op_trf=> null() ! --
c$$$      end if
c$$$
c$$$      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
c$$$         op_in => me_special(4)%mel%op
c$$$      else
c$$$         op_in => me_special(1)%mel%op
c$$$      endif
c$$$
c$$$
c$$$      call  switch_mel_record(me_in(iopt)%mel,iroot)
c$$$
c$$$      call switch_mel_record(me_out(iopt)%mel,iroot)
c$$$      
c$$$c dbg     
c$$$c      call print_list('trial vector before back transformation:',
c$$$c     &     me_in,"LIST",
c$$$c     &     -1d0,0d0,
c$$$c     &     orb_info,str_info)
c$$$c dbg end 
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      call change_basis_old(flist, depend,
c$$$     &     me_in(iopt)%mel, op_in,
c$$$     &     me_out(iopt)%mel, me_out(iopt)%mel%op, xnrm,
c$$$     &     me_trf, op_trf,                         ! 
c$$$     &     me_tgt(iopt)%mel,
c$$$     &     op_info, str_info, strmap_info, orb_info)
c$$$
c$$$c dbg     
c$$$c      call print_list('trial vector after back transformation:',
c$$$c     &     me_opt(iopt)%mel,"LIST",
c$$$c     &     -1d0,0d0,
c$$$c     &     orb_info,str_info)
c$$$c dbg end 
c$$$
c$$$      return
c$$$      end subroutine
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$*----------------------------------------------------------------------*
c$$$!>    subroutine for the transformation into the orthogonal basis
c$$$!!
c$$$!!    uses the old convention where the transformation formula is a subset of the whole formula
c$$$!!    me_tgt and determines the me_list that was bound to the target operator as the dependencies where 
c$$$!!    evaluated
c$$$*----------------------------------------------------------------------*
c$$$      subroutine change_basis_old(flist, depend,
c$$$     &     me_in, op_in, 
c$$$     &     me_out, op_out, outnrm,
c$$$     &     me_trf, op_trf,
c$$$     &     me_tgt,
c$$$     &     op_info, str_info, strmap_info, orb_info)
c$$$*----------------------------------------------------------------------*
c$$$      implicit none
c$$$      character(len=*),parameter::
c$$$     &     i_am="change_basis_old"
c$$$      integer,parameter::
c$$$     &     ntest=1000
c$$$      
c$$$      type(formula_item)::
c$$$     &     flist
c$$$
c$$$      type(me_list)::
c$$$     &     me_in,              !> list to be transformed 
c$$$     &     me_out,             !> result list
c$$$     &     me_tgt             !>list with the transformation operator
c$$$      type(me_list),pointer::     
c$$$     &     me_trf               !> target list definition to specify which subformula should actually be evaluated (stupid design decision)
c$$$      type(operator)::
c$$$     &     op_in,
c$$$     &     op_out
c$$$      type(operator),pointer::
c$$$     &     op_trf
c$$$      real(8),intent(out)::
c$$$     &     outnrm               !norm of the output list
c$$$      type(dependency_info)::
c$$$     &     depend               !>dependency info for the formula 
c$$$      type(orbinf), intent(in) ::
c$$$     &     orb_info
c$$$      type(operator_info), intent(inout) ::
c$$$     &     op_info
c$$$      type(strinf), intent(in) ::
c$$$     &     str_info
c$$$      type(strmapinf) ::
c$$$     &     strmap_info
c$$$
c$$$      integer,dimension(:),allocatable::
c$$$     &     idxselect
c$$$      real(8),dimension(:),allocatable::
c$$$     &     xret
c$$$      integer::
c$$$     &     nselect
c$$$
c$$$      if (ntest.ge.100)then
c$$$         call write_title(lulog,wst_dbg_subr,i_am)
c$$$         write(lulog,*) "out:",me_out%label,"bound to:",op_out%name
c$$$         write(lulog,*) " in:",me_in%label,"bound to:",op_in%name
c$$$         write(lulog,*) "tgt:",me_tgt%label,"bound to:",me_tgt%op%name
c$$$         if (associated(me_trf))
c$$$     &        write(lulog,*) "trf:",me_trf%label
c$$$         if (associated(op_trf))
c$$$     &        write(lulog,*) "trf bound:",op_trf%name
c$$$      end if 
c$$$      call assign_me_list(me_out%label,
c$$$     &     op_out%name, op_info)
c$$$      call assign_me_list(me_in%label,
c$$$     &     op_in%name,op_info)
c$$$      if(associated(me_trf).and. associated(op_trf)) then
c$$$         call assign_me_list(me_trf%label, op_trf%name, op_info)
c$$$      else if(associated(me_trf).or. associated(op_trf))then
c$$$         call quit(1,i_am,
c$$$     &        "please make sure that either both (transformation list"//
c$$$     &        "and operator) or neither is associated") 
c$$$      end if
c$$$
c$$$      allocate(xret(depend%ntargets),idxselect(depend%ntargets))
c$$$      nselect=0
c$$$      call select_formula_target(idxselect,nselect,
c$$$     &     me_tgt%label,depend,op_info)
c$$$! pretend that me_tgt is not up to date
c$$$      call reset_file_rec(me_tgt%fhand)
c$$$      call frm_sched(xret,flist,depend,idxselect, nselect,
c$$$     &     .true.,.false.,op_info,str_info,strmap_info,orb_info)
c$$$      ! actually it stays up to date
c$$$      call touch_file_rec(me_tgt%fhand)
c$$$      outnrm=xret(idxselect(1))
c$$$      deallocate(xret,idxselect)
c$$$      return
c$$$
c$$$      end subroutine

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
     &     ntest=1000
      

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
      pure function dvdsbsp_get_nfree(dvdsbsp)
      
      integer::
     &     dvdsbsp_get_nfree

      type(davidson_subspace_t),intent(in)::
     &     dvdsbsp

      dvdsbsp_get_nfree=dvdsbsp%nmaxsub-dvdsbsp%ncursub
      end function
      end subroutine
