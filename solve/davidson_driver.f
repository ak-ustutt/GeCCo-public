
!>    driver for the davidson algoritmus
!!    
!!    Davidson      
!!    on entry: Mv-product on me_mvp me
      
      subroutine davidson_driver(
     &     dvdsbsp,
     &     iter,  xrsnrm , xeig,
     &     use_s, update,
     &     me_opt,me_trv,me_mvp,me_dia,me_met,me_scr,me_ext,
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

      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer, intent(inout) ::
     &     iter
      integer, intent(in) ::
     &     nincore, lenbuf, nspecial, nspcfrm

      type(me_list_array), intent(in) ::
     &     me_opt(*), me_dia(*),
     &     me_mvp(*), me_special(*), me_scr(*),me_ext(*)
      type(me_list_array), intent(inout) ::
     &     me_met(*), me_trv(*)


      type(optimize_info), intent(in) ::
     &     opti_info
      type(optimize_status), intent(inout), target ::
     &     opti_stat
      
      real(8), intent(inout) ::
     &     xrsnrm(opti_info%nroot,opti_info%nopt),
     &     xeig(opti_info%nroot,2)

      real(8), intent(inout) ::
     &     xbuf1(*), xbuf2(*), xbuf3(*)
      logical,intent(inout)::
     &     use_s(*),
     &     update(*)
      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend
      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)


      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      

*----------------------------------------------------------------------*
!     locals
 
      integer::
     &     iopt, nopt,
     &     iroot, nroot,
     &     irecscr
      real(8)::
     &     xnrm                !temporary variable for some norms
      logical ::
     &     trafo                !is a transformation necessary for this optimization?
      
      type(me_list),pointer::
     &     me_mvort, me_vort     ! lists for the mv-product and vector in the orthogonal space
      real(8),external::
     &     xnormop

      
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,'evpc_core entered')
      if (ntest.ge.100) write (lulog,*) "iteration ",iter

      nopt = opti_info%nopt
      nroot = opti_info%nroot
      
      irecscr = 1
      do iroot=1,nroot
         do iopt=1,nopt
            
            if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
               me_vort => me_special(1)%mel
               trafo = .true.
            elseif (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
               if (nspecial .lt. 4) call quit(1,'evpc_core',
     &              'TR0 with <4 me_lists')
               me_vort => me_special(4)%mel
               trafo = .true.
            else
               me_vort => me_trv(iopt)%mel
               trafo = .false.
            end if
            if (trafo) then
               
               call transform_forward_wrap(flist,depend,
     &              me_special,me_scr,me_mvp,
     &              xrsnrm, nroot, 
     &              iroot, iopt, iroot,
     &              op_info, str_info, strmap_info, orb_info, opti_info)
               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                  call set_blks(me_scr(iopt)%mel,"P,H|P,V|V,H|V,V",0d0)
               endif
               
               call transform_forward_wrap(flist,depend,
     &              me_special,me_vort,me_trv,
     &              xrsnrm, nroot,
     &              iroot, iopt, iroot,
     &              op_info, str_info, strmap_info, orb_info, opti_info)

               if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
                  call set_blks(me_vort,"P,H|P,V|V,H|V,V",0d0)
                  xnrm=xnormop(me_vort)
                  print *,"dbg: backtransformed trial_vector:",xnrm
               else
                  print *, "dbg: backtransformed trial_vector:",
     &                 xrsnrm(iroot,iopt)
               endif
            end if

            call dvdsbsp_update(dvdsbsp,
     &           me_vort,
     &           me_scr,nopt,
     &           xbuf1, xbuf2, xbuf3, nincore, lenbuf)

         end do
      end do
      call davidson_update_residuals(dvdsbsp,
     &     xeig,
     &     me_scr, nopt, nroot, xrsnrm,
     &     update,
     &     xbuf1, xbuf2, nincore, lenbuf) 

    
!     not yet converged? increase record counter
!      if (.not.is_converged(xrsnrm, iroot, nroot, nopt, 
!     &     opti_info%thrgrd)
!     &     .and. iter.lt.opti_info%maxmacit ) then
        
!      end if
         
      contains
      pure function is_converged(xrsnrm, iroot, nroot , nopt, thrsh)
*----------------------------------------------------------------------*
      implicit none
      logical::
     &     is_converged

      integer, intent(in) ::
     &     iroot, nroot, nopt 

      real(8),dimension(nroot,nopt),intent(in)::
     &     xrsnrm

      real(8),dimension(nopt),intent(in)::
     &     thrsh
      integer ::
     &     iopt
      is_converged=.true.

      do iopt=1,nopt
         is_converged = is_converged
     &        .and. ( xrsnrm(iroot,iopt).lt.thrsh(iopt) )
      end do
      return
      end function
      

      end subroutine
