
*----------------------------------------------------------------------*
!>       driver routines for the transformation between orthogonal and normal basis
!!
!######################################################################
! subroutines for the transformation 
!######################################################################
*----------------------------------------------------------------------*
!>    wrapper for forward transformation to encapsulate some stupid decisions
!!
!!    transforms a list into the orthogonal basis
*----------------------------------------------------------------------*
      subroutine transform_forward_wrap(flist,depend,
     &     me_special,me_in,me_out,
     &     xrsnrm, 
     &     nroot, iroot, iopt, irecscr, nspecial,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info, opti_info)
*----------------------------------------------------------------------*
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
      integer, parameter::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="transform_forward_wrap"
 
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
     &     irecscr,
     &     nspecial

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
     &     op_trf,
     &     op_out,
     &     op_in_save
      
      real(8) ::
     &     xnrm
      logical::
     &     trf
      
      if (nspecial.ge.3)then 
         me_trf=> me_special(3)%mel
         op_trf=> me_special(3)%mel%op
         trf=.true.
      else
         me_trf=> null() ! leave it associated as it is now
         op_trf=> null() ! --
         trf=.false.
      end if
      if (associated(op_trf))
     &     write(lulog,*) "trf:",op_trf%name
      if (associated(me_trf))
     &     write(lulog,*) "trf:",me_trf%label

      
      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
         op_in => me_special(4)%mel%op
      else
         op_in => me_special(1)%mel%op
      endif
      op_out => me_out(iopt)%mel%op
      op_in_save => me_in(iopt)%mel%op

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
     &     me_out(iopt)%mel,  me_tgt(iopt)%mel%op, xnrm,
     &     me_trf, op_trf, trf,                      
     &     me_tgt(iopt)%mel,
     &     op_info, str_info, strmap_info, orb_info)
      
      if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
         call assign_me_list(me_special(4)%mel%label,
     &        op_in%name, op_info)
      else
        call assign_me_list(me_special(1)%mel%label,
     &        op_in%name, op_info)
      endif
      call assign_me_list(me_out(iopt)%mel%label,
     &     op_out%name, op_info)
      call assign_me_list(me_in(iopt)%mel%label,
     &     op_in_save%name, op_info)

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
     &     iroot, iopt, nspecial,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info, opti_info)
*----------------------------------------------------------------------*
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

      type(me_list_array), dimension(*)::
     &     me_special,me_in, me_out, me_tgt
      type(formula_item),intent(in)::
     &     flist
      type(dependency_info),intent(in)::
     &     depend
      integer, intent(in)::
     &     iroot, 
     &     iopt,
     &     nspecial


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
      logical::
     &     trf

      if (nspecial.ge.3)then ! who thought it would be a good idea to determine the algorithm by the number of arguments? 
         me_trf=> me_special(2)%mel
         op_trf=> me_special(2)%mel%op
         trf=.true.
      else
         me_trf=> null() ! leave it associated as it is now
         op_trf=> null() ! --
         trf=.false.
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
     &     me_out(iopt)%mel, me_tgt(iopt)%mel%op , xnrm,
     &     me_trf, op_trf, trf,                   ! 
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
     &     me_trf, op_trf, trf,
     &     me_tgt,
     &     op_info, str_info, strmap_info, orb_info)
*----------------------------------------------------------------------*
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

      character(len=*),parameter::
     &     i_am="change_basis_old"
      integer,parameter::
     &     ntest=0000
      
      type(formula_item)::
     &     flist

      logical, intent(in)::
     &     trf

      type(me_list)::
     &     me_in,              !> list to be transformed 
     &     me_out,             !> result list
     &     me_tgt             !>list with the transformation operator
      type(me_list)::     
     &     me_trf               !> target list definition to specify which subformula should actually be evaluated (stupid design decision)
      type(operator)::
     &     op_in,
     &     op_out
      type(operator)::
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
         if (trf               )
     &        write(lulog,*) "trf:",me_trf%label
         if (trf)
     &        write(lulog,*) "trf bound:",op_trf%name
      end if 
      call assign_me_list(me_out%label,
     &     op_out%name, op_info)
      call assign_me_list(me_in%label,
     &     op_in%name,op_info)
      if(trf) then
         call assign_me_list(me_trf%label, op_trf%name, op_info)
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
