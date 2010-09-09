*----------------------------------------------------------------------*
      subroutine update_metric(me_dia,me_special,nspecial,
     &     flist,depend,orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     update metric, transformation matrices, and preconditioner
*
*     matthias, sep. 2010
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_dependency_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_dia

      type(formula_item), intent(inout) ::
     &     flist
      type(dependency_info) ::
     &     depend

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info


      integer ::
     &     nselect

      integer, pointer ::
     &     idxselect(:)
      real(8), pointer ::
     &     xret(:)

      allocate(xret(depend%ntargets),idxselect(depend%ntargets))

      ! calculate metric (if not up to date)
        nselect = 0
        call select_formula_target(idxselect,nselect,
     &              'ME_DENS',depend,op_info)
        call select_formula_target(idxselect,nselect,
     &              me_special(5)%mel%label,depend,op_info)
        call frm_sched(xret,flist,depend,idxselect,nselect,
     &              op_info,str_info,strmap_info,orb_info)

      ! get half-transform of square root of inverted metric
      ! and projector matrix
      call inv_op(trim(me_special(5)%mel%label),
     &            trim(me_special(6)%mel%label),
     &             'invsqrt',
     &             op_info,orb_info,str_info,strmap_info)
      ! reorder to transformation matrix ...
      call reo_mel(trim(me_special(2)%mel%label),
     &             trim(me_special(6)%mel%label),
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.false.)  ! dirty: reo vtx. 1 --> 3
      ! ... and to adjoint of transformation matrix
      call reo_mel(trim(me_special(3)%mel%label),
     &             trim(me_special(6)%mel%label),
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.true.)   ! dirty: reo vtx. 1 --> 3

c      ! now get projection matrix
c      call inv_op(trim(me_special(5)%mel%label),
c     &            trim(me_special(6)%mel%label),
c     &             'invsqrt',
c     &             op_info,orb_info,str_info,strmap_info)
      ! reorder projector ...
      call reo_mel(trim(me_special(4)%mel%label),
     &             trim(me_special(5)%mel%label),
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.false.)  ! dirty: reo vtx. 1 --> 3

      ! update preconditioner if requested
      if (nspecial.ge.7) then
        call assign_me_list(me_special(2)%mel%label,
     &                     me_special(2)%mel%op%name,op_info)
        nselect = 0
        call select_formula_target(idxselect,nselect,
     &              me_special(7)%mel%label,depend,op_info)
        call frm_sched(xret,flist,depend,idxselect,nselect,
     &              op_info,str_info,strmap_info,orb_info)

        ! put diagonal of Jacobian to preconditioner
        call dia_from_op(trim(me_dia%label),
     &                   trim(me_special(7)%mel%label),
     &                   op_info,str_info,orb_info)
      end if
      deallocate(xret,idxselect)

      return
      end
