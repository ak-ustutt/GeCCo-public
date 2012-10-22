*----------------------------------------------------------------------*
      subroutine update_metric(me_dia,me_special,nspecial,
     &     fspc,nspcfrm,orb_info,op_info,str_info,strmap_info,
     &     prcupdate)
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
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspecial, nspcfrm
      logical, intent(in) ::
     &     prcupdate
      type(me_list_array), intent(inout) ::
     &     me_special(nspecial)
      type(me_list), intent(in) ::
     &     me_dia

      type(formula_item), intent(in) ::
     &     fspc(nspcfrm)
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

      real(8) ::
     &     xdum
      integer ::
     &     gno

      if (nspcfrm.lt.2) call quit(1,'update_metric',
     &      'No formula for metric passed.')
      if (nspecial.lt.6) call quit(1,'update_metric',
     &      'Not enough special lists passed.')
      call get_argument_value('method.MR','GNO',ival=gno)

      ! calculate metric (if not up to date)
      call evaluate2(fspc(2),.true.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)

      ! get half-transform of square root of inverted metric
      ! and projector matrix
      call inv_op(trim(me_special(5)%mel%label),
     &            1,trim(me_special(6)%mel%label),
     &            'invsqrt',
     &            op_info,orb_info,str_info,strmap_info)

      ! Trafo into GNO basis required?
      if (gno.gt.0) then
        if (nspcfrm.lt.4) call quit(1,'update_metric',
     &             'Not enough formulas for trafo into GNO basis given')
        ! apply GNO trafo to transformation matrix
        call assign_me_list(me_special(6)%mel%label,
     &                     me_special(5)%mel%op%name,op_info)
        call evaluate2(fspc(3),.false.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)
        ! apply GNO trafo to projector
        call assign_me_list(me_special(5)%mel%label,
     &                     me_special(5)%mel%op%name,op_info)
        call evaluate2(fspc(4),.false.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)
      end if

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

      ! reorder projector ...
      call reo_mel(trim(me_special(4)%mel%label),
     &             trim(me_special(5)%mel%label),
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.false.)  ! dirty: reo vtx. 1 --> 3

      ! update preconditioner if requested
      if (prcupdate) then
        if (nspcfrm.lt.3) call quit(1,'update_metric',
     &        'No formula for Jacobian passed.')
        if (nspecial.lt.7) call quit(1,'update_metric',
     &        'Special list for Jacobian missing.')

        call assign_me_list(me_special(2)%mel%label,
     &                     me_special(2)%mel%op%name,op_info)
        call assign_me_list(me_special(7)%mel%label,
     &                     me_special(7)%mel%op%name,op_info)

        call evaluate2(fspc(3),.true.,.false.,
     &              op_info,str_info,strmap_info,orb_info,xdum,.false.)

        ! put diagonal of Jacobian to preconditioner
        call dia_from_op(trim(me_dia%label),
     &                   trim(me_special(7)%mel%label),.false.,.false.,
     &                   op_info,str_info,orb_info)
      end if

      return
      end
