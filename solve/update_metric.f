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
      include 'routes.h'

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
     &     xdum, prc_min, prc_impfac
      integer ::
     &     gno, prc_type, prc_iter, idx_jac, project,
     &     idx_jacuni
      logical ::
     &     prc_traf,fix_met
      character(len_opname) ::
     &     dia_label

      if (nspcfrm.lt.2) call quit(1,'update_metric',
     &      'No formula for metric passed.')
      if (nspecial.lt.6) call quit(1,'update_metric',
     &      'Not enough special lists passed.')
      call get_argument_value('method.MR','GNO',ival=gno)
      call get_argument_value('method.MR','project',ival=project)
      call get_argument_value('method.MR','prc_traf',lval=prc_traf)
      call get_argument_value('calculate.routes','fix_met',lval=fix_met)

      ! fix_met option is used to fix the metric matrix throughout the
      ! calculation by hacking into update part of it. This, along with
      ! oldref=T, can be used to fix the metric matrix even for
      ! different calculations. This can be used to get energies for
      ! different field strengths while ignoring the relaxation effect
      ! of metric matrix (To check the analytic values of properties)

      ! FIXME: The routines should not assume any naming
      ! conventions. 
      ! ME_C00 has the initial C0 copied
       if (fix_met)
     &  call assign_me_list('ME_C00','C0',op_info)

      ! calculate metric (if not up to date)
      call evaluate2(fspc(2),.true.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)

      ! assigning  back C0 to its exact ME list
       if (fix_met)
     &  call assign_me_list('ME_C0','C0',op_info)

      if (gno.gt.0) then
        ! perform spin projection?
        if (spinadapt.ge.2)
     &     call spin_prj_list(1d0,me_special(5)%mel,me_special(5)%mel,0,
     &                        xdum,.false.,
     &                        op_info,str_info,strmap_info,orb_info)
        ! transform to GNO basis
        call evaluate2(fspc(3),.false.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)

        if (project.eq.1) then
          if (nspecial.lt.7.or.nspcfrm.lt.6) call quit(1,
     &         'update_metric','not enough arguments for GNO/seq.orth')
          ! set up special matrices for sequential orthogonalization
          ! (a) evaluate formula
          call evaluate2(fspc(6),.true.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)
          ! (b) substract unity from second block
          call add_unity(-1d0,1,me_special(7)%mel,2,orb_info,str_info)
          ! (c) invert
          call inv_op(1,trim(me_special(7)%mel%label),
     &                1,trim(me_special(7)%mel%label),
     &                'pseudoinv',
     &                op_info,orb_info,str_info,strmap_info)
        end if
      end if

      ! get half-transform of square root of inverted metric
      ! and projector matrix
      if (gno.gt.0.and.project.eq.1) then
        call inv_op(2,(/trim(me_special(5)%mel%label),
     &                  trim(me_special(7)%mel%label)/),
     &              1,trim(me_special(6)%mel%label),
     &              'invsqrt',
     &              op_info,orb_info,str_info,strmap_info)
      else
        call inv_op(1,trim(me_special(5)%mel%label),
     &              1,trim(me_special(6)%mel%label),
     &              'invsqrt',
     &              op_info,orb_info,str_info,strmap_info)
      end if

      ! Trafo into GNO basis required?
      if (gno.gt.0) then
        if (nspcfrm.lt.5) call quit(1,'update_metric',
     &             'Not enough formulas for trafo into GNO basis given')
        ! apply GNO trafo to transformation matrix
        call assign_me_list(me_special(6)%mel%label,
     &                     me_special(5)%mel%op%name,op_info)
        call evaluate2(fspc(4),.false.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)
        ! apply GNO trafo to projector
        call assign_me_list(me_special(5)%mel%label,
     &                     me_special(5)%mel%op%name,op_info)
        call evaluate2(fspc(5),.false.,.false.,
     &               op_info,str_info,strmap_info,orb_info,xdum,.false.)
      end if

      ! reorder to transformation matrix ...
      call reo_mel(trim(me_special(2)%mel%label),
     &             trim(me_special(6)%mel%label),.false.,
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.false.)  ! dirty: reo vtx. 1 --> 3
c      ! ... and to adjoint of transformation matrix
c      call reo_mel(trim(me_special(3)%mel%label),
c     &             trim(me_special(6)%mel%label),.false.,
c     &             op_info,str_info,strmap_info,orb_info,
c     &             13,.true.)   ! dirty: reo vtx. 1 --> 3

      ! reorder projector ...
      call reo_mel(trim(me_special(4)%mel%label),
     &             trim(me_special(5)%mel%label),.false.,
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.false.)  ! dirty: reo vtx. 1 --> 3

      ! update preconditioner if requested
      if (prcupdate.or.prc_traf) then
        idx_jac = 7
        if (gno.gt.0.and.project.eq.1) idx_jac = 8
        if (nspcfrm.lt.3) call quit(1,'update_metric',
     &        'No formula for Jacobian passed.')
        if (nspecial.lt.idx_jac) call quit(1,'update_metric',
     &        'Special list for Jacobian missing.')

        call assign_me_list(me_special(2)%mel%label,
     &                     me_special(2)%mel%op%name,op_info)
        call assign_me_list(me_special(idx_jac)%mel%label,
     &                     me_special(idx_jac)%mel%op%name,op_info)

        call evaluate2(fspc(3),.true.,.false.,
     &              op_info,str_info,strmap_info,orb_info,xdum,.false.)
      end if

      if (prc_traf) then
        idx_jacuni = 8
        if (gno.gt.0.and.project.eq.1) idx_jacuni = 9
        if (nspecial.lt.idx_jacuni) call quit(1,'update_metric',
     &        'Special list for unitary matrix missing.')
        call inv_op(1,trim(me_special(idx_jac)%mel%label),
     &              2,(/trim(me_special(idx_jacuni)%mel%label),
     &                  trim(me_special(6)%mel%label)/),
     &              'invdiagmult',
     &              op_info,orb_info,str_info,strmap_info)
        ! reorder (again) to transformation matrix
        call reo_mel(trim(me_special(2)%mel%label),
     &               trim(me_special(6)%mel%label),.false.,
     &               op_info,str_info,strmap_info,orb_info,
     &               13,.false.)  ! dirty: reo vtx. 1 --> 3
      end if

      ! FIXME: No assumptions on list names here (DIAxxxxx,ME_FREF,...)
      if (prcupdate) then
        ! first add inactive elements?
        call get_argument_value('method.MR','prc_type',
     &       ival=prc_type)
        call get_argument_value('method.MR','prc_iter',
     &       ival=prc_iter)
        call get_argument_value('method.MR','prc_min',
     &       xval=prc_min)
        call me_list_label(dia_label,'DIA',1,0,0,0,.false.)
        dia_label = trim(dia_label)//'_T'
        if (prc_type.gt.0) then
          if (prc_type.lt.3) call quit(1,'update_metric',
     &       'preconditioner update only for type 0 and 3')
          call set_prc4op(trim(dia_label),'dia-F',0d0,
     &         (/'ME_FREF'/),1,-huge(1d1),
     &         op_info,str_info,orb_info)
        end if

        ! put diagonal of Jacobian to preconditioner
        if (prc_type.ge.3) then
          call dia_from_op(trim(me_dia%label),
     &                     trim(me_special(idx_jac)%mel%label),'extend',
     &                     op_info,str_info,orb_info)
        else
          call dia_from_op(trim(me_dia%label),
     &                     trim(me_special(idx_jac)%mel%label),'---',
     &                     op_info,str_info,orb_info)
        end if

        ! restrict elements to minimum value?
        call scale_copy_op(trim(dia_label),trim(dia_label),prc_min,1,
     &                     'prc_thresh',0,op_info,orb_info,str_info)

        ! prepare Aoff for iterative improvement?
        if (prc_iter.ge.1) then
          call get_argument_value('method.MR','prc_impfac',
     &         xval=prc_impfac)
          ! change sign of A
          call scale_copy_op(trim(me_special(idx_jac)%mel%label),
     &                       trim(me_special(idx_jac)%mel%label),
     &                       -prc_impfac,1,
     &                       '---',0,op_info,orb_info,str_info)
          ! zero the diagonal => Aoff
          call dia_from_op(trim(me_dia%label), !will not be changed
     &                     trim(me_special(idx_jac)%mel%label),
     &                     'zero_dia',
     &                     op_info,str_info,orb_info)
          ! reorder to form of transformation matrix
          call reo_mel(trim(me_special(idx_jac+1)%mel%label),
     &                 trim(me_special(idx_jac)%mel%label),.true.,
     &                 op_info,str_info,strmap_info,orb_info,
     &                 13,.false.)  ! dirty: reo vtx. 1 --> 3
        end if
      end if

      ! now reorder to adjoint of transformation matrix
      call reo_mel(trim(me_special(3)%mel%label),
     &             trim(me_special(6)%mel%label),.false.,
     &             op_info,str_info,strmap_info,orb_info,
     &             13,.true.)   ! dirty: reo vtx. 1 --> 3

      return
      end
