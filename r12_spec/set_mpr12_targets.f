*----------------------------------------------------------------------*
      subroutine set_mpr12_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set up targets and dependencies for MP(n)-R12 calculations
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_actions.h'
      include 'def_orbinf.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type

      integer ::
     &     mode
      logical ::
     &     fixed

      call set_ccmp_general_targets(tgt_info,orb_info)

      call get_argument_value('method.R12','fixed',lval=fixed)
      call get_argument_value('method.R12','r12op',ival=mode)

      if (fixed.or.mode.gt.0) then
        call set_r12f_general_targets(tgt_info,orb_info,env_type)
        call set_mpr12f_special_targets(tgt_info,orb_info)
      else
        call set_r12_general_targets(tgt_info,orb_info,env_type)
        call set_mpr12_special_targets(tgt_info,orb_info)
      end if

c      call set_ccr12_debug_targets(tgt_info,orb_info)

      return
      end
