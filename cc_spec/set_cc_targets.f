*----------------------------------------------------------------------*
      subroutine set_cc_targets(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set up targets and dependencies for CC calculations
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

      call set_ccmp_general_targets(tgt_info,orb_info)
c      call set_cc_general_targets(tgt_info,orb_info)
      call set_cc_special_targets(tgt_info,orb_info)
      call set_cc_pt_targets(tgt_info,orb_info,env_type)
      call set_cc_gsrsp_targets(tgt_info,orb_info)
      call set_cc_exst_targets(tgt_info,orb_info)
      call set_cc_ipst_targets(tgt_info,orb_info)

      return
      end
