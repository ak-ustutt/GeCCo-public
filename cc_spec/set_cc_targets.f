*----------------------------------------------------------------------*
      subroutine set_cc_targets(tgt_info,orb_info)
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

      call set_cc_operator_targets(tgt_info)
      call set_cc_formula_targets(tgt_info)
      call set_cc_me_list_targets(tgt_info,orb_info)

      return
      end
