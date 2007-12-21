*----------------------------------------------------------------------*
      subroutine set_target_list(tgt_info,orb_info)
*----------------------------------------------------------------------*
*     set up targets
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info

      if (is_keyword_set('method.CC').gt.0) then
        call set_cc_targets(tgt_info,orb_info)
      end if

      call set_idx4deps(tgt_info)

      if (iprlvl.ge.10) then
        call print_target_list(tgt_info)
      end if

      return
      end
