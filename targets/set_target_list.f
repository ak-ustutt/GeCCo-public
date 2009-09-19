*----------------------------------------------------------------------*
      subroutine set_target_list(tgt_info,orb_info,env_type)
*----------------------------------------------------------------------*
*     set up targets
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'
      include 'mdef_target_info.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type

      type(filinf) ::
     &     fftgt

      call set_general_targets(tgt_info,orb_info,env_type)

      ! normal CC
      if (is_keyword_set('method.CC').gt.0 .and.
     &    is_keyword_set('method.R12').eq.0 ) then
        call set_cc_targets(tgt_info,orb_info,env_type)
      end if

      ! extended CC
      if (is_keyword_set('method.ECC').gt.0 .and.
     &    is_keyword_set('method.R12').eq.0 ) then
        call set_ecc_targets(tgt_info,orb_info)
      end if

      ! MP-R12
      if (is_keyword_set('method.MP').gt.0 .and.
     &    is_keyword_set('method.R12').gt.0 ) then
        call set_mpr12_targets(tgt_info,orb_info,env_type)
      end if

      ! CC-R12
      if (is_keyword_set('method.CC').gt.0 .and.
     &    is_keyword_set('method.R12').gt.0 ) then
        call set_ccr12_targets(tgt_info,orb_info,env_type)
      end if

      call set_r12_test_targets(tgt_info,orb_info,env_type)

      ! response section
      if (is_keyword_set('calculate.response').gt.0)
     &    call set_response_targets(tgt_info,orb_info,env_type)

      ! experimental section
      if (is_keyword_set('calculate.experimental').gt.0) then
        call set_experimental_targets(tgt_info,orb_info,env_type)
      end if

      call set_idx4deps(tgt_info)

      if (iprlvl.ge.20) then
        call print_target_list(luout,tgt_info)
      end if
      
      call file_init(fftgt,'target_list',ftyp_sq_frm,0)
      call file_open(fftgt)
      call print_target_list(fftgt%unit,tgt_info)
      call file_close_keep(fftgt)

      return
      end
