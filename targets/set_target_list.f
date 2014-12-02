*----------------------------------------------------------------------*
      subroutine set_target_list(tgt_info,orb_info,env_type,
     &     name_infile,fforbinf)
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
      type(filinf), intent(in) ::
     &     fforbinf
      character(*), intent(in) ::
     &     name_infile

      type(filinf) ::
     &     fftgt
      integer ::
     &     len
      character(len=256) ::
     &     gecco_path

      call get_environment_variable( "GECCO_DIR", value=gecco_path,
     &     length = len)
      if (len.EQ.0)
     &     call quit(1,'set_mr_targets',
     &     "Please, set the GECCO_DIR environment variable.")
      ! Dummy python interface for warnings and general checks
      call set_python_targets(tgt_info,
     &     trim(gecco_path)//'/interfaces/setting_up_python_interface',
     &     name_infile,fforbinf%name)

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
      if (is_keyword_set('method.response').gt.0)
     &    call set_response_targets(tgt_info,orb_info,env_type)

      ! multireference section
      if (is_keyword_set('method.MR').gt.0)
     &     call set_mr_targets(tgt_info,orb_info,env_type,
     &     name_infile,fforbinf)

      ! experimental section
      if (is_keyword_set('calculate.experimental').gt.0) then
        call set_experimental_targets(tgt_info,orb_info,env_type)
      end if

      ! Interfaces
      if (is_keyword_set('calculate.interfaces').gt.0) then
        call set_interface_targets(tgt_info,name_infile,fforbinf%name)
      end if

      call set_idx4deps(tgt_info)

      if (iprlvl.ge.20) then
        call print_target_list(lulog,tgt_info)
      end if
      
      call file_init(fftgt,'target_list',ftyp_sq_frm,0)
      call file_open(fftgt)
      call print_target_list(fftgt%unit,tgt_info)
      call file_close_keep(fftgt)

      return
      end
