*----------------------------------------------------------------------*
      subroutine set_cc_targets(tgt_info,orb_info,env_type,
     &     name_infile,fforbinf)
*----------------------------------------------------------------------*
*     set up targets and dependencies for CC calculations
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      include 'stdunit.h'

      include 'ifc_input.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_actions.h'
      include 'def_filinf.h'
      include 'def_orbinf.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      type(orbinf), intent(in) ::
     &     orb_info
      character(*), intent(in) ::
     &     env_type
      ! needed to parse python-based input:
      type(filinf), intent(in) ::
     &     fforbinf
      character(*), intent(in) ::
     &     name_infile

      integer ::
     &     len
      character(len=256) ::
     &     gecco_path, filename

      if (is_keyword_set('method.CC.NEW').gt.0) then
        write(lulog,*) 'entering new setup for CC ...'
        ! branch out for the python based setup of all computations
        call get_environment_variable( "GECCO_DIR", value=gecco_path,
     &                                length = len)

        if (len.eq.0)
     &     call quit(1,'set_cc_targets',
     &     "Please set the GECCO_DIR environment variable.")
        if (len.eq.0) call quit(1,'set_cc_targets','GECCO_DIR not set')

        ! print orbital informations (after get it changed from input)
        call put_orbinfo(orb_info, fforbinf)

        call set_python_targets(tgt_info,
     &       trim(gecco_path)//"/python_spec/sr_cc/ground_state.py",
     &       name_infile,fforbinf%name)

        if (is_keyword_set('calculate.ionization').gt.0) then
          write(lulog,*) 'reading EOMIP targets'
          call set_python_targets(tgt_info,
     &       trim(gecco_path)//"/python_spec/sr_cc/eomip.py",
     &       name_infile,fforbinf%name)
        end if


        write(lulog,*) 'DONE with this section'

        return
      end if
      call set_ccmp_general_targets(tgt_info,orb_info)
c      call set_cc_general_targets(tgt_info,orb_info)
      call set_cc_special_targets(tgt_info,orb_info)
      call set_cc_pt_targets(tgt_info,orb_info,env_type)
      call set_cc_gsrsp_targets(tgt_info,orb_info,env_type)
      call set_cc_exst_targets(tgt_info,orb_info)
      call set_cc_ipst_targets(tgt_info,orb_info)

      return
      end
