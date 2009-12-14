*----------------------------------------------------------------------*
      subroutine define_r12_hhat(tgt_info,screen)
*----------------------------------------------------------------------*
*     define Hhat operator in the case of r12 response
*
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'par_opnames_gen.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      logical, intent(in) ::
     &     screen

      integer ::
     &     occ_def(ngastp,2,60), ndef
      character(len_command_par) ::
     &     parameters(2)

      occ_def = 0
      ! normal hamiltonian except unnecessary (deexcitation only) blocks
      occ_def(IHOLE,1,2) = 1
      occ_def(IHOLE,2,2) = 1

      occ_def(IPART,1,3) = 1
      occ_def(IHOLE,2,3) = 1

      occ_def(IPART,1,4) = 1
      occ_def(IPART,2,4) = 1

      occ_def(IHOLE,1,5) = 2
      occ_def(IHOLE,2,5) = 2

      occ_def(IHOLE,1,6) = 2
      occ_def(IHOLE,2,6) = 1
      occ_def(IPART,2,6) = 1

      occ_def(IHOLE,1,7) = 1
      occ_def(IHOLE,2,7) = 2
      occ_def(IPART,1,7) = 1

      occ_def(IHOLE,1,8) = 1
      occ_def(IHOLE,2,8) = 1
      occ_def(IPART,1,8) = 1
      occ_def(IPART,2,8) = 1

      occ_def(IHOLE,1,9) = 1
      occ_def(IPART,1,9) = 1
      occ_def(IPART,2,9) = 2

      occ_def(IPART,1,10) = 2
      occ_def(IHOLE,2,10) = 2

      occ_def(IPART,1,11) = 2
      occ_def(IHOLE,2,11) = 1
      occ_def(IPART,2,11) = 1

      occ_def(IPART,1,12) = 2
      occ_def(IPART,2,12) = 2

      ndef = 12

      if (.not.screen) then
        ! blocks with only one external line (except deexcitation only blocks)
        occ_def(IPART,1,13) = 1
        occ_def(IEXTR,2,13) = 1

        occ_def(IEXTR,1,14) = 1
        occ_def(IHOLE,2,14) = 1

        occ_def(IEXTR,1,15) = 1
        occ_def(IPART,2,15) = 1

        occ_def(IHOLE,1,16) = 2
        occ_def(IHOLE,2,16) = 1
        occ_def(IEXTR,2,16) = 1

        occ_def(IHOLE,1,17) = 2
        occ_def(IPART,2,17) = 1
        occ_def(IEXTR,2,17) = 1

        occ_def(IHOLE,1,18) = 1
        occ_def(IHOLE,2,18) = 1
        occ_def(IPART,1,18) = 1
        occ_def(IEXTR,2,18) = 1

        occ_def(IHOLE,1,19) = 1
        occ_def(IPART,2,19) = 1
        occ_def(IPART,1,19) = 1
        occ_def(IEXTR,2,19) = 1

        occ_def(IHOLE,1,20) = 1
        occ_def(IHOLE,2,20) = 2
        occ_def(IEXTR,1,20) = 1

        occ_def(IHOLE,1,21) = 1
        occ_def(IHOLE,2,21) = 1
        occ_def(IEXTR,1,21) = 1
        occ_def(IPART,2,21) = 1

        occ_def(IHOLE,1,22) = 1
        occ_def(IPART,2,22) = 2
        occ_def(IEXTR,1,22) = 1

        occ_def(IPART,1,23) = 2
        occ_def(IHOLE,2,23) = 1
        occ_def(IEXTR,2,23) = 1

        occ_def(IPART,1,24) = 2
        occ_def(IPART,2,24) = 1
        occ_def(IEXTR,2,24) = 1

        occ_def(IPART,1,25) = 1
        occ_def(IHOLE,2,25) = 2
        occ_def(IEXTR,1,25) = 1

        occ_def(IPART,1,26) = 1
        occ_def(IHOLE,2,26) = 1
        occ_def(IEXTR,1,26) = 1
        occ_def(IPART,2,26) = 1

        occ_def(IPART,1,27) = 1
        occ_def(IPART,2,27) = 2
        occ_def(IEXTR,1,27) = 1

        ndef = 27
      end if
      call op_from_occ_parameters(-1,parameters,2,
     &                            occ_def,ndef,1,(/     0,     0/),ndef)
      call set_rule(op_hhat,ttype_op,DEF_OP_FROM_OCC,
     &              op_hhat,1,1,
     &              parameters,2,tgt_info)

      return
      end
