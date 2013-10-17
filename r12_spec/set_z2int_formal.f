*----------------------------------------------------------------------*
      subroutine set_z2int_formal(active_orbs,tgt_info)
*----------------------------------------------------------------------*
*     new direct definition of Z2 via targets
*     --- at the moment only for the case when active orbitals are present ---
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'def_orbinf.h'

      include 'ifc_input.h'
      include 'ifc_targets.h'

      include 'par_opnames_gen.h'
      include 'par_formnames_gen.h'
      include 'par_gen_targets.h'
      include 'par_actions.h'

      logical, intent(in) ::
     &     active_orbs
      type(target_info), intent(inout) ::
     &     tgt_info


      if (.not.active_orbs) call quit(0,'set_z2int_formal',
     &    'only for cases with active orbitals present (at the moment)')

      call add_target3([character(len=80) ::
     & 'target Z2INT_R12(',
     & '  depend Z2-INT',
     & '  depend R12',
     & '  depend H',
     & '  #',
     & '  # for Z([HV][HV],PP)',
     & '  #',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=T,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,",  "1,,,V","1,,,VV",     ',
     & '              "5,,PP,HH","5,,PP,H","5,,PP,",',
     & '              "3,5,,P","4,5,,P",',
     & '              "2,3,,[PX]","2,4,,X","3,4,,[PX]") ',
     & '  ) ',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=F,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,",  "1,,,V","1,,,VV",       ',
     & '              "5,,PP,HH","5,,PP,H","5,,PP,",',
     & '              "3,5,,P","4,5,,P",',
     & '              "2,3,,X","2,4,,P","3,4,,X") ',
     & '  )',
     & '  #',
     & '  # for Z([HV][HV],PV)',
     & '  # ',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=F,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,",  "1,,,V","1,,,VV",     ',
     & '              "5,,PV,HH","5,,PV,H","5,,PV,",',
     & '              "3,5,,V","4,5,,P",',
     & '              "2,3,,[PX]","2,4,,X","3,4,,[PX]") ',
     & '  ) ',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=F,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,",  "1,,,V","1,,,VV",       ',
     & '              "5,,PV,HH","5,,PV,H","5,,PV,",',
     & '              "3,5,,V","4,5,,P",',
     & '              "2,3,,X","2,4,,P","3,4,,X") ',
     & '  ) ',
     & '  #',
     & '  # for Z([HV][HV]V,[PV]PV)',
     & '  #',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=F,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,V",  "1,,,VV","1,,,VVV",',
     & '              "5,,[PV]PV,HH","5,,[PV]PV,H","5,,[PV]PV,",',
     & '              "4,5,,P",',
     & '              "2,3,,[PX]","2,4,,X","3,4,,[PX]")',
     & '  )',
     & '  EXPAND_OP_PRODUCT(label=Z2INT_R12,op_res=Z2-INT,new=F,',
     & '   operators=(Z2-INT,R12^+,H,R12,Z2-INT),',
     & '   idx_sv   =(1,2,3,4,1),',
     & '   descr    =("1,,,V",  "1,,,VV","1,,,VVV",',
     & '              "5,,[PV]PV,HH","5,,[PV]PV,H","5,,[PV]PV,",',
     & '              "4,5,,P",',
     & '              "2,3,,X","2,4,,P","3,4,,X")',
     & '  )',
     & '  #',
     & '  # make sure that formula is ordered wrt. result blocks:',
     & '  #',
     & '  REORDER_FORMULA(label_res=Z2INT_R12,label_in=Z2INT_R12)',
     & ')'],tgt_info)

      end

