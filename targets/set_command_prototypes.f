*----------------------------------------------------------------------*
      subroutine set_command_prototypes(tgt_info,env_type)
*----------------------------------------------------------------------*
*     set all commands and possible arguments with their default values
*----------------------------------------------------------------------*
     
      implicit none

      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'ifc_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     env_type

      ! note:
      !   def (for default_provided) is .false. by default  :-)
      !   req (for required) is also .false. by default

      ! operators
*----------------------------------------------------------------------*
      call add_command_proto(DEF_GENERAL_OPERATOR,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_OP_FROM_OCC,tgt_info)
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/),def=.true.)
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'CORE',10,tgt_info,
     &     val_int=(/0,0,0,0,0,0,0,0,0,0/),def=.true.) ! works up to njoined=5
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'FORMAL',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'DESCR',0,tgt_info,
     &     val_str='',req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_SCALAR,tgt_info)
      call set_arg('_PROTO_',DEF_SCALAR,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_HAMILTONIAN,tgt_info)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/),def=.true.)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'FORMAL',1,tgt_info,
     &     val_int=(/10/),def=.true.)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'SET_X',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'X_SPCS',1,tgt_info,
     &     val_int=(/IEXTR/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_EXCITATION,tgt_info)
      call set_arg('_PROTO_',DEF_EXCITATION,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_EXCITATION,'MIN_RANK',1,tgt_info,
     &     val_int=(/1/),def=.true.)
      call set_arg('_PROTO_',DEF_EXCITATION,'MAX_RANK',0,tgt_info,
     &     val_int=(/-1/),req=.true.)
      call set_arg('_PROTO_',DEF_EXCITATION,'FORMAL',1,tgt_info,
     &     val_int=(/10/),def=.true.)
      call set_arg('_PROTO_',DEF_EXCITATION,'CHARGE',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_EXCITATION,'ADJOINT',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
*----------------------------------------------------------------------*
      ! deprecated
      call add_command_proto(DEF_DENSITY,tgt_info)
*----------------------------------------------------------------------*
      ! deprecated
      call add_command_proto(DEF_CC_HBAR_OP,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12GEMINAL,tgt_info)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'MIN_RANK',1,tgt_info,
     &     val_int=(/2/),def=.true.)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/),def=.true.)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'ANSATZ',1,tgt_info,
     &     val_int=(/3/),def=.true.)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'N_PART',1,tgt_info,
     &     val_int=(/0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12COEFF,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INT,tgt_info)
      call set_arg('_PROTO_',DEF_R12INT,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_R12INT,'MIN_RANK',1,tgt_info,
     &     val_int=(/2/),def=.true.)
      call set_arg('_PROTO_',DEF_R12INT,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/),def=.true.)
      call set_arg('_PROTO_',DEF_R12INT,'N_PART',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_R12INT,'CHARGE',1,tgt_info,
     &     val_int=(/0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INTERM,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(CLONE_OP,tgt_info)
      call set_arg('_PROTO_',CLONE_OP,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',CLONE_OP,'TEMPLATE',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SET_ORDER,tgt_info)
      call set_arg('_PROTO_',SET_ORDER,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SET_ORDER,'ORDER',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',SET_ORDER,'IDX_FREQ',0,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',SET_ORDER,'SPECIES',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SET_HERMIT,tgt_info)
      call set_arg('_PROTO_',SET_HERMIT,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SET_HERMIT,'CA_SYMMETRY',0,tgt_info,
     &     val_int=(/-1/),req=.true.)
! formulae:
*----------------------------------------------------------------------*
      call add_command_proto(CHECK_FORMGEN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CC_LAGRANGIAN,tgt_info)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'OP_H',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'OP_L',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'OP_T',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='CC Lagrange functional',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_ECC_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CCPT_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_MRCC_LAGRANGIAN,tgt_info)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'OPERATORS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,
     &                                          'MAXCOM_RES',0,tgt_info,
     &     val_int=(/-1/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'MAXCOM_EN',0,tgt_info,
     &     val_int=(/-1/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'MODE',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_MRCC_LAGRANGIAN,'TITLE',1,tgt_info,
     &     val_str='CC Lagrange functional',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_HHAT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CC_HBAR,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INTM_FORMAL,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INTM_CABS,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_MPR12_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CCR12_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CCR12_METRIC,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(SPLIT_R12EXC_FORMULA,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_EXP_FORMULA,tgt_info)
      call set_arg('_PROTO_',DEF_EXP_FORMULA,'SWITCH',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_EXP_FORMULA,'MODE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',DEF_EXP_FORMULA,'TITLE',1,tgt_info,
     &     val_str='My unnamed experimental formula',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_FORMULA,tgt_info)
      call set_arg('_PROTO_',DEF_FORMULA,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_FORMULA,'FORMULA',0,tgt_info,
     &     val_str=' ',req=.true.)
      call set_arg('_PROTO_',DEF_FORMULA,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'OPERATORS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'IDX_SV',1,tgt_info,
     &     val_int=(/-1/),req=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'BLK_MIN',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'BLK_MAX',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'CONNECT',0,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_CONNECT',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'AVOID',0,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'INPROJ',0,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_INPROJ',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'DESCR',0,tgt_info,
     &     (/'-'/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_DESCR',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.true./),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/1d0/),def=.true.)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'FIX_VTX',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(FACTOR_OUT,tgt_info)
      call set_arg('_PROTO_',FACTOR_OUT,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',FACTOR_OUT,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',FACTOR_OUT,'INTERM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',FACTOR_OUT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EXPAND,tgt_info)
      call set_arg('_PROTO_',EXPAND,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND,'INTERM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXPAND,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',EXPAND,'IMODE',1,tgt_info,
     &     val_int=(/0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(REPLACE,tgt_info)
      call set_arg('_PROTO_',REPLACE,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REPLACE,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REPLACE,'OP_LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REPLACE,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SUM_HERMIT,tgt_info)
      call set_arg('_PROTO_',SUM_HERMIT,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SUM_HERMIT,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SUM_HERMIT,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SUM_HERMIT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(REORDER_FORMULA,tgt_info)
      call set_arg('_PROTO_',REORDER_FORMULA,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REORDER_FORMULA,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REORDER_FORMULA,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(INVARIANT,tgt_info)
      call set_arg('_PROTO_',INVARIANT,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVARIANT,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVARIANT,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVARIANT,'OPERATORS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVARIANT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',INVARIANT,'REORDER',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DERIVATIVE,tgt_info)
      call set_arg('_PROTO_',DERIVATIVE,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DERIVATIVE,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DERIVATIVE,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DERIVATIVE,'OP_DERIV',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DERIVATIVE,'OP_MULT',1,tgt_info,
     &     val_label=(/' '/),def=.true.)
      call set_arg('_PROTO_',DERIVATIVE,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(LEQ_SPLIT,tgt_info)
      call set_arg('_PROTO_',LEQ_SPLIT,'LABEL_TRF',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'LABEL_RHS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'LABEL_RAW',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'OP_TRF',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'OP_RHS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'OP_X',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'TITLE_TRF',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',LEQ_SPLIT,'TITLE_RHS',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(OPTIMIZE,tgt_info)
      call set_arg('_PROTO_',OPTIMIZE,'LABEL_OPT',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',OPTIMIZE,'LABELS_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',OPTIMIZE,'INTERM',0,tgt_info,
     &     val_label=(/''/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_FORMULA,tgt_info)
      call set_arg('_PROTO_',PRINT_FORMULA,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',PRINT_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(TEX_FORMULA,tgt_info)
      call set_arg('_PROTO_',TEX_FORMULA,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',TEX_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(KEEP_TERMS,tgt_info)
      call set_arg('_PROTO_',KEEP_TERMS,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',KEEP_TERMS,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',KEEP_TERMS,'TERMS',1,tgt_info,
     &     val_int=(/-1/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_TERMS,tgt_info)
      call set_arg('_PROTO_',SELECT_TERMS,'OP_INCL',0,tgt_info,
     &     val_label=(/''/),def=.true.)
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_INCL',0,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',SELECT_TERMS,'OP_INCL_OR',0,tgt_info,
     &     val_label=(/''/),def=.true.)
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_INCL_OR',0,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',SELECT_TERMS,'OP_EXCL',0,tgt_info,
     &     val_label=(/''/),def=.true.)
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_EXCL',0,tgt_info,
     &     val_int=(/-1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_SPECIAL,tgt_info)
      call set_arg('_PROTO_',SELECT_SPECIAL,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_SPECIAL,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_SPECIAL,'OPERATORS',0,tgt_info,
     &     val_label=(/''/),def=.true.)
      call set_arg('_PROTO_',SELECT_SPECIAL,'TYPE',0,tgt_info,
     &     val_str='',req=.true.)
      call set_arg('_PROTO_',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEL_TERMS,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(MODIFY_FACTORIZATION,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_ORDER,tgt_info)
      call set_arg('_PROTO_',EXTRACT_ORDER,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_FREQ,tgt_info)
      call set_arg('_PROTO_',EXTRACT_FREQ,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(CLASS_FORMULA,tgt_info)
      call set_arg('_PROTO_',CLASS_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_HERMIT,tgt_info)
      call set_arg('_PROTO_',SELECT_HERMIT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_LINE,tgt_info)
      call set_arg('_PROTO_',SELECT_LINE,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'OP_INCL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'IGAST',0,tgt_info,
     &     val_int=(/1/),def=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='keep',def=.true.)
      call set_arg('_PROTO_',SELECT_LINE,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CUMULANTS,tgt_info)
      call set_arg('_PROTO_',DEF_CUMULANTS,'LABEL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CUMULANTS,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CUMULANTS,'OPERATORS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_CUMULANTS,'MODE',1,tgt_info,
     &     val_str='DENSITY',def=.true.)
      call set_arg('_PROTO_',DEF_CUMULANTS,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',DEF_CUMULANTS,'LEVEL',1,tgt_info,
     &     val_int=(/0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(INSERT,tgt_info)
      call set_arg('_PROTO_',INSERT,'LABEL_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INSERT,'LABEL_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INSERT,'OP_INS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INSERT,'OP_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INSERT,'OP_INCL',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INSERT,'TITLE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_MRCC_INTM,tgt_info)
      call set_arg('_PROTO_',DEF_MRCC_INTM,'FAC',1,tgt_info,
     &     val_rl8=(/0d0/),def=.true.)
      call set_arg('_PROTO_',DEF_MRCC_INTM,'MAXCOM',1,tgt_info,
     &     val_int=(/0/),def=.true.)
!     ME-lists:
*----------------------------------------------------------------------*
      call add_command_proto(DEF_ME_LIST,tgt_info)
      call set_arg('_PROTO_',DEF_ME_LIST,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'OPERATOR',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'IRREP',0,tgt_info,
     &     val_int=(/0/),req=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'2MS',0,tgt_info,
     &     val_int=(/0/),req=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'S2',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &     val_int=(/999/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',DEF_ME_LIST,'REC',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(RES_ME_LIST,tgt_info)
      call set_arg('_PROTO_',RES_ME_LIST,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(DELETE_ME_LIST,tgt_info)
      call set_arg('_PROTO_',DELETE_ME_LIST,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(ASSIGN_ME2OP,tgt_info)
      call set_arg('_PROTO_',ASSIGN_ME2OP,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',ASSIGN_ME2OP,'OPERATOR',0,tgt_info,
     &     val_label=(/''/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(IMPORT,tgt_info)
      call set_arg('_PROTO_',IMPORT,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',IMPORT,'TYPE',0,tgt_info,
     &     val_str='-',req=.true.)
      call set_arg('_PROTO_',IMPORT,'ENV',1,tgt_info,
     &     val_str=trim(env_type),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SPIN_PROJECT,tgt_info)
      call set_arg('_PROTO_',SPIN_PROJECT,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SPIN_PROJECT,'S2',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(PRECONDITIONER,tgt_info)
      call set_arg('_PROTO_',PRECONDITIONER,'LIST_PRC',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',PRECONDITIONER,'LIST_INP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',PRECONDITIONER,'MODE',1,tgt_info,
     &     val_str='dia-F',def=.true.)
      call set_arg('_PROTO_',PRECONDITIONER,'SHIFT',1,tgt_info,
     &     val_rl8=(/0d0/),def=.true.)
      call set_arg('_PROTO_',PRECONDITIONER,'THRES',1,tgt_info,
     &     val_rl8=(/-huge(1d1)/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(ADD,tgt_info)
      call set_arg('_PROTO_',ADD,'LIST_SUM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',ADD,'LISTS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',ADD,'FAC',0,tgt_info,
     &     val_rl8=(/1d0/),req=.true.)
      call set_arg('_PROTO_',ADD,'REPLACE',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SCALE,tgt_info)
      call set_arg('_PROTO_',SCALE,'LIST_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SCALE,'LIST_INP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SCALE,'FAC',0,tgt_info,
     &     val_rl8=(/0d0/),req=.true.)
      call set_arg('_PROTO_',SCALE,'LIST_SCAL',1,tgt_info,
     &     val_label=(/'-'/),def=.true.)
!      call set_arg('_PROTO_',SCALE,'NFAC',1,tgt_info,
!     &     val_int=(/-1/),def=.true.)
      call set_arg('_PROTO_',SCALE,'IDX_LIST',1,tgt_info,
     &     val_int=(/0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SCALE_COPY,tgt_info)
      call set_arg('_PROTO_',SCALE_COPY,'LIST_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SCALE_COPY,'LIST_INP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SCALE_COPY,'FAC',0,tgt_info,
     &     val_rl8=(/0d0/),req=.true.)
      call set_arg('_PROTO_',SCALE_COPY,'MODE',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',SCALE_COPY,'LIST_SHAPE',0,tgt_info,
     &     (/'-'/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(INVERT,tgt_info)
      call set_arg('_PROTO_',INVERT,'LIST_INV',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVERT,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',INVERT,'MODE',0,tgt_info,
     &     val_str='-',req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EVAL,tgt_info)
      call set_arg('_PROTO_',EVAL,'FORM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EVAL,'INIT',1,tgt_info,
     &     val_log=(/.true./),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EVALPROP,tgt_info)
      call set_arg('_PROTO_',EVALPROP,'DENS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EVALPROP,'RANK',0,tgt_info,
     &     val_int=(/1/),req=.true.)
      call set_arg('_PROTO_',EVALPROP,'ENV',1,tgt_info,
     &     val_str=env_type,def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SOLVENLEQ,tgt_info)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_OPT',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_RESID',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_PRC',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_E',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'FORM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'MODE',0,tgt_info,
     &     val_str='',req=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
      call set_arg('_PROTO_',SOLVENLEQ,'FORM_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SOLVELEQ,tgt_info)
      call set_arg('_PROTO_',SOLVELEQ,'LIST_OPT',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'LIST_PRC',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'OP_MVP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'OP_SVP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'OP_RHS',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'FORM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'MODE',0,tgt_info,
     &     val_str='',req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'N_ROOTS',0,tgt_info,
     &     val_int=(/0/),req=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'LIST_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
      call set_arg('_PROTO_',SOLVELEQ,'FORM_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SOLVEEVP,tgt_info)
      call set_arg('_PROTO_',SOLVEEVP,'LIST_OPT',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'LIST_PRC',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'OP_MVP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'OP_SVP',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'FORM',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'MODE',0,tgt_info,
     &     val_str='',req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'N_ROOTS',0,tgt_info,
     &     val_int=(/0/),req=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'LIST_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'FORM_SPC',0,tgt_info,
     &     (/'-'/),def=.true.)
      call set_arg('_PROTO_',SOLVEEVP,'TARG_ROOT',1,tgt_info,
     &     val_int=(/-1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(UNITY,tgt_info)
      call set_arg('_PROTO_',UNITY,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',UNITY,'FAC',1,tgt_info,
     &     val_rl8=(/1d0/),def=.true.)
      call set_arg('_PROTO_',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
      call set_arg('_PROTO_',UNITY,'MIN_BLK',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',UNITY,'MAX_BLK',1,tgt_info,
     &     val_int=(/0/),def=.true.)
      call set_arg('_PROTO_',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &     val_int=(/1/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SET_FREQ,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_RES,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_MEL,tgt_info)
      call set_arg('_PROTO_',PRINT_MEL,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',PRINT_MEL,'FORMAT',1,tgt_info,
     &     val_str='LIST',def=.true.)
      call set_arg('_PROTO_',PRINT_MEL,'COMMENT',1,tgt_info,
     &     val_str='---',def=.true.)
      call set_arg('_PROTO_',PRINT_MEL,'CHECK_THRESH',1,tgt_info,
     &     val_rl8=(/-1d0/),def=.true.)
      call set_arg('_PROTO_',PRINT_MEL,'EXPECTED',1,tgt_info,
     &     val_rl8=(/0d0/),def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(SET_MEL,tgt_info)
      call set_arg('_PROTO_',SET_MEL,'LIST',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',SET_MEL,'VAL_LIST',0,tgt_info,
     &     val_rl8=(/0d0/),req=.true.)
      call set_arg('_PROTO_',SET_MEL,'IDX_LIST',0,tgt_info,
     &     val_int=(/0/),req=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_DIAG,tgt_info)
      call set_arg('_PROTO_',EXTRACT_DIAG,'LIST_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXTRACT_DIAG,'LIST_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',EXTRACT_DIAG,'MODE',1,tgt_info,
     &     val_str='---',def=.true.)
*----------------------------------------------------------------------*
      call add_command_proto(REORDER_MEL,tgt_info)
      call set_arg('_PROTO_',REORDER_MEL,'LIST_RES',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REORDER_MEL,'LIST_IN',0,tgt_info,
     &     val_label=(/''/),req=.true.)
      call set_arg('_PROTO_',REORDER_MEL,'FROMTO',0,tgt_info,
     &     val_int=(/1/),req=.true.)
      call set_arg('_PROTO_',REORDER_MEL,'ADJOINT',1,tgt_info,
     &     val_log=(/.false./),def=.true.)
      call set_arg('_PROTO_',REORDER_MEL,'SEARCH',1,tgt_info,
     &     val_log=(/.false./),def=.true.)

      end
