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

      ! operators
*----------------------------------------------------------------------*
      call add_command_proto(DEF_GENERAL_OPERATOR,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_OP_FROM_OCC,tgt_info)
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'JOIN',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'CORE',10,tgt_info,
     &     val_int=(/0,0,0,0,0,0,0,0,0,0/)) ! works up to njoined=5
      call set_arg('_PROTO_',DEF_OP_FROM_OCC,'FORMAL',1,tgt_info,
     &     val_int=(/-1/))
*----------------------------------------------------------------------*
      call add_command_proto(DEF_SCALAR,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_HAMILTONIAN,tgt_info)
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'MIN_RANK',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'FORMAL',1,tgt_info,
     &     val_int=(/10/))
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'SET_X',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('_PROTO_',DEF_HAMILTONIAN,'X_SPCS',1,tgt_info,
     &     val_int=(/IEXTR/))
*----------------------------------------------------------------------*
      call add_command_proto(DEF_EXCITATION,tgt_info)
      call set_arg('_PROTO_',DEF_EXCITATION,'MIN_RANK',1,tgt_info,
     &     val_int=(/1/))
      call set_arg('_PROTO_',DEF_EXCITATION,'FORMAL',1,tgt_info,
     &     val_int=(/10/))
      call set_arg('_PROTO_',DEF_EXCITATION,'CHARGE',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_EXCITATION,'ADJOINT',1,tgt_info,
     &                                            val_log=(/.false./))
*----------------------------------------------------------------------*
      call add_command_proto(DEF_DENSITY,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CC_HBAR_OP,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12GEMINAL,tgt_info)
      call set_arg('_PROTO_',DEF_R12GEMINAL,'MIN_RANK',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('_PROTO_',DEF_R12GEMINAL,'MAX_RANK',1,tgt_info,
     &     val_int=(/2/))
      call set_arg('_PROTO_',DEF_R12GEMINAL,'ANSATZ',1,tgt_info,
     &     val_int=(/3/))
      call set_arg('_PROTO_',DEF_R12GEMINAL,'N_PART',1,tgt_info,
     &     val_int=(/0/))
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12COEFF,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_R12INTERM,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(CLONE_OP,tgt_info)
      call set_arg('_PROTO_',CLONE_OP,'ADJOINT',1,tgt_info,
     &     val_log=(/.false./))
*----------------------------------------------------------------------*
      call add_command_proto(SET_ORDER,tgt_info)
      call set_arg('_PROTO_',SET_ORDER,'ORDER',1,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',SET_ORDER,'IDX_FREQ',0,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',SET_ORDER,'SPECIES',1,tgt_info,
     &     val_int=(/-1/))
*----------------------------------------------------------------------*
      call add_command_proto(SET_HERMIT,tgt_info)
! formulae:
*----------------------------------------------------------------------*
      call add_command_proto(CHECK_FORMGEN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CC_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_ECC_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CCPT_LAGRANGIAN,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DEF_MRCC_LAGRANGIAN,tgt_info)
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
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_EXP_FORMULA,'MODE',1,tgt_info,
     &     val_str='---')
      call set_arg('_PROTO_',DEF_EXP_FORMULA,'TITLE',1,tgt_info,
     &     val_str='My unnamed experimental formula')
*----------------------------------------------------------------------*
      call add_command_proto(DEF_FORMULA,tgt_info)
      call set_arg('_PROTO_',DEF_FORMULA,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(EXPAND_OP_PRODUCT,tgt_info)
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'BLK_MIN',1,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'BLK_MAX',1,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'CONNECT',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_CONNECT',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'AVOID',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_AVOID',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'INPROJ',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'N_INPROJ',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'TITLE',1,tgt_info,
     &     val_str='---')
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'NEW',1,tgt_info,
     &     val_log=(/.true./))
      call set_arg('_PROTO_',EXPAND_OP_PRODUCT,'FAC',1,tgt_info,
     &     val_rl8=(/1d0/))
*----------------------------------------------------------------------*
      call add_command_proto(FACTOR_OUT,tgt_info)
      call set_arg('_PROTO_',FACTOR_OUT,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(EXPAND,tgt_info)
      call set_arg('_PROTO_',EXPAND,'TITLE',1,tgt_info,
     &     val_str='---')
      call set_arg('_PROTO_',EXPAND,'IMODE',1,tgt_info,
     &     val_int=(/0/))
*----------------------------------------------------------------------*
      call add_command_proto(REPLACE,tgt_info)
      call set_arg('_PROTO_',REPLACE,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(SUM_HERMIT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(INVARIANT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DERIVATIVE,tgt_info)
      call set_arg('_PROTO_',DERIVATIVE,'OP_MULT',1,tgt_info,
     &     val_label=(/' '/))
*----------------------------------------------------------------------*
      call add_command_proto(LEQ_SPLIT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(OPTIMIZE,tgt_info)
      call set_arg('_PROTO_',OPTIMIZE,'INTERM',0,tgt_info,
     &     val_label=(/''/))
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_FORMULA,tgt_info)
      call set_arg('_PROTO_',PRINT_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout')
*----------------------------------------------------------------------*
      call add_command_proto(TEX_FORMULA,tgt_info)
      call set_arg('_PROTO_',TEX_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout')
*----------------------------------------------------------------------*
      call add_command_proto(KEEP_TERMS,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_TERMS,tgt_info)
      call set_arg('_PROTO_',SELECT_TERMS,'OP_INCL',0,tgt_info,
     &     val_label=(/''/))
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_INCL',0,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',SELECT_TERMS,'OP_INCL_OR',0,tgt_info,
     &     val_label=(/''/))
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_INCL_OR',0,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',SELECT_TERMS,'OP_EXCL',0,tgt_info,
     &     val_label=(/''/))
      call set_arg('_PROTO_',SELECT_TERMS,'BLK_EXCL',0,tgt_info,
     &     val_int=(/-1/))
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_SPECIAL,tgt_info)
      call set_arg('_PROTO_',SELECT_SPECIAL,'MODE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(DEL_TERMS,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(MODIFY_FACTORIZATION,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_ORDER,tgt_info)
      call set_arg('_PROTO_',EXTRACT_ORDER,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_FREQ,tgt_info)
      call set_arg('_PROTO_',EXTRACT_FREQ,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(CLASS_FORMULA,tgt_info)
      call set_arg('_PROTO_',CLASS_FORMULA,'OUTPUT',1,tgt_info,
     &     val_str='stdout')
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_HERMIT,tgt_info)
      call set_arg('_PROTO_',SELECT_HERMIT,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(SELECT_LINE,tgt_info)
      call set_arg('_PROTO_',SELECT_LINE,'MODE',1,tgt_info,
     &     val_str='keep')
      call set_arg('_PROTO_',SELECT_LINE,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(DEF_CUMULANTS,tgt_info)
      call set_arg('_PROTO_',DEF_CUMULANTS,'MODE',1,tgt_info,
     &     val_str='DENSITY')
      call set_arg('_PROTO_',DEF_CUMULANTS,'TITLE',1,tgt_info,
     &     val_str='---')
      call set_arg('_PROTO_',DEF_CUMULANTS,'LEVEL',1,tgt_info,
     &     val_int=(/0/))
*----------------------------------------------------------------------*
      call add_command_proto(INSERT,tgt_info)
      call set_arg('_PROTO_',INSERT,'TITLE',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(DEF_MRCC_INTM,tgt_info)
!     ME-lists:
*----------------------------------------------------------------------*
      call add_command_proto(DEF_ME_LIST,tgt_info)
      call set_arg('_PROTO_',DEF_ME_LIST,'AB_SYM',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_ME_LIST,'CA_SYM',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_ME_LIST,'S2',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_ME_LIST,'MS_FIX',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_TYPE',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_IRREP',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',DEF_ME_LIST,'DIAG_MS',1,tgt_info,
     &     val_int=(/999/))
      call set_arg('_PROTO_',DEF_ME_LIST,'MIN_REC',1,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',DEF_ME_LIST,'MAX_REC',1,tgt_info,
     &     val_int=(/-1/))
*----------------------------------------------------------------------*
      call add_command_proto(RES_ME_LIST,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(DELETE_ME_LIST,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(ASSIGN_ME2OP,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(IMPORT,tgt_info)
      call set_arg('_PROTO_',IMPORT,'ENV',1,tgt_info,
     &     val_str=trim(env_type))
*----------------------------------------------------------------------*
      call add_command_proto(PRECONDITIONER,tgt_info)
      call set_arg('_PROTO_',PRECONDITIONER,'SHIFT',1,tgt_info,
     &     val_rl8=(/0d0/))
*----------------------------------------------------------------------*
      call add_command_proto(ADD,tgt_info)
      call set_arg('_PROTO_',ADD,'REPLACE',1,tgt_info,
     &     val_log=(/.false./))
*----------------------------------------------------------------------*
      call add_command_proto(SCALE,tgt_info)
      call set_arg('_PROTO_',SCALE,'LIST_SCAL',1,tgt_info,
     &     val_label=(/'-'/))
      call set_arg('_PROTO_',SCALE,'NFAC',1,tgt_info,
     &     val_int=(/-1/))
      call set_arg('_PROTO_',SCALE,'IDX_LIST',1,tgt_info,
     &     val_int=(/0/))
*----------------------------------------------------------------------*
      call add_command_proto(SCALE_COPY,tgt_info)
      call set_arg('_PROTO_',SCALE_COPY,'MODE',1,tgt_info,
     &     val_str='---')
      call set_arg('_PROTO_',SCALE_COPY,'LIST_SHAPE',0,tgt_info,
     &     val_label=(/'-'/))
*----------------------------------------------------------------------*
      call add_command_proto(INVERT,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(EVAL,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(EVALPROP,tgt_info)
      call set_arg('_PROTO_',EVALPROP,'ENV',1,tgt_info,
     &     val_str=env_type)
*----------------------------------------------------------------------*
      call add_command_proto(SOLVENLEQ,tgt_info)
      call set_arg('_PROTO_',SOLVENLEQ,'LIST_SPC',0,tgt_info,
     &     (/'-'/))
      call set_arg('_PROTO_',SOLVENLEQ,'FORM_SPC',0,tgt_info,
     &     (/'-'/))
*----------------------------------------------------------------------*
      call add_command_proto(SOLVELEQ,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(SOLVEEVP,tgt_info)
      call set_arg('_PROTO_',SOLVEEVP,'LIST_SPC',0,tgt_info,
     &     (/'-'/))
*----------------------------------------------------------------------*
      call add_command_proto(UNITY,tgt_info)
      call set_arg('_PROTO_',UNITY,'FAC',1,tgt_info,
     &     val_rl8=(/1d0/))
      call set_arg('_PROTO_',UNITY,'INIT',1,tgt_info,
     &     val_log=(/.false./))
      call set_arg('_PROTO_',UNITY,'MIN_BLK',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',UNITY,'MAX_BLK',1,tgt_info,
     &     val_int=(/0/))
      call set_arg('_PROTO_',UNITY,'MS_SYM_SIGN',1,tgt_info,
     &     val_int=(/1/))
*----------------------------------------------------------------------*
      call add_command_proto(SET_FREQ,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_RES,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(PRINT_MEL,tgt_info)
      call set_arg('_PROTO_',PRINT_MEL,'FORMAT',1,tgt_info,
     &     val_str='LIST')
      call set_arg('_PROTO_',PRINT_MEL,'COMMENT',1,tgt_info,
     &     val_str='---')
*----------------------------------------------------------------------*
      call add_command_proto(SET_MEL,tgt_info)
*----------------------------------------------------------------------*
      call add_command_proto(EXTRACT_DIAG,tgt_info)
      call set_arg('_PROTO_',EXTRACT_DIAG,'EXTEND',1,tgt_info,
     &     val_log=(/.false./))
*----------------------------------------------------------------------*
      call add_command_proto(REORDER_MEL,tgt_info)

      end
