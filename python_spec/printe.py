from gecco_interface import *

new_target('PRINTE',True)

depend('SOLVE_MRCC')

FACTOR_OUT({LABEL_RES:'F_E(RDMs)',
            LABEL_IN:'F_MRCC_E',
            INTERM:'F_DENS0'})

OPTIMIZE({LABEL_OPT:'FOPT_Erdm',
          LABELS_IN:'F_E(RDMs)'})

RESET_ME_LIST({LIST:'ME_E(MR)'})

EVALUATE({FORM:'FOPT_Erdm'})

PRINT_MEL({LIST:'ME_E(MR)',
           FORMAT:'SCAL F20.12',
           COMMENT:'RDM-based final energy:'})

export_targets()
