from gecco_interface import *
from gecco_modules.NoticeUtil import *

#------------------------------------------------------------------
#------------------------------------------------------------------
#SOLVE it
#------------------------------------------------------------------
#------------------------------------------------------------------

new_target('SOLVE_MRCCPT2')
depend('DEF_FORM_PT_LAG')
depend('BUILD_PRECON')

debug_MEL('PRECON_LST')

#SOLVE_NLEQ uses FOPT_T2_orth for both transformations and needs X_TRM_LIST_DAG also bound to that opperator


ASSIGN_ME2OP({
        LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM'})

debug_FORM('FORM_T2_orth')


SOLVE_NLEQ({
        LIST_OPT:'ME_T2g',
        LIST_RESID:'ME_O2g',
        LIST_PRC:'ME_PRECON2g',
        LIST_E:'PT_LAG_LST',
        FORM:'FOPT_PT_LAG',
        MODE:'TRF',
        FORM_SPC:['FOPT_T2_orth'],
        LIST_SPC:['ME_T2_orth','ME_X_TRM','ME_X_TRM_DAG']
        })

PUSH_RESULT({LIST:'PT_LAG_LST',COMMENT:"MRCCPT2", FORMAT:"SCAL F20.14"})
