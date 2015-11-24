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
        LIST:'X_TRM_LIST_DAG',
        OPERATOR:'X_TRM'})


debug_FORM('FORM_T2_orth')


SOLVE_NLEQ({
        LIST_OPT:'T2_ca_LST',
        LIST_RESID:'Oges_LST',
        LIST_PRC:'PRECON_LST',
        LIST_E:'PT_LAG_LST',
        FORM:'FOPT_PT_LAG',
        MODE:'TRF',
        FORM_SPC:['FOPT_T2_orth'],
        LIST_SPC:['T2_orth_LIST','X_TRM_LIST','X_TRM_LIST_DAG']
        })

