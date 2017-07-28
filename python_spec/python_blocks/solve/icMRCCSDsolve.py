from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

new_target('SOLVE_MRCC')

depend('DEF_FORM_MRCC_LAG')
depend('EVAL_E0')
depend('BUILD_PRECON')


SOLVE_NLEQ({
        LIST_OPT:['ME_T1','ME_T2g'],
        LIST_RESID:['ME_O1','ME_O2g'],
        LIST_PRC:['ME_PRECON1','ME_PRECON2g'],
        LIST_E:'MRCC_LAG_LST',
        FORM:'FOPT_MRCC_LAG2',
        MODE:'TRF TR0',
        FORM_SPC:['FOPT_T1_orth','FOPT_T2_orth'],
        LIST_SPC:['ME_T1_orth','ME_X_TRM','ME_X_TRM_DAG','ME_T2_orth']
        })


PUSH_RESULT({LIST:'MRCC_LAG_LST',COMMENT:"MRCC", FORMAT:"SCAL F20.14"})
