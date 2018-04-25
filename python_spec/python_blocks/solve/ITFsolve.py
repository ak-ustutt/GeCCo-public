"""Solves the icMRCCSD equations used in ITF translator

History:

Based on icMRCCSDsolve.py

"""
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


# Don't care about optref for now...
optref_n = int(keywords.get("calculate.solve.non_linear.optref"))
optref = False
if keywords.get("calculate.solve.non_linear.optref") is not None:
    if optref_n != 0:
        quit_error('calculate solve non_linear optref should be 0.')


new_target('SOLVE_ITF')
heading('Solving the icMRCC equations for ITF')

depend('DEF_FORM_ITF_LAG')
depend('EVAL_E0')
depend('BUILD_PRECON')

ASSIGN_ME2OP({
        LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM'})

SOLVE_NLEQ({
        LIST_OPT:['ME_T2g'],
        LIST_RESID:['ME_O2g'],
        LIST_PRC:['ME_PRECON2g'],
        LIST_E:'MRCC_LAG_LST',
        FORM:'FOPT_MRCC_LAG',
        MODE:'TRF',
        FORM_SPC:['FOPT_T2_orth'],
        LIST_SPC:['ME_T2_orth','ME_X_TRM','ME_X_TRM_DAG']
        })

PUSH_RESULT({LIST:'MRCC_LAG_LST',COMMENT:"MRCC", FORMAT:"SCAL F20.14"})
