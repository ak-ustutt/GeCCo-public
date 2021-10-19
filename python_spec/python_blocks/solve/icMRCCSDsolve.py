"""Solves the icMRCCSD equations

Implementation based on separated ME-lists and operators for T1 and T2

History:

Yuri August 2017: Creation based on MRCC2solve.py.

"""
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


optref_n = int(keywords.get("calculate.solve.non_linear.optref"))
optref = False
if keywords.get("calculate.solve.non_linear.optref") is not None:
    if optref_n == -3:
        optref = True
    elif optref_n != 0:
        quit_error('calculate solve non_linear optref should be 0 or -3.')

new_target('SOLVE_MRCC')
heading('Solving the icMRCC equations')

depend('DEF_FORM_MRCC_LAG')
depend('EVAL_E0')
depend('BUILD_PRECON')

ASSIGN_ME2OP({
        LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM'})

if (True): # switch if for debug purposes only solving for T2g component is requested set to False
    if optref:
        depend("FOPT_OMG_C0")
        _form_spc = ['FOPT_T1_orth','FOPT_GAM_S','FOPT_T2_orth']
        _list_spc = ['ME_T1_orth','ME_X_TRM','ME_X_TRM_DAG',"ME_P_PROJ","ME_GAM_S","ME_GAM_S_ISQ",'ME_T2_orth']

    else:
        _form_spc = ['FOPT_T1_orth','FOPT_T2_orth']
        _list_spc = ['ME_T1_orth','ME_X_TRM','ME_X_TRM_DAG','ME_T2_orth']

    SOLVE_NLEQ({
            LIST_OPT:['ME_T1','ME_T2g'],
            LIST_RESID:['ME_O1','ME_O2g'],
            LIST_PRC:['ME_PRECON1','ME_PRECON2g'],
            LIST_E:'MRCC_LAG_LST',
            FORM:'FOPT_MRCC_LAG',
            MODE:'TRF TR0',
            #MODE:'DIA DIA',
            FORM_SPC:_form_spc,
            LIST_SPC:_list_spc
            })
else:
    if optref:
        depend("FOPT_OMG_C0")
        _form_spc = ['FOPT_T2_orth','FOPT_GAM_S']
        _list_spc = ['ME_T2_orth','ME_X_TRM','ME_X_TRM_DAG',"ME_P_PROJ","ME_GAM_S","ME_GAM_S_ISQ"]

    else:
        _form_spc = ['FOPT_T2_orth']
        _list_spc = ['ME_T2_orth','ME_X_TRM','ME_X_TRM_DAG']

    SOLVE_NLEQ({
            LIST_OPT:['ME_T2g'],
            LIST_RESID:['ME_O2g'],
            LIST_PRC:['ME_PRECON2g'],
            LIST_E:'MRCC_LAG_LST',
            FORM:'FOPT_MRCC_LAG',
            MODE:'TRF',
            FORM_SPC:_form_spc,
            LIST_SPC:_list_spc
            })
    

PUSH_RESULT({LIST:'MRCC_LAG_LST',COMMENT:"MRCC", FORMAT:"SCAL F20.14"})
