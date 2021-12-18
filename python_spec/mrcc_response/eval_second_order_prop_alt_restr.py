#This is python interface to evaluate the second order properties using the all the zeroeth 
#and  first order parameters. This code only work when the ci coeffients are not relaxed at all.

from gecco_interface import *
from python_spec.mrcc_response.get_response_data import _response_data

#Getting the frequency
_freq=_response_data['freq']
#Getting the value of the restart option
_restart=_response_data['restart']

#n_par tells how many version of the same operator has to be defined 
#depending on whether we are doing static or dynamic property calcualtion
if (_freq == 0.0):
    n_par = 1
else:
    n_par = 2

new_target('F_RSPNS(3)')

depend('F_preRSPNS(2)')

for i in range(0,n_par):

    i_par = str(i+1)

#   depend('SOLVE_T(1)'+i_par)

    DEF_SCALAR({LABEL:'RSPNS(3)_1'+i_par})
    DEF_SCALAR({LABEL:'RSPNS(3)_2'+i_par})

    DERIVATIVE({LABEL_RES:'F_RSPNS(3)_1'+i_par,
                LABEL_IN:'F_preRSPNS(2)_1',
                OP_RES:'RSPNS(3)_1'+i_par,
                OP_DERIV:'T',
                OP_MULT:'T(1)'+i_par})

    DERIVATIVE({LABEL_RES:'F_RSPNS(3)_2'+i_par,
                LABEL_IN:'F_preRSPNS(2)_1',
                OP_RES:'RSPNS(3)_2'+i_par,
                OP_DERIV:'L',
                OP_MULT:'L(1)'+i_par})

DEF_SCALAR({LABEL:'RSPNS(3)'})

if (n_par == 1):

    DEF_FORMULA({LABEL:'F_RSPNS(3)',
                 FORMULA:'RSPNS(3)=RSPNS(3)_11+RSPNS(3)_21'})

    EXPAND({LABEL_RES:'F_RSPNS(3)',
            LABEL_IN:'F_RSPNS(3)',
            INTERM:['F_RSPNS(3)_11','F_RSPNS(3)_21']})

else:
    print('Not done yet:')

new_target('OPT_RSPNS(2)_alt')

depend('F_RSPNS(3)')

_op_list={'RSPNS(3)':'ME_RSPNS(3)'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(3)',
          LABELS_IN:'F_RSPNS(3)'})


PRINT_FORMULA({LABEL:'F_RSPNS(3)'})

new_target('EVAL_RSPNS(2)_alt',True)

depend('OPT_RSPNS(2)_alt')

EVALUATE({FORM:'FOPT_RSPNS(3)'})
PRINT_MEL({LIST:'ME_RSPNS(3)',
           FORMAT:'SCAL',
           COMMENT:'Total value of the second order property: '})
