#This is python interface to evaluate the second order properties using the all the zeroeth 
#and  first order parameters. This code only work when the ci coeffients are not relaxed at all.

from gecco_interface import *
from get_response_data import _response_data

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

for i in range(0,n_par):

    i_par = str(i+1)

    new_target('ORTH_C0_bar(1)'+i_par)

    depend('DEF_ME_C0_bar(1)'+i_par)

    _op_list={'C0_bar(1)_orth_'+i_par:'C0_bar'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op]})

    EXPAND_OP_PRODUCT({LABEL:'F_ORTH_C0_bar(1)'+i_par,
                       OP_RES:'C0_bar(1)_orth_'+i_par,
                       NEW:True,
                       OPERATORS:['C0_bar(1)_orth_'+i_par,'C0(1)'+i_par,'C0_bar(1)_orth_'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_ORTH_C0_bar(1)'+i_par,
                       OP_RES:'C0_bar(1)_orth_'+i_par,
                       NEW:False,
                       OPERATORS:['C0_bar(1)_orth_'+i_par,'C0_bar(1)'+i_par,'C0_bar^+','C0_bar','C0_bar(1)_orth_'+i_par],
                       IDX_SV:[1,2,3,4,1],
                       AVOID:[1,3,2,5],
                       FAC:-1.0})

    _op_list={'C0_bar(1)_orth_'+i_par:'ME_C0_bar(1)_orth_'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    OPTIMIZE({LABEL_OPT:'FOPT_ORTH_C0_bar(1)'+i_par,
              LABELS_IN:['F_ORTH_C0_bar(1)'+i_par]})

    EVALUATE({FORM:'FOPT_ORTH_C0_bar(1)'+i_par})

    PRINT_MEL({LIST:'ME_C0_bar(1)_orth_'+i_par})
#   SCALE_COPY({LIST_RES:'ME_C0_bar(1)'+i_par,LIST_INP:'ME_C0_bar(1)_orth_'+i_par,FAC:1.0})

new_target('F_RSPNS(3)')

depend('F_preRSPNS(2)')
depend('IMPORT_V(1)')

for i in range(0,n_par):

    i_par = str(i+1)

    depend('DEF_ME_L(1)'+i_par)
    depend('DEF_ME_C0_bar(1)'+i_par)
    depend('ORTH_C0_bar(1)'+i_par)

    CLONE_OPERATOR({LABEL:'C0_bar(1)_x'+i_par,
                    TEMPLATE:'C0_bar'})

    DEF_SCALAR({LABEL:'RSPNS(3)_1'+i_par})
    DEF_SCALAR({LABEL:'RSPNS(3)_2'+i_par})

    _deriv_arg={}
    _deriv_arg[LABEL_RES] = 'F_RSPNS(3)_1'+i_par
    _deriv_arg[LABEL_IN] = 'F_preRSPNS(2)' 
    _deriv_arg[OP_RES] = 'RSPNS(3)_1'+i_par
    _deriv_arg[OP_DERIV] =['T','C0','C0^+']
    _deriv_arg[OP_MULT]  = ['T(1)'+i_par,'C0(1)'+i_par,'C0_bar(1)_x'+i_par]

    _replace_arg = {}
    _replace_arg[LABEL_RES] = 'F_RSPNS(3)_1'+i_par
    _replace_arg[LABEL_IN]  = 'F_RSPNS(3)_1'+i_par
    _replace_arg[OP_LIST]   = ['C0_bar(1)_x'+i_par,'C0(1)'+i_par+'^+']

    DERIVATIVE(_deriv_arg)

    REPLACE(_replace_arg)

    _deriv_arg={}
    _deriv_arg[LABEL_RES] = 'F_RSPNS(3)_2'+i_par
    _deriv_arg[LABEL_IN] = 'F_preRSPNS(2)' 
    _deriv_arg[OP_RES] = 'RSPNS(3)_2'+i_par
    _deriv_arg[OP_DERIV] =['L','C0_bar']
    _deriv_arg[OP_MULT]  = ['L(1)'+i_par,'C0_bar(1)'+i_par]

    DERIVATIVE(_deriv_arg)

DEF_SCALAR({LABEL:'RSPNS(3)'})

if (n_par == 1):

    DEF_FORMULA({LABEL:'F_RSPNS(3)',
                 FORMULA:'RSPNS(3)=RSPNS(3)_11+RSPNS(3)_21'})

    EXPAND({LABEL_RES:'F_RSPNS(3)',
            LABEL_IN:'F_RSPNS(3)',
            INTERM:['F_RSPNS(3)_11','F_RSPNS(3)_21']})

else:
    print 'Not done yet:'

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

PRINT_MEL({LIST:'ME_C0_bar(1)'+i_par})

EVALUATE({FORM:'FOPT_RSPNS(3)'})
PRINT_MEL({LIST:'ME_RSPNS(3)',
           FORMAT:'SCAL',
           COMMENT:'Total value of the second order property: '})
