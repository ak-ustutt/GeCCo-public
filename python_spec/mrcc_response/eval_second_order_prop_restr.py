#This is python interface to evaluate the second order properties using the all the zeroeth 
#order paramaters and the first order cluster amplitudes. This code only work when the ci 
#coeffients are not relaxed at all.

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


#Getting the part of the Lagrangian after replacing the electronic hamiltonian with the perturbation.
#NOTE: This perturbation is still the previous one. It should be a new one for a non-diagonal component of the properties

new_target('F_preRSPNS(2)')

depend('F_MRCC_LAG')

DEF_SCALAR({LABEL:'preRSPNS(2)_1'})

REPLACE({LABEL_RES:'F_preRSPNS(2)_1',
         LABEL_IN:'F_MRCC_LAG',
         OP_LIST:['H','V(1)']})

INVARIANT({LABEL_RES:'F_preRSPNS(2)_1',
           LABEL_IN:'F_preRSPNS(2)_1',
           OP_RES:'preRSPNS(2)_1',
           OPERATORS:'H'})

#Getting two different parts of the second order Response function after replacing one of the zeroeth order T with the first order T
for i in range(0,n_par):

    i_par = str(i+1)

    new_target('F_RSPNS(2)'+i_par)

    depend('SOLVE_T(1)'+i_par)
    depend('F_preRSPNS(2)')

    DEF_SCALAR({LABEL:'RSPNS(2)'+i_par})

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)'+i_par,
                LABEL_IN:'F_preRSPNS(2)_1',
                OP_RES:'RSPNS(2)'+i_par,
                OP_DERIV:'T',
                OP_MULT:'T(1)'+i_par})

#Getting the part of the Response function which is the product of two first order T
new_target('F_RSPNS(2)')

depend('F_RSPNS(2)1')

if (n_par == 2):
    depend('F_RSPNS(2)2')

DEF_SCALAR({LABEL:'RSPNS(2)_coup'})

DERIVATIVE({LABEL_RES:'F_intRSPNS(2)',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'RSPNS(2)_coup',
            OP_DERIV:'T',
            OP_MULT:'T(1)1'})

if (n_par == 1):

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup',
                LABEL_IN:'F_intRSPNS(2)',
                OP_RES:'RSPNS(2)_coup',
                OP_DERIV:'T',
                OP_MULT:'T(1)1'})
else:

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup',
                LABEL_IN:'F_intRSPNS(2)',
                OP_RES:'RSPNS(2)_coup',
                OP_DERIV:'T',
                OP_MULT:'T(1)2'})
    
DEF_SCALAR({LABEL:'RSPNS(2)'})

#Adding different parts of the second order response function:

if (n_par == 1):

    DEF_FORMULA({LABEL:'F_RSPNS(2)',
                 FORMULA:'RSPNS(2)=RSPNS(2)1+RSPNS(2)1+RSPNS(2)_coup'})

    EXPAND({LABEL_RES:'F_RSPNS(2)',
            LABEL_IN:'F_RSPNS(2)',
            INTERM:['F_RSPNS(2)1','F_RSPNS(2)_coup']})

else:

    DEF_FORMULA({LABEL:'F_RSPNS(2)',
                 FORMULA:'RSPNS(2)=RSPNS(2)1+RSPNS(2)2+RSPNS(2)_coup'})

    EXPAND({LABEL_RES:'F_RSPNS(2)',
            LABEL_IN:'F_RSPNS(2)',
            INTERM:['F_RSPNS(2)1','F_RSPNS(2)2','F_RSPNS(2)_coup']})

new_target('OPT_RSPNS(2)')

depend('F_RSPNS(2)')
depend('DEF_ME_L')

#Getting the me-list for RSPNS(2)
_op_list={'RSPNS(2)':'ME_RSPNS(2)'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)',
          LABELS_IN:'F_RSPNS(2)'})

#evaluating the response function:
if (_restart<3):
    new_target('EVAL_RSPNS(2)',True)
else:
    new_target('EVAL_RSPNS(2)')


depend('OPT_RSPNS(2)')

EVALUATE({FORM:'FOPT_RSPNS(2)'})

PRINT_MEL({LIST:'ME_RSPNS(2)'})
