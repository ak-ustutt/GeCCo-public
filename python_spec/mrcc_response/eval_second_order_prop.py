#This is python interface to evaluate the second order properties using the all the zeroeth 
#order paramaters and the first order cluster amplitudes and ci coefficients.

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

#Getting the option from response data
_option=_response_data['option']

#preprocessing the second order response function. First changing the GeCCo specific Lagrangian 
#to the actual Lagrangian of ic-MRCC and then replacing the electronic Hamiltonian with the Perturbation
#NOTE: This perturbation is still the previous one. It should be a new one for a non-diagonal component of the properties

new_target('F_preRSPNS(2)')

depend('F_MRCC_LAG')

DEF_SCALAR({LABEL:'preRSPNS(2)'})
DEF_SCALAR({LABEL:'preRSPNS(2)_1'})
DEF_SCALAR({LABEL:'preRSPNS(2)_2'})
DEF_SCALAR({LABEL:'preRSPNS(2)_3'})

DERIVATIVE({LABEL_RES:'F_preRSPNS(2)_1',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'preRSPNS(2)_1',
            OP_DERIV:'L',
            OP_MULT:'L'})

INVARIANT({LABEL_RES:'F_TEMP(2)',
           LABEL_IN:'F_MRCC_LAG',
           OP_RES:'preRSPNS(2)_2',
           OPERATORS:'L'})

DERIVATIVE({LABEL_RES:'F_preRSPNS(2)_2',
            LABEL_IN:'F_TEMP(2)',
            OP_RES:'preRSPNS(2)_2',
            OP_DERIV:'C0^+',
            OP_MULT:'C0_bar'})

DEF_FORMULA({LABEL:'F_preRSPNS(2)',
             FORMULA:'preRSPNS(2)=preRSPNS(2)_1+preRSPNS(2)_2'})

EXPAND({LABEL_RES:'F_preRSPNS(2)',
        LABEL_IN:'F_preRSPNS(2)',
        INTERM:['F_preRSPNS(2)_1','F_preRSPNS(2)_2']})

REPLACE({LABEL_RES:'F_preRSPNS(2)_3',
         LABEL_IN:'F_preRSPNS(2)',
         OP_LIST:['H','V(1)']})

INVARIANT({LABEL_RES:'F_preRSPNS(2)_3',
           LABEL_IN:'F_preRSPNS(2)_3',
           OP_RES:'preRSPNS(2)_3',
           OPERATORS:'H'})


for i in range(0,n_par):

    i_par = str(i+1)

#Getting the part of the second order response function which is a product of perturbation and 
#first order wave function paramaeters, T(1) and C(1). It Corresponds the scalar RSPNS(2)_1

    new_target('F_RSPNS(2)'+i_par)

    depend('F_preRSPNS(2)')

    DEF_SCALAR({LABEL:'RSPNS(2)'+i_par})

    _deriv_arg = {}
    _deriv_arg[LABEL_RES] = 'F_RSPNS(2)'+i_par
    _deriv_arg[LABEL_IN]  = 'F_preRSPNS(2)_3'
    _deriv_arg[OP_RES]    = 'RSPNS(2)'+i_par

    _replace_arg = {}
    _replace_arg[LABEL_RES] = 'F_RSPNS(2)'+i_par
    _replace_arg[LABEL_IN]  = 'F_RSPNS(2)'+i_par
    _replace_arg[OP_LIST]   = ['C0(1)_x'+i_par,'C0(1)'+i_par+'^+']

    if _option == 1: 
# option=1 has extra contribution from the 'C0(1)^+'
        _deriv_arg[OP_DERIV] = ['T','C0','C0^+']
        _deriv_arg[OP_MULT]  = ['T(1)'+i_par,'C0(1)'+i_par,'C0(1)_x'+i_par]

        DERIVATIVE(_deriv_arg)

        REPLACE(_replace_arg)

    elif _option == 2: 
        _deriv_arg[OP_DERIV] = ['T','C0']
        _deriv_arg[OP_MULT]  = ['T(1)'+i_par,'C0(1)'+i_par]

        DERIVATIVE(_deriv_arg)

    else:
        quit_error('Input error: unknown option for ic-MRCC properties') 

    
    #PRINT_FORMULA({LABEL:'F_preRSPNS(2)_3'})
    #PRINT_FORMULA({LABEL:'F_RSPNS(2)'+i_par})

#Getting the part of the second order response function which is the product of two first order wave function 
#parameters, T(1) and C0(1). This also includes the cross terms between T(1) and C0(1). It corresponds the scalar RSPNS(2)_2.
#npar=2 would involve the product of T(1) and C0(1) corresponding to different frequencies. 

new_target('F_RSPNS(2)')

depend('F_RSPNS(2)1')

if (n_par == 2):
    depend('F_RSPNS(2)2')

DEF_SCALAR({LABEL:'RSPNS(2)_coup_1'})
DEF_SCALAR({LABEL:'RSPNS(2)_coup_2'})
DEF_SCALAR({LABEL:'RSPNS(2)_coup_3'})
DEF_SCALAR({LABEL:'RSPNS(2)_coup_4'})


if _option == 1:

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup_3',
                LABEL_IN:'F_preRSPNS(2)_1',
                OP_RES:'RSPNS(2)_coup_3',
                OP_DERIV:['T','C0'],
                OP_MULT:['T(1)1','C0(1)1']})

    REPLACE({LABEL_RES:'F_RSPNS(2)_coup_3',
             LABEL_IN:'F_RSPNS(2)_coup_3',
             OP_LIST:['C0^+','C0(1)1'+'^+']})

    if n_par == 2:
 
        DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup_4',
                    LABEL_IN:'F_preRSPNS(2)_1',
                    OP_RES:'RSPNS(2)_coup_4',
                    OP_DERIV:['T','C0'],
                    OP_MULT:['T(1)2','C0(1)2']})

        REPLACE({LABEL_RES:'F_RSPNS(2)_coup_4',
                 LABEL_IN:'F_RSPNS(2)_coup_4',
                 OP_LIST:['C0^+','C0(1)2'+'^+']})

_deriv_arg = {}
_deriv_arg[LABEL_RES] = 'F_intRSPNS(2)1' 
_deriv_arg[LABEL_IN]  = 'F_preRSPNS(2)' 
_deriv_arg[OP_RES]    = 'RSPNS(2)_coup_1'
_deriv_arg[OP_DERIV]  =  ['T','C0']

if n_par == 1:
    _deriv_arg[OP_MULT] = ['T(1)1','C0(1)1']
else:
    _deriv_arg[OP_MULT] = ['T(1)2','C0(1)2']


DERIVATIVE(_deriv_arg)

DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup_1',
            LABEL_IN:'F_intRSPNS(2)1',
            OP_RES:'RSPNS(2)_coup_1',
            OP_DERIV:'T',
            OP_MULT:'T(1)1'})
    
if n_par == 2:

    DERIVATIVE({LABEL_RES:'F_intRSPNS(2)2',
                LABEL_IN:'F_preRSPNS(2)',
                OP_RES:'RSPNS(2)_coup_2',
                OP_DERIV:'C0',
                OP_MULT:'C0(1)1'})
    
    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_coup_2',
                LABEL_IN:'F_intRSPNS(2)2',
                OP_RES:'RSPNS(2)_coup_2',
                OP_DERIV:'T',
                OP_MULT:'T(1)2'})
    
#PRINT_FORMULA({LABEL:'F_RSPNS(2)_coup_1'})
#PRINT_FORMULA({LABEL:'F_RSPNS(2)_coup_2'})

DEF_SCALAR({LABEL:'RSPNS(2)'})
DEF_SCALAR({LABEL:'RSPNS(2)_1'})
DEF_SCALAR({LABEL:'RSPNS(2)_2'})


if n_par ==1:
    DEF_FORMULA({LABEL:'F_RSPNS(2)_1',
                 FORMULA:'RSPNS(2)_1=RSPNS(2)1+RSPNS(2)1'})
    EXPAND({LABEL_RES:'F_RSPNS(2)_1',
            LABEL_IN:'F_RSPNS(2)_1',
            INTERM:['F_RSPNS(2)1']})

else:
    DEF_FORMULA({LABEL:'F_RSPNS(2)_1',
                 FORMULA:'RSPNS(2)_1=RSPNS(2)1+RSPNS(2)2'})
    EXPAND({LABEL_RES:'F_RSPNS(2)_1',
            LABEL_IN:'F_RSPNS(2)_1',
            INTERM:['F_RSPNS(2)1','F_RSPNS(2)2']})

_def_form_arg = {}
_def_form_arg[LABEL] = 'F_RSPNS(2)_2'

_expand_form = {}
_expand_form[LABEL_RES] = 'F_RSPNS(2)_2'
_expand_form[LABEL_IN] = 'F_RSPNS(2)_2'

#Here we are adding different parts of the RSPNS(2)_2 together. This would depend on 
#both of _option and n_par. 
if _option ==1:
    if n_par == 1:
        _def_form_arg[FORMULA] = 'RSPNS(2)_2=RSPNS(2)_coup_1+RSPNS(2)_coup_3'
        _expand_form[INTERM] = ['F_RSPNS(2)_coup_1','F_RSPNS(2)_coup_3']
    if n_par == 2:
        _def_form_arg[FORMULA] = 'RSPNS(2)_2=RSPNS(2)_coup_1+RSPNS(2)_coup_2+RSPNS(2)_coup_3+RSPNS(2)_coup_4' #Should have a factor half adjoint
        _expand_form[INTERM] = ['F_RSPNS(2)_coup_1','F_RSPNS(2)_coup_2','F_RSPNS(2)_coup_3','F_RSPNS(2)_coup_4']
if _option ==2:
    if n_par == 1:
        _def_form_arg[FORMULA] = 'RSPNS(2)_2=RSPNS(2)_coup_1'
        _expand_form[INTERM] = ['F_RSPNS(2)_coup_1']
    if n_par == 2:
        _def_form_arg[FORMULA] = 'RSPNS(2)_2=RSPNS(2)_coup_1+RSPNS(2)_coup_2'
        _expand_form[INTERM] = ['F_RSPNS(2)_coup_1','F_RSPNS(2)_coup_2']

#PRINT_FORMULA({LABEL:'F_RSPNS(2)_2'})

DEF_FORMULA(_def_form_arg)

EXPAND(_expand_form)

#RSPNS(2)_1 and RSPNS(2)_2 are then summed up to calculate the total response function RSPNS(2)
DEF_FORMULA({LABEL:'F_RSPNS(2)',
             FORMULA:'RSPNS(2)=RSPNS(2)_1+RSPNS(2)_2'})

EXPAND({LABEL_RES:'F_RSPNS(2)',
        LABEL_IN:'F_RSPNS(2)',
        INTERM:'F_RSPNS(2)_coup_1'})

new_target('OPT_RSPNS(2)')

depend('F_RSPNS(2)')

#Defining the me-list for all the response function
_op_list={'RSPNS(2)':'ME_RSPNS(2)',
          'RSPNS(2)_1':'ME_RSPNS(2)_1',
          'RSPNS(2)_2':'ME_RSPNS(2)_2'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

#Optimizing all the formula
OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)_1',
          LABELS_IN:'F_RSPNS(2)_1'})

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)_2',
          LABELS_IN:'F_RSPNS(2)_2'})

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)',
          LABELS_IN:'F_RSPNS(2)'})

#evaluating the response function:
if (_restart<3):
    new_target('EVAL_RSPNS(2)',True)
else:
    new_target('EVAL_RSPNS(2)')

depend('OPT_RSPNS(2)')

EVALUATE({FORM:'FOPT_RSPNS(2)_1'})
PRINT_MEL({LIST:'ME_RSPNS(2)_1',
           FORMAT:'SCAL',
           COMMENT:'Asymmetric part of the second order property: '})

EVALUATE({FORM:'FOPT_RSPNS(2)_2'})
PRINT_MEL({LIST:'ME_RSPNS(2)_2',
           FORMAT:'SCAL',
           COMMENT:'Symmetric part of the second order property: '})

EVALUATE({FORM:'FOPT_RSPNS(2)'})
PRINT_MEL({LIST:'ME_RSPNS(2)',
           FORMAT:'SCAL',
           COMMENT:'Total value of the second order property: '})
