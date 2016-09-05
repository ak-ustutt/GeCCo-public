#This is a python interface for GeCCo to solve the first order T(\omega) to 
#calculate static and dynamic polarizability using the linear response method.
#This particular interface started to be written on September, 2015

from gecco_interface import *
from get_response_data import _response_data
import math

_inp = GeCCo_Input()

# Get the name of the package GeCCo uses the integrals from 
_env = _inp.env

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

print 'frequency for the second order proeprties is set to:', _freq

#Define the perturbation operator. This is the first operator that has been 
#defined. For second order property one more pertubation is involved that 
#takes part while during the evaluation of the property. This code now 
#works only for diagonal elements of the second order properties. So 
#Defining one pertubation is still okay. 

new_target('V(1)')

DEF_HAMILTONIAN({LABEL:'V(1)',MAX_RANK:1})

new_target('IMPORT_V(1)')

_op_list={'V(1)':'ME_V(1)'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:1}) # For second order properties, the perturbations are mainly singlet. 

#Importing the integrals: It is hard coded to be a dipole operator in z-direction. 
##TODO##  ###NEED TO BE GENERALIZED###

IMPORT({LIST:'ME_V(1)',TYPE:'ZDIPLEN',ENV:_env})

#Loop over n_par to define operators corresponding to both negative and positive frequencies, if necessary.
#Operators with structure of the cluster operators are necessary here.

for i in range(0,n_par):

    i_par = str(i+1)

    new_target('T(1)'+i_par)
    depend('T')

    CLONE_OPERATOR({LABEL:'T(1)'+i_par,
                    TEMPLATE:'T'})

#Setting the order of the T(1) operator to be first

    SET_ORDER({LABEL:'T(1)'+i_par,
              ORDER:1,
              SPECIES:1})

    new_target('RSPNS(1)_OP'+i_par)

    depend('T','C0','OMG','A_C0')

    _op_list={'Ttr(1)'+i_par:'T',
              'ST(1)'+i_par:'OMG',
              'O(1)_T'+i_par:'OMG',
              'O(1)LT'+i_par:'OMG',
              'O(1)RT'+i_par:'OMG',
              'DIAG_T(1)'+i_par:'T'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op]})

    DEF_SCALAR({LABEL:'RED_LAG(1)_T'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RT'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LT'+i_par})
    DEF_SCALAR({LABEL:'den12'+i_par})

#Setting formula corresponding to the left hand side of the equations: 
#\bra{{\Psi}_0^{(0)}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \frac{\partial \bar{H}_0}
#{\partial \boldsymbol{t^{\prime0}}} \ket{\Psi_0^{(0)}} \boldsymbol{t^{\prime}}(\omega)

    new_target('F_O(1)LT'+i_par)

    depend('T(1)'+i_par)
    depend('RSPNS(1)_OP'+i_par,'F_MRCC_LAG')

    DERIVATIVE({LABEL_RES:'F_preRED_LAG(1)LT'+i_par,
                LABEL_IN:'F_MRCC_LAG',
                OP_RES:'RED_LAG(1)LT'+i_par,
                OP_DERIV:'L',
                OP_MULT:'L'})

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LT'+i_par,
                LABEL_IN:'F_preRED_LAG(1)LT'+i_par,
                OP_RES:'RED_LAG(1)LT'+i_par,
                OP_DERIV:'T',
                OP_MULT:'T(1)'+i_par})

    DERIVATIVE({LABEL_RES:'F_O(1)LT'+i_par,
                LABEL_IN:'F_RED_LAG(1)LT'+i_par,
                OP_RES:'O(1)LT'+i_par,
                OP_DERIV:'L'})

#Setting formula corresponding to the right hand side of the equations: 
#\bra{{\Psi}_0^{(0)}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \bar{V}(\omega) \ket{\Psi_0^{(0)}}

    new_target('F_O(1)RT'+i_par)

    depend ('F_O(1)LT'+i_par)
    depend('V(1)')

    REPLACE({LABEL_RES:'F_RED_LAG(1)RT'+i_par,
             LABEL_IN:'F_preRED_LAG(1)LT'+i_par,
             OP_LIST:['H','V(1)']})

    INVARIANT({LABEL_RES:'F_RED_LAG(1)RT'+i_par,
               LABEL_IN:'F_RED_LAG(1)RT'+i_par,
               OP_RES:'RED_LAG(1)RT'+i_par,
               OPERATORS:'H'})

    DERIVATIVE({LABEL_RES:'F_O(1)RT'+i_par,
                LABEL_IN:'F_RED_LAG(1)RT'+i_par,
                OP_RES:'O(1)RT'+i_par,
                OP_DERIV:'L'})

#Setting formula for transforming the frist order cluster amplitude between linearly independent and 
#dependent basis. Here the formula for the zeroeth order T (F_T) has been used to start with. 

    new_target('F_T(1)'+i_par)

    depend('F_T')

    INVARIANT({LABEL_RES:'F_T(1)'+i_par,
               LABEL_IN:'F_T',
               OP_RES:'T(1)'+i_par,
               OPERATORS:'H'})

    REPLACE({LABEL_RES:'F_T(1)'+i_par,LABEL_IN:'F_T(1)'+i_par,
              OP_LIST:['Ttr','Ttr(1)'+i_par]})

#Setting formula for the metric matrix:
#\bra{{\Psi}_0^{(0)}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \boldsymbol{\tilde{\tau}^{\prime}} \ket{\Psi_0^{(0)}} \boldsymbol{t^{\prime}}(\omega)
#The {\tilde{\tau}^{\prime} operator is then expanded in terms of simple cluster operator \tau
#This is the same definition of the metric matrix that has been used for the excitation energy calculation

    new_target('F_ST(1)'+i_par)

    _form = 'F_den12'+i_par
    _den = 'den12'+i_par
    _op_t = 'T(1)'+i_par

    _expand_product_basis={LABEL:_form,
                           NEW:True,
                           OP_RES:_den}
    
    _ops_contract={OPERATORS:[_den,'C0^+','L',_op_t,'C0',_den],
                   IDX_SV:[1,2,3,4,5,1]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))
    
    _expand_product_basis={LABEL:_form,
                           NEW:False,
                           OP_RES:_den}

    _ops_contract={OPERATORS:[_den,'C0^+','L',_op_t,'T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,1],
                   FAC:0.5}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))
    
    _ops_contract={OPERATORS:[_den,'C0^+','L','T',_op_t,'C0',_den],
                   IDX_SV:[1,2,3,4,5,6,1],
                   FAC:-0.5}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L',_op_t,'T','T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:1.0/6}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L','T',_op_t,'T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:-0.33333333333333330}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))
    
    _ops_contract={OPERATORS:[_den,'C0^+','L','T','T',_op_t,'C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:0.16666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L',_op_t,'T','T','T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:0.041666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L','T',_op_t,'T','T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:-0.125}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L','T','T',_op_t,'T','C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:0.125}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:[_den,'C0^+','L','T','T','T',_op_t,'C0',_den],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:-0.041666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))
    
    SELECT_SPECIAL({LABEL_RES:_form,
                    LABEL_IN:_form,
                    TYPE:'nonzero',
                    MODE:'sum'})

    DERIVATIVE({LABEL_RES:'F_ST(1)'+i_par,
                LABEL_IN:_form,
                OP_RES:'ST(1)'+i_par,
                OP_DERIV:'L'})


#setting frequencies for the first order cluster operators. This feature is still not used later, but can be used later to make life lot simpler
    _freq_sign = _freq*math.pow(-1,i)

#Defining the me-list for T(1)
    new_target('DEF_ME_T(1)'+i_par)
    depend('T(1)'+i_par)

    DEF_ME_LIST({LIST:'ME_T(1)'+i_par,
                 OPERATOR:'T(1)'+i_par,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:1})

    SET_FREQ({LIST:'ME_T(1)'+i_par,
              FREQ:_freq_sign})

#defining the me-list for rest of the new operators.
    new_target('LIST_RSPNS(1)_OP'+i_par)

    depend ('F_O(1)RT'+i_par)

    _op_list={'ST(1)'+i_par:'ME_ST(1)'+i_par,
              'O(1)LT'+i_par:'ME_O(1)LT'+i_par,
              'O(1)RT'+i_par:'ME_O(1)RT'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    _op_list={'Ttr(1)'+i_par:'ME_Ttr(1)'+i_par,
              'DIAG_T(1)'+i_par:'ME_DIAG_T(1)'+i_par}
    
    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0})

#optimizing the formula:

    new_target('OPT_RSPNS_T(1)'+i_par)

    depend('H0','DEF_ME_T','DEF_ME_E(MR)')
    depend('DEF_ME_T(1)'+i_par)
    depend('LIST_RSPNS(1)_OP'+i_par)
    depend('F_T(1)'+i_par,'DEF_ME_Dtrdag')
    depend('F_ST(1)'+i_par)
    depend('IMPORT_V(1)')

    OPTIMIZE({LABEL_OPT:'FOPT_RSPNS_T(1)'+i_par,
              LABELS_IN:['F_O(1)LT'+i_par,'F_O(1)RT'+i_par,'F_ST(1)'+i_par]})
  
    OPTIMIZE({LABEL_OPT:'FOPT_T(1)'+i_par,
              LABELS_IN:'F_T(1)'+i_par})

#Getting the preconditioner by copying the one of the zeroeth order

    new_target('DIAG_T(1)'+i_par)

    depend('DIAG1SxxM00_T')

    SCALE_COPY({LIST_RES:'ME_DIAG_T(1)'+i_par,LIST_INP:'DIAG1SxxM00_T',FAC:1.0})

#Solving the linear equation:

    if (_restart<3):
        new_target('SOLVE_T(1)'+i_par,True)
    else:
        new_target('SOLVE_T(1)'+i_par)

    depend('OPT_RSPNS_T(1)'+i_par)
    depend('DIAG_T(1)'+i_par)

    SOLVE_LEQ({LIST_OPT:'ME_T(1)'+i_par,
               LIST_PRC:'ME_DIAG_T(1)'+i_par,
               OP_MVP:'O(1)LT'+i_par,
               OP_SVP:'ST(1)'+i_par,
               OP_RHS:'O(1)RT'+i_par,
               FORM:'FOPT_RSPNS_T(1)'+i_par,
               LIST_SPC:['ME_T(1)'+i_par,'ME_Ttr(1)'+i_par,'ME_Dtr','ME_Dtrdag'],
               FORM_SPC:'FOPT_T(1)'+i_par,
               MODE:'TRF',
               N_ROOTS:1})
