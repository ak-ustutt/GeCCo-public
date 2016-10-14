#This is a python interface for GeCCo to solve the first order T(\omega) to 
#calculate static and dynamic polarizability using the linear response method.
#This particular interface started to be written on September, 2015

from gecco_interface import *
from get_response_data import _response_data, _pop_data, _cmp_data
import math

_inp = GeCCo_Input()

# Get the name of the package GeCCo uses the integrals from 
_env = _inp.env

# _multd2h will be required later to find out the spatial symmetry of 'T(1)' ####

_multd2h = [[1,2,3,4,5,6,7,8],
[2,1,4,3,6,5,8,7],
[3,4,1,2,7,8,5,6],
[4,3,2,1,8,7,6,5],
[5,6,7,8,1,2,3,4],
[6,5,8,7,2,1,4,3],
[7,8,5,6,3,4,1,2],
[8,7,6,5,4,3,2,1]]

_orb = Orb_Info()

_s2 = _orb.get('imult')

_ms = _orb.get('ims')

_isym_0 = _orb.get('lsym')
_isym= _isym_0

if ((_ms == 0) and ((_s2-1 % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and ((_s2+1 % 4) == 0)):
    _msc = -1
else:
    _msc = 0

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

if _option == 1:
    print 'Including the disputed term :D'
elif _option == 2:
    print 'Excluding the disputed term :D'
else:
    quit_error('Input error: unknown option for ic-MRCC properties') 
    
print 'frequency for the second order proeprties is set to:', _freq

#Define the perturbation operator. This is the first operator that has been 
#defined. For second order property one more pertubation is involved that 
#takes part while during the evaluation of the property. This code now 
#works only for diagonal elements of the second order properties. So 
#Defining one pertubation is still okay. 

new_target('V(1)')

DEF_HAMILTONIAN({LABEL:'V(1)',MAX_RANK:1})

new_target('IMPORT_V(1)')

depend('V(1)')

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


# We define all the parameters from earlier calculation that will be needed here

new_target('DEF_RSPNS(1)')
depend('DEF_ME_L')
depend('DEF_ME_C0_bar')
depend('IMPORT_V(1)')

DEF_SCALAR({LABEL:'RSPNS(1)'})

_op_list={'RSPNS(1)':'ME_RSPNS(1)'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

new_target('GET_RSPNS(1)_FORM')

depend('F_MRCC_LAG')
depend('DEF_RSPNS(1)')

DEF_SCALAR({LABEL:'preRSPNS(1)_1'})
DEF_SCALAR({LABEL:'preRSPNS(1)_2'})
DEF_SCALAR({LABEL:'preRSPNS(1)'})

INVARIANT({LABEL_RES:'F_preRSPNS(1)_1',
           LABEL_IN:'F_MRCC_LAG',
           OP_RES:'preRSPNS(1)_1',
           OPERATORS:'L'})

REPLACE({LABEL_RES:'F_preRSPNS(1)_1',
         LABEL_IN:'F_preRSPNS(1)_1',
         OP_LIST:['C0^+','C0_bar']})

DERIVATIVE({LABEL_RES:'F_preRSPNS(1)_2',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'preRSPNS(1)_2',
            OP_DERIV:'L',
            OP_MULT:'L'})

DEF_FORMULA({LABEL:'F_preRSPNS(1)',
             FORMULA:'preRSPNS(1)=preRSPNS(1)_1+preRSPNS(1)_2'})

EXPAND({LABEL_RES:'F_preRSPNS(1)',
        LABEL_IN:'F_preRSPNS(1)',
        INTERM:['F_preRSPNS(1)_1','F_preRSPNS(1)_2']})

REPLACE({LABEL_RES:'F_RSPNS(1)',
         LABEL_IN:'F_preRSPNS(1)',
         OP_LIST:['H','V(1)']})

INVARIANT({LABEL_RES:'F_RSPNS(1)',
           LABEL_IN:'F_RSPNS(1)',
           OP_RES:'RSPNS(1)',
           OPERATORS:'H'})

new_target('EVAL_RSPNS(1)')

depend('H0','DEF_ME_T','DEF_ME_C0')
depend('GET_RSPNS(1)_FORM')

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(1)',
          LABELS_IN:'F_RSPNS(1)'})

EVALUATE({FORM:'FOPT_RSPNS(1)'})

PRINT_MEL({LIST:'ME_RSPNS(1)'})

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

    new_target('C0(1)'+i_par)
    depend('C0')

    CLONE_OPERATOR({LABEL:'C0(1)'+i_par,
                    TEMPLATE:'C0'})

#Setting the order of the C0(1) operator to be first

    SET_ORDER({LABEL:'C0(1)'+i_par,
              ORDER:1,
              SPECIES:1})

#Defining all the other operators needed for this calculation
    new_target('RSPNS(1)_OP'+i_par)

    depend('T','C0','OMG','A_C0')

    _op_list={'Ttr(1)'+i_par:'T',
              'ST(1)'+i_par:'OMG',
              'O(1)_T'+i_par:'OMG',
              'O(1)LT'+i_par:'OMG',
              'O(1)RT'+i_par:'OMG',
              'O(1)LC'+i_par:'A_C0',
              'O(1)LCM'+i_par:'A_C0',
              'O(1)RC'+i_par:'A_C0',
              'O(1)RCM'+i_par:'A_C0',
              'DIAG_T(1)'+i_par:'T'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op]})

    _op_list={'C0(1)_x'+i_par:'C0'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op],
                        ADJOINT:True})

    DEF_SCALAR({LABEL:'RED_LAG(1)_T'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)_C'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RT'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LT'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RC'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LC'+i_par})
    DEF_SCALAR({LABEL:'den12'+i_par})
    DEF_SCALAR({LABEL:'FREQ'+i_par})

# setting the formula for the left hand side of the linear equation to solve T(1):
# \boldsymbol{c}^*(-\omega) \bra{\boldsymbol{\Phi}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \bar{H}_0 \ket{\Psi_0^{(0)}} 
# + \bra{{\Psi}_0^{(0)}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \bar{H}_0 \ket{\boldsymbol{\Phi}} \boldsymbol{c}(\omega)
# + \bra{{\Psi}_0^{(0)}} \boldsymbol{\hat{\tau}^{\prime\dagger}} \frac{\partial \bar{H}_0} {\partial \boldsymbol{t^{\prime0}}} \ket{\Psi_0^{(0)}} \boldsymbol{t^{\prime}}(\omega)

    new_target('F_O(1)LT'+i_par)

    depend('T(1)'+i_par,'C0(1)'+i_par)
    depend('RSPNS(1)_OP'+i_par,'F_MRCC_LAG')

    DERIVATIVE({LABEL_RES:'F_preRED_LAG(1)LT'+i_par,
                LABEL_IN:'F_MRCC_LAG',
                OP_RES:'RED_LAG(1)LT'+i_par,
                OP_DERIV:'L',
                OP_MULT:'L'})

    if _option == 1:

# option=1 has all the three terms from last comment

        DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LT'+i_par,
                    LABEL_IN:'F_preRED_LAG(1)LT'+i_par,
                    OP_RES:'RED_LAG(1)LT'+i_par,
                    OP_DERIV:['T','C0','C0^+'],
                    OP_MULT:['T(1)'+i_par,'C0(1)'+i_par,'C0(1)_x'+i_par]})

        REPLACE({LABEL_RES:'F_RED_LAG(1)LT'+i_par,
                 LABEL_IN:'F_RED_LAG(1)LT'+i_par,
                 OP_LIST:['C0(1)_x'+i_par,'C0(1)'+i_par+'^+']})

    elif _option == 2:

# option=2 excludes the first term

        DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LT'+i_par,
                    LABEL_IN:'F_preRED_LAG(1)LT'+i_par,
                    OP_RES:'RED_LAG(1)LT'+i_par,
                    OP_DERIV:['T','C0'],
                    OP_MULT:['T(1)'+i_par,'C0(1)'+i_par]})

    DERIVATIVE({LABEL_RES:'F_O(1)LT'+i_par,
                LABEL_IN:'F_RED_LAG(1)LT'+i_par,
                OP_RES:'O(1)LT'+i_par,
                OP_DERIV:'L'})

# setting the formula for the left hand side of the linear equation to solve C0(1):
# \bra{\boldsymbol{\Phi}} (\bar{H}_0 - E_0) \ket{\boldsymbol{\Phi}} \boldsymbol{c}(\omega)
# + \bra{\boldsymbol{\Phi}} \frac{\partial \bar{H}_0}{\partial \boldsymbol{t^{\prime0}}} \ket{\Psi_0^{(0)}} \boldsymbol{t^{\prime }}(\omega)

    new_target('F_O(1)LC'+i_par)

    depend('RSPNS(1)_OP'+i_par,'F_MRCC_LAG','E(MR)')

    INVARIANT({LABEL_RES:'F_preRED_LAG(1)LC'+i_par,
               LABEL_IN:'F_MRCC_LAG',
               OP_RES:'RED_LAG(1)_C'+i_par,
               OPERATORS:'L'})

    EXPAND_OP_PRODUCT({LABEL:'F_preRED_LAG(1)LC'+i_par,NEW:False,
                       OP_RES:'RED_LAG(1)LC'+i_par,
                       OPERATORS:['E(MR)','C0^+','C0'],
                       IDX_SV:[1,2,3],
                       FAC:-1.0})

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LC'+i_par,
                LABEL_IN:'F_preRED_LAG(1)LC'+i_par,
                OP_RES:'RED_LAG(1)LC'+i_par,
                OP_DERIV:['T','C0'],
                OP_MULT:['T(1)'+i_par,'C0(1)'+i_par]})

    DERIVATIVE({LABEL_RES:'F_O(1)LC'+i_par,
                LABEL_IN:'F_RED_LAG(1)LC'+i_par,
                OP_RES:'O(1)LC'+i_par,
                OP_DERIV:['C0^+']})

    #PRINT_FORMULA({LABEL:'F_RED_LAG(1)LC'+i_par})
    #PRINT_FORMULA({LABEL:'F_O(1)LC'+i_par})

# The component of C0 are projected out from the left hand side of C0(1)

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_LC(1)'+i_par,
                       OP_RES:'O(1)LCM'+i_par,
                       NEW:True,
                       OPERATORS:['O(1)LCM'+i_par,'O(1)LC'+i_par,'O(1)LCM'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_LC(1)'+i_par,
                       OP_RES:'O(1)LCM'+i_par,
                       NEW:False,
                       OPERATORS:['O(1)LCM'+i_par,'C0','C0^+','O(1)LC'+i_par,'O(1)LCM'+i_par],
                       IDX_SV:[1,2,3,4,1],
                       AVOID:[1,4,3,5],
                       FAC:-1.0})

    EXPAND({LABEL_RES:'F_PRJ_LC(1)'+i_par,
            LABEL_IN:'F_prePRJ_LC(1)'+i_par,
            INTERM:'F_O(1)LC'+i_par})

    #PRINT_FORMULA({LABEL:'F_PRJ_LC(1)'+i_par})

# setting the formula for the right hand side of the linear equation to solve T(1):
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

# setting the formula for the right hand side of the linear equation to solve C0(1):
#\bra{\boldsymbol{\Phi}} \bar{V}(\omega) \ket{\Psi_0^{(0)}}

    new_target('F_O(1)RC'+i_par)

    depend('EVAL_RSPNS(1)')

    depend ('F_O(1)LC'+i_par)
    depend('V(1)')

    REPLACE({LABEL_RES:'F_RED_LAG(1)RC'+i_par,
             LABEL_IN:'F_preRED_LAG(1)LC'+i_par,
             OP_LIST:['H','V(1)']})

    INVARIANT({LABEL_RES:'F_RED_LAG(1)RC'+i_par,
               LABEL_IN:'F_RED_LAG(1)RC'+i_par,
               OP_RES:'RED_LAG(1)RC'+i_par,
               OPERATORS:['H','E(MR)']})

    EXPAND_OP_PRODUCT({LABEL:'F_RED_LAG(1)RC'+i_par,NEW:False,
                       OP_RES:'RED_LAG(1)RC'+i_par,
                       OPERATORS:['RSPNS(1)','C0^+','C0'],
                       IDX_SV:[1,2,3],
                       FAC:-1.0})

    DERIVATIVE({LABEL_RES:'F_O(1)RC'+i_par,
                LABEL_IN:'F_RED_LAG(1)RC'+i_par,
                OP_RES:'O(1)RC'+i_par,
                OP_DERIV:['C0^+']})

    #PRINT_FORMULA({LABEL:'F_O(1)RC'+i_par})

# The component of C0 are projected out from the right hand side of C0(1)

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_RC(1)'+i_par,
                       OP_RES:'O(1)RCM'+i_par,
                       NEW:True,
                       OPERATORS:['O(1)RCM'+i_par,'O(1)RC'+i_par,'O(1)RCM'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_RC(1)'+i_par,
                       OP_RES:'O(1)RCM'+i_par,
                       NEW:False,
                       OPERATORS:['O(1)RCM'+i_par,'C0','C0^+','O(1)RC'+i_par,'O(1)RCM'+i_par],
                       IDX_SV:[1,2,3,4,1],
                       AVOID:[1,4,3,5],
                       FAC:-1.0})

    EXPAND({LABEL_RES:'F_PRJ_RC(1)'+i_par,
            LABEL_IN:'F_prePRJ_RC(1)'+i_par,
            INTERM:'F_O(1)RC'+i_par})

    #PRINT_FORMULA({LABEL:'F_PRJ_RC(1)'+i_par})

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

#This is another way to project out the C0 components, but this time from 
#C0(1) itself while solving it though linear solver

    new_target('F_PRJ_C(1)'+i_par)

    EXPAND_OP_PRODUCT({LABEL:'F_PRJ_C(1)'+i_par,
                       OP_RES:'C0(1)'+i_par,
                       OPERATORS:['C0(1)'+i_par,'C0','C0^+','C0(1)'+i_par,'C0(1)'+i_par],
                       IDX_SV:[1,2,3,4,1],
                       AVOID:[1,4,3,5],
                       FAC:-1.0})

    #PRINT_FORMULA({LABEL:'F_PRJ_C(1)'+i_par})

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

#Getting the signed frequency that would be needed to initialize frequency for the 
#parameters we are going to solve

    _freq_sign = _freq*math.pow(-1,i)

    new_target('DEF_ME_T(1)'+i_par)
    depend('T(1)'+i_par)

    DEF_ME_LIST({LIST:'ME_T(1)'+i_par,
                 OPERATOR:'T(1)'+i_par,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:1})

    SET_FREQ({LIST:'ME_T(1)'+i_par,
              FREQ:_freq_sign})

    new_target('DEF_ME_C0(1)'+i_par)
    depend('C0(1)'+i_par)

    DEF_ME_LIST({LIST:'ME_C0(1)'+i_par,
                 OPERATOR:'C0(1)'+i_par,
                 IRREP:_isym,
                 '2MS':_ms,
                 AB_SYM:_msc})

    SET_FREQ({LIST:'ME_C0(1)'+i_par,
              FREQ:_freq_sign})

    new_target('LIST_RSPNS(1)_OP'+i_par)

    depend ('F_O(1)RT'+i_par)
    depend ('F_O(1)RC'+i_par)

    _op_list={'ST(1)'+i_par:'ME_ST(1)'+i_par,
              'O(1)LT'+i_par:'ME_O(1)LT'+i_par,
              'O(1)RT'+i_par:'ME_O(1)RT'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    _op_list={'O(1)LC'+i_par:'ME_O(1)LC'+i_par,
              'O(1)LCM'+i_par:'ME_O(1)LCM'+i_par,
              'O(1)RC'+i_par:'ME_O(1)RC'+i_par,
              'O(1)RCM'+i_par:'ME_O(1)RCM'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_isym,
                     '2MS':_ms,
                     AB_SYM:_msc})

    _op_list={'Ttr(1)'+i_par:'ME_Ttr(1)'+i_par,
              'DIAG_T(1)'+i_par:'ME_DIAG_T(1)'+i_par,
              'FREQ'+i_par:'ME_FREQ'+i_par}
    
    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0})

# Initializing the ME_FREQ for this i_par
    SET_MEL({LIST:'ME_FREQ'+i_par,
             IDX_LIST:1,
             VAL_LIST:-1.0*_freq_sign})

    new_target('OPT_RSPNS_T(1)'+i_par)

    depend('H0','DEF_ME_T','DEF_ME_E(MR)')
    depend('DEF_ME_T(1)'+i_par)
    depend('DEF_ME_C0(1)'+i_par)
    depend('LIST_RSPNS(1)_OP'+i_par)
    depend('F_T(1)'+i_par,'DEF_ME_Dtrdag')
    depend('F_ST(1)'+i_par)
    depend('IMPORT_V(1)')
    depend('F_PRJ_C(1)'+i_par)

    OPTIMIZE({LABEL_OPT:'FOPT_RSPNS_T(1)'+i_par,
 #            LABELS_IN:['F_O(1)LT'+i_par,'F_O(1)RT'+i_par,'F_ST(1)'+i_par]})
 #            LABELS_IN:['F_O(1)LT'+i_par,'F_O(1)RT'+i_par,'F_O(1)LC'+i_par,'F_O(1)RC'+i_par,'F_ST(1)'+i_par]})
              LABELS_IN:['F_O(1)LT'+i_par,'F_O(1)RT'+i_par,'F_PRJ_LC(1)'+i_par,'F_PRJ_RC(1)'+i_par,'F_ST(1)'+i_par]})
  
    OPTIMIZE({LABEL_OPT:'FOPT_T(1)'+i_par,
              LABELS_IN:'F_T(1)'+i_par})

    OPTIMIZE({LABEL_OPT:'FOPT_PRJ_C(1)'+i_par,
              LABELS_IN:'F_PRJ_C(1)'+i_par})

    new_target('DIAG_T(1)'+i_par)

    depend('DIAG1SxxM00_T')

    SCALE_COPY({LIST_RES:'ME_DIAG_T(1)'+i_par,LIST_INP:'DIAG1SxxM00_T',FAC:1.0})


    new_target('DIAG_C0(1)'+i_par)

    depend('H0', 'DEF_ME_E(MR)')

    _op_list={'DIAG_C0(1)'+i_par:'C0'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op]})

    _op_list={'DIAG_C0(1)'+i_par:'ME_DIAG_C0(1)'+i_par}
    
    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_isym,
                     '2MS':_ms})

    DEF_ME_LIST({LIST:'ME_MINEN'+i_par,
                 OPERATOR:'E(MR)',
                 IRREP:1,
                 '2MS':0})

    ASSIGN_ME2OP({LIST:'ME_E(MR)',OPERATOR:'E(MR)'})

    PRECONDITIONER({LIST_PRC:'ME_DIAG_C0(1)'+i_par,
                    LIST_INP:'H0',
                    MODE:'dia-H'})

    SCALE_COPY({LIST_RES:'ME_MINEN'+i_par,
                LIST_INP:'ME_E(MR)',
                FAC:-1.0})
 
    EXTRACT_DIAG({LIST_RES:'ME_DIAG_C0(1)'+i_par,
                  LIST_IN:'ME_MINEN'+i_par,
                  MODE:'ext_act'})
    #PRINT_MEL({LIST:'ME_DIAG_C0(1)'+i_par})

    #EXTRACT_DIAG({LIST_RES:'ME_DIAG_C0(1)'+i_par,
    #              LIST_IN:'ME_FREQ'+i_par,
    #              MODE:'ext_act'})
    #PRINT_MEL({LIST:'ME_DIAG_C0(1)'+i_par})

#Solving the linear equation
    if (_restart<3):
        new_target('SOLVE_T(1)'+i_par,True)
    else:
        new_target('SOLVE_T(1)'+i_par)

    depend('OPT_RSPNS_T(1)'+i_par)
    depend('DIAG_T(1)'+i_par)
    depend('DIAG_C0(1)'+i_par)

    _solve_leq_arg = {}
    _solve_leq_arg[LIST_OPT] = ['ME_T(1)'+i_par,'ME_C0(1)'+i_par]
    _solve_leq_arg[LIST_PRC] = ['ME_DIAG_T(1)'+i_par,'ME_DIAG_C0(1)'+i_par]
    #_solve_leq_arg[OP_MVP]   = ['O(1)LT'+i_par,'O(1)LC'+i_par]
    _solve_leq_arg[OP_MVP]   = ['O(1)LT'+i_par,'O(1)LCM'+i_par]
    _solve_leq_arg[OP_SVP]   = ['ST(1)'+i_par,'C0(1)'+i_par]
    #_solve_leq_arg[OP_RHS]   = ['O(1)RT'+i_par,'O(1)RC'+i_par]
    _solve_leq_arg[OP_RHS]   = ['O(1)RT'+i_par,'O(1)RCM'+i_par]
    _solve_leq_arg[FORM]     = 'FOPT_RSPNS_T(1)'+i_par
    _solve_leq_arg[LIST_SPC] = ['ME_T(1)'+i_par,'ME_Ttr(1)'+i_par,'ME_Dtr','ME_Dtrdag']
    _solve_leq_arg[FORM_SPC] = ['FOPT_T(1)'+i_par,'FOPT_PRJ_C(1)'+i_par]
    #_solve_leq_arg[MODE]     = 'TRF PRJ'
    _solve_leq_arg[MODE]     = 'TRF DIA'
    _solve_leq_arg[N_ROOTS]  = 1

    SOLVE_LEQ(_solve_leq_arg)

    PRINT_MEL({LIST:'ME_T(1)'+i_par,
               COMMENT:'First order T amplitude'})

    PRINT_MEL({LIST:'ME_C0(1)'+i_par,
               COMMENT:'First order C0 coefficient'})

#orthogonalisation of the CO(1) respect to C(0), need to check if this is redundant!
    if (_restart<3):
        new_target('ORTH_C0(1)'+i_par,True)
    else:
        new_target('ORTH_C0(1)'+i_par)


    _op_list={'C0(1)_orth_'+i_par:'C0'}

    for _op in _op_list:
        CLONE_OPERATOR({LABEL:_op,
                        TEMPLATE:_op_list[_op]})

    EXPAND_OP_PRODUCT({LABEL:'F_ORTH_C(1)'+i_par,
                       OP_RES:'C0(1)_orth_'+i_par,
                       NEW:True,
                       OPERATORS:['C0(1)_orth_'+i_par,'C0(1)'+i_par,'C0(1)_orth_'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_ORTH_C(1)'+i_par,
                       OP_RES:'C0(1)_orth_'+i_par,
                       NEW:False,
                       OPERATORS:['C0(1)_orth_'+i_par,'C0','C0^+','C0(1)'+i_par,'C0(1)_orth_'+i_par],
                       IDX_SV:[1,2,3,4,1],
                       AVOID:[1,4,3,5],
                       FAC:-1.0})

    _op_list={'C0(1)_orth_'+i_par:'ME_C0(1)_orth_'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    OPTIMIZE({LABEL_OPT:'FOPT_ORTH_C(1)'+i_par,
              LABELS_IN:['F_ORTH_C(1)'+i_par]})

    EVALUATE({FORM:'FOPT_ORTH_C(1)'+i_par})

    PRINT_MEL({LIST:'ME_C0(1)_orth_'+i_par})
    SCALE_COPY({LIST_RES:'ME_C0(1)'+i_par,LIST_INP:'ME_C0(1)_orth_'+i_par,FAC:1.0})
    #COPY_LIST({LIST_RES:'ME_C0(1)'+i_par,LIST_INP:'ME_C0(1)_orth_'+i_par})

