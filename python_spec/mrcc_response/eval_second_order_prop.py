#This is python interface to evaluate the second order properties using the all the zeroeth 
#order paramaters and the first order cluster amplitudes and ci coefficients.

from python_interface.gecco_interface import *
from python_spec.mrcc_response.get_response_data import _response_data, _pop_data, _cmp_data, _calc_data

#Getting the frequency
_freq=_response_data['freq']
#Getting the value of the restart option
_restart=_response_data['restart']

#Getting the total number of perturbation operator need to be defined 
_npop=_response_data['nPop']
#Getting total number of response calculation
_ncnt=_response_data['nCnt']
#Getting the maximum order of the properties that will be calculated
_maxord=_response_data['maxorder']

#Getting the option from response data
_option=_response_data['option']

#preprocessing the second order response function. First changing the GeCCo specific Lagrangian 
#to the actual Lagrangian of ic-MRCC and then replacing the electronic Hamiltonian with the Perturbation
#NOTE: This perturbation is still the previous one. It should be a new one for a non-diagonal component of the properties

new_target('F_MRCC_LAG_PROP')

depend('F_MRCC_LAG')
depend('L')
depend('C0_bar')

DEF_SCALAR({LABEL:'LAG_PROP'})
DEF_SCALAR({LABEL:'preRSPNS(2)_1'})
DEF_SCALAR({LABEL:'preRSPNS(2)_2'})

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

DEF_FORMULA({LABEL:'F_MRCC_LAG_PROP',
             FORMULA:'LAG_PROP=preRSPNS(2)_1+preRSPNS(2)_2'})

EXPAND({LABEL_RES:'F_MRCC_LAG_PROP',
        LABEL_IN:'F_MRCC_LAG_PROP',
        INTERM:['F_preRSPNS(2)_1','F_preRSPNS(2)_2']})

new_target('F_preRSPNS(2)')

depend('F_MRCC_LAG_PROP')
depend('EVAL_RSPNS(1)')

DEF_SCALAR({LABEL:'preRSPNS(2)'})
REPLACE({LABEL_RES:'F_preRSPNS(2)',
         LABEL_IN:'F_MRCC_LAG_PROP',
         OP_LIST:['H','V(1)']})

INVARIANT({LABEL_RES:'F_preRSPNS(2)',
           LABEL_IN:'F_preRSPNS(2)',
           OP_RES:'preRSPNS(2)',
           OPERATORS:'H'})



_formula_to_add_1=[]
_formula_to_add_2=[]
_formula_to_add_3=[]
_formula_to_add_4=[]

_interm_to_add_1=[]
_interm_to_add_2=[]
_interm_to_add_3=[]
_interm_to_add_4=[]

n_par=len(_calc_data[0]['conj_comp'])

for i in range(0,n_par):

    i_par = str(i+1)

    _pop_idx = _calc_data[0]['conj_prop'][i]-1

    _cur_ext=_pop_data['name'][_pop_idx]+_pop_data['comp'][_pop_idx]
    _pop_name='V'+_cur_ext

    ext_1=str(_calc_data[0]['prop_comp'][i])
    t_op_1='T(1)'+ext_1
    c_op_1='C0(1)'+ext_1
    ext_2=str(_calc_data[0]['conj_comp'][i])
    t_op_2='T(1)'+ext_2
    c_op_2='C0(1)'+ext_2

#Getting the part of the second order response function which is a product of perturbation and 
#first order wave function paramaeters, T(1) and C(1). It Corresponds the scalar RSPNS(2)_1

    new_target('F_RSPNS(2)_1'+i_par)

    depend('F_preRSPNS(2)')
    depend('EVAL_RSPNS(1)')
    depend('DEF_ME_C0_bar')
    depend('DEF_ME_C0')
    depend('DEF_ME_T')
    depend('DEF_ME_'+t_op_1)
    depend('DEF_ME_'+t_op_2)
    depend('DEF_ME_'+c_op_1)
    depend('DEF_ME_'+c_op_2)

    DEF_SCALAR({LABEL:'RSPNS(2)_1'+i_par})

    _deriv_arg = {}
    _deriv_arg[LABEL_RES] = 'F_RSPNS(2)_1'+i_par
    _deriv_arg[LABEL_IN]  = 'F_preRSPNS(2)'
    _deriv_arg[OP_RES]    = 'RSPNS(2)_1'+i_par

    _replace_arg = {}
    _replace_arg[LABEL_RES] = 'F_RSPNS(2)_1'+i_par
    _replace_arg[LABEL_IN]  = 'F_RSPNS(2)_1'+i_par
    _replace_arg[OP_LIST]   = ['C0(1)_x1',c_op_1+'^+']

    if _option == 1: 
# option=1 has extra contribution from the 'C0(1)^+'
        _deriv_arg[OP_DERIV] = ['T','C0','C0^+']
        _deriv_arg[OP_MULT]  = [t_op_1,c_op_1,'C0(1)_x1']

        DERIVATIVE(_deriv_arg)

        REPLACE(_replace_arg)

    elif _option == 2: 
        _deriv_arg[OP_DERIV] = ['T','C0']
        _deriv_arg[OP_MULT]  = [t_op_1,c_op_1]

        DERIVATIVE(_deriv_arg)

    else:
        quit_error('Input error: unknown option for ic-MRCC properties') 

    EXPAND_OP_PRODUCT({LABEL:'F_RSPNS(2)_1'+i_par,NEW:False,
                       OP_RES:'RSPNS(2)_1'+i_par,
                       OPERATORS:['RSPNS(1)'+_cur_ext,'C0_bar',c_op_1],
                       IDX_SV:[1,2,3],
                       FAC:-1.0})

    REPLACE({LABEL_RES:'F_RSPNS(2)_1'+i_par,
             LABEL_IN:'F_RSPNS(2)_1'+i_par,
             OP_LIST:['V(1)',_pop_name]})

    PRINT_MEL({LIST:'ME_'+_pop_name})

    _formula_to_add_1.append('F_RSPNS(2)_1'+i_par)
    _interm_to_add_1.append('RSPNS(2)_1'+i_par)
    
#Getting the part of the second order response function which is the product of two first order wave function 
#parameters, T(1) and C0(1). This also includes the cross terms between T(1) and C0(1). It corresponds the scalar RSPNS(2)_2.
#npar=2 would involve the product of T(1) and C0(1) corresponding to different frequencies. 

    new_target('F_RSPNS(2)_2'+i_par)

    DEF_SCALAR({LABEL:'RSPNS(2)_2'+i_par})

    _deriv_arg = {}
    _deriv_arg[LABEL_RES] = 'F_intRSPNS(2)_1'+i_par 
    _deriv_arg[LABEL_IN]  = 'F_MRCC_LAG_PROP' 
    _deriv_arg[OP_RES]    = 'RSPNS(2)_2'+i_par
    _deriv_arg[OP_DERIV]  =  'C0'
    _deriv_arg[OP_MULT] = c_op_1


    DERIVATIVE(_deriv_arg)

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_2'+i_par,
                LABEL_IN:'F_intRSPNS(2)_1'+i_par,
                OP_RES:'RSPNS(2)_2'+i_par,
                OP_DERIV:'T',
                OP_MULT:t_op_2})
    
    _formula_to_add_2.append('F_RSPNS(2)_2'+i_par)
    _interm_to_add_2.append('RSPNS(2)_2'+i_par)

    new_target('F_RSPNS(2)_3'+i_par)

    DEF_SCALAR({LABEL:'RSPNS(2)_3'+i_par})

    _deriv_arg = {}
    _deriv_arg[LABEL_RES] = 'F_intRSPNS(2)_2'+i_par 
    _deriv_arg[LABEL_IN]  = 'F_MRCC_LAG_PROP' 
    _deriv_arg[OP_RES]    = 'RSPNS(2)_3'+i_par
    _deriv_arg[OP_DERIV]  =  'T'
    _deriv_arg[OP_MULT] = t_op_1


    DERIVATIVE(_deriv_arg)

    DERIVATIVE({LABEL_RES:'F_RSPNS(2)_3'+i_par,
                LABEL_IN:'F_intRSPNS(2)_2'+i_par,
                OP_RES:'RSPNS(2)_3'+i_par,
                OP_DERIV:'T',
                OP_MULT:t_op_2})
    
    _formula_to_add_3.append('F_RSPNS(2)_3'+i_par)
    _interm_to_add_3.append('RSPNS(2)_3'+i_par)

    if _option == 1:

        new_target('F_RSPNS(2)_4'+i_par)

        DEF_SCALAR({LABEL:'RSPNS(2)_4'+i_par})

        DERIVATIVE({LABEL_RES:'F_RSPNS(2)_4'+i_par,
                    LABEL_IN:'F_preRSPNS(2)_1',
                    OP_RES:'RSPNS(2)_4'+i_par,
                    OP_DERIV:['T','C0'],
                    OP_MULT:[t_op_1,c_op_1]})

        REPLACE({LABEL_RES:'F_RSPNS(2)_4'+i_par,
                 LABEL_IN:'F_RSPNS(2)_4'+i_par,
                 OP_LIST:['C0^+',c_op_2+'^+']})

        _formula_to_add_4.append('F_RSPNS(2)_4'+i_par)
        _interm_to_add_4.append('RSPNS(2)_4'+i_par)


if _option==1:
    _fac=[1,1,0.5,1]
    _formula_to_add=[_formula_to_add_1, _formula_to_add_2, _formula_to_add_3, _formula_to_add_4]
    _interm_to_add=[_interm_to_add_1, _interm_to_add_2, _interm_to_add_3, _interm_to_add_4]
    _n_formula=4
else:
    _fac=[1,1,0.5]
    _formula_to_add=[_formula_to_add_1, _formula_to_add_2, _formula_to_add_3]
    _interm_to_add=[_interm_to_add_1, _interm_to_add_2, _interm_to_add_3]
    _n_formula=3
    

_target_to_depend=[]

new_target('RSPNS(2)')

DEF_SCALAR({LABEL:'RSPNS(2)'})

DEF_ME_LIST({LIST:'ME_RSPNS(2)',
             OPERATOR:'RSPNS(2)',
             IRREP:1,
             '2MS':0})
    
for i in range (0,_n_formula):

    j = str(i+1)
    new_target('EVAL_RSPNS(2)_'+j)

    depend('RSPNS(2)')
    DEF_SCALAR({LABEL:'RSPNS(2)_'+j})

    _def_formula='RSPNS(2)_'+j+'='

    for k in range (0,n_par):
        i_par = str(k+1)

        print(_formula_to_add[i][k])
        depend(_formula_to_add[i][k])

        _def_formula= _def_formula+'RSPNS(2)_'+j+i_par
        if (k!=(n_par-1)):
            _def_formula=_def_formula+'+'

    DEF_FORMULA({LABEL:'F_RSPNS(2)_'+j,
                 FORMULA:_def_formula})

    EXPAND({LABEL_RES:'F_RSPNS(2)_'+j,
            LABEL_IN:'F_RSPNS(2)_'+j,
            INTERM:_formula_to_add[i]})
    

    DEF_ME_LIST({LIST:'ME_RSPNS(2)_'+j,
                 OPERATOR:'RSPNS(2)_'+j,
                 IRREP:1,
                 '2MS':0})
    
    OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)_'+j,
              LABELS_IN:'F_RSPNS(2)_'+j})

    EVALUATE({FORM:'FOPT_RSPNS(2)_'+j})


    PRINT_MEL({LIST:'ME_RSPNS(2)_'+j,
               FORMAT:'SCAL',
               COMMENT:'Part of the second order property: '})

    EXPAND_OP_PRODUCT({LABEL:'F_RSPNS(2)',
                       NEW:bool(i==0),
                       OP_RES:'RSPNS(2)',
                       OPERATORS:'RSPNS(2)_'+j,
                       IDX_SV:1,
                       FAC:_fac[i]})

    _target_to_depend.append('EVAL_RSPNS(2)_'+j)

new_target('EVAL_RSPNS(2)', True)

for ele in _target_to_depend: 
    depend(ele)

OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(2)',
          LABELS_IN:'F_RSPNS(2)'})

EVALUATE({FORM:'FOPT_RSPNS(2)'})


PRINT_MEL({LIST:'ME_RSPNS(2)',
           FORMAT:'SCAL',
           COMMENT:'Total second order property: '})
