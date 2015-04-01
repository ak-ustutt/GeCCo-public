#
#This is a python interface for GeCCo to calculate excitation energies using ic-MRCEOM-LRT.
#This code is originally written as an experimental file by Pradipta Samanta on October, 2012.
#This particular interface is written on November, 2014 by Pradipta Samanta.

from gecco_interface import *
import numpy as np

_inp = GeCCo_Input()
_orb = Orb_Info()

# _multd2h will be required later to find out the spatial symmetry of what will be called 'R_q' ####

_multd2h = np.matrix('1 2 3 4 5 6 7 8;2 1 4 3 6 5 8 7;3 4 1 2 7 8 5 6;4 3 2 1 8 7 6 5;5 6 7 8 1 2 3 4;6 5 8 7 2 1 4 3;7 8 5 6 3 4 1 2;8 7 6 5 4 3 2 1')

########### constructed _multd2h ############

### Now the necessary information will be obtained from the input and orb_info ###

_method= _inp.get('method.MRCC.excite.method')

if (_method == 'LR' or _method == None):
    _lr_opt = 1
elif (_method == 'EOM'):
    _lr_opt = 2
else:
    quit_error('Input Error: specified wrong method to calculate excitation energy')

_choice= _inp.get('calculate.solve.eigen.guess')
if (_choice == None):
    _choice = 0


_s2_0 = _orb.get('imult')
_isym_0 = _orb.get('lsym')
_ms_0 = _orb.get('ims')
_nsym = _orb.get('nsym')


if ((_ms_0 == 0) and ((_s2_0-1 % 4) == 0)):
    _msc_0 = 1
elif ((_ms_0 == 0) and ((_s2_0+1 % 4) == 0)):
    _msc_0 = -1
else:
    _msc_0 = 0

_sym = _inp.get('calculate.excitation.sym')
_mult = _inp.get('calculate.excitation.mult')

if (_sym == None):
    quit_error('Input missing: spatial symmetry of the excited state has not been set in the input')
else:
    _ncnt = len(_sym)

if (_mult == None):
    quit_error('Input missing: spin multiplicity of the excited state has not been set in the input')

### Some informations to be printed in the output ####

new_target('INPUT_INFO')

PRINT({STRING: 'Doing icMRCC response in' + '    ' +  str(_method) + '      ' + 'framework' })

PRINT({STRING: 'IRREP, S2, Ms of the reference state:' + str(_isym_0) + 
       '       ' + str(_s2_0) + '          ' + str(_msc_0)  })

PRINT({STRING: 'factor for spin-combination' + '    ' + str(_msc_0)})

################################

#################################
# Defining the operators required 
#################################
# 'R_q' is the response vector corresponding to the cluster amplitude T.
# 'R_mu' is the response vector corresponding to the reference coefficient C_mu.
# 'R_prime_q' is the response of T in the independent excitation basis.
# 'AR_rspns_q' and 'AR_rspns_mu' are the matrix vector product contributing to the equations to solve
# 'R_q' and 'R_mu' respectively. 
# 'SR_rspns_q' and 'SR_rspns_mu' are the metric vector product contributing to the equations to solve
# 'R_q' and 'R_mu' respectively. 

new_target('RSPNS_OP')

depend('OMG','D','L','DENS','C0','A_C0','T')

_op_list={'R_prime_q':'T',
          'R_q':'T',
          'R_mu':'C0',
          'AR_rspns_q':'OMG',
          'AR_rspns_mu':'A_C0',
          'SR_rspns_q':'OMG',
          'SR_rspns_mu':'C0'}

for _op in _op_list:
    CLONE_OPERATOR({LABEL:_op,
                    TEMPLATE:_op_list[_op]})

DEF_SCALAR({LABEL:'den12'})


#################################
# Setting up the formula for matrix vector product 'AR_rspns_q' 
#################################

new_target('FORM_AR_RSPNS_Q')

depend('RSPNS_OP')

if (_lr_opt == 1):

    depend('F_OMG')

    DERIVATIVE({LABEL_RES:'F_AR_rspns_q',
                LABEL_IN:'F_OMG',
                OP_RES:'AR_rspns_q',
                OP_DERIV:['T','C0'],
                OP_MULT:['R_q','R_mu']})

elif (_lr_opt == 2):

    depend('H0')

### calculate <Psi_0|L H_bar|Phi> R_c and then diff. w.r.t. L ####

    _expand_product_basis={LABEL:'F_AR_rspns_t',
                           NEW:True,
                           OP_RES:'den12'}

    _ops_contract={OPERATORS:['C0^+','L','H','R_mu'],
                   IDX_SV:[1,2,3,4]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _expand_product_basis={LABEL:'F_AR_rspns_t',
                           NEW:False,
                           OP_RES:'den12'}

    _ops_contract={OPERATORS:['C0^+','L','H','T','R_mu'],
                   IDX_SV:[1,2,3,4,5],
                   CONNECT:[3,4]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))
    

    _ops_contract={OPERATORS:['C0^+','L','T','H','R_mu'],
                   IDX_SV:[1,2,3,4,5],
                   CONNECT:[3,4],
                   FAC:-1.0}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','H','T','T','R_mu'],
                   IDX_SV:[1,2,3,4,5,6],
                   FIX_VTX:True,
                   FAC:0.50}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','T','H','T','R_mu'],
                   IDX_SV:[1,2,3,4,5,6],
                   FIX_VTX:True,
                   FAC:-1.0}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','T','T','H','R_mu'],
                   IDX_SV:[1,2,3,4,5,6],
                   FIX_VTX:True,
                   FAC:0.50}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

### calculate <Psi_0|L[H_bar,tau]|Psi_0>R_t and then diff. w.r.t L ###

    _ops_contract={OPERATORS:['C0^+','L','H','R_q','C0'],
                   IDX_SV:[1,2,3,4,5],
                   CONNECT:[3,4]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','R_q','H','C0'],
                   IDX_SV:[1,2,3,4,5],
                   CONNECT:[3,4],
                   FAC:-1.0}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','H','T','R_q','C0'],
                   IDX_SV:[1,2,3,4,5,6]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','R_q','H','T','C0'],
                   IDX_SV:[1,2,3,4,5,6],
                   FAC:-1.0}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','T','H','R_q','C0'],
                   IDX_SV:[1,2,3,4,5,6],
                   FAC:-1.0}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    _ops_contract={OPERATORS:['C0^+','L','R_q','T','H','C0'],
                   IDX_SV:[1,2,3,4,5,6]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))


    SELECT_SPECIAL({LABEL_RES:'F_AR_rspns_t',
                    LABEL_IN:'F_AR_rspns_t',
                    TYPE:'nonzero',
                    MODE:'sum'})

    DERIVATIVE({LABEL_RES:'F_AR_rspns_q',
                LABEL_IN:'F_AR_rspns_t',
                OP_RES:'AR_rspns_q',
                OP_DERIV:'L'})



#### The next two targets has been used to form an intermediate for 'FORM_AR_RSPNS_Q'
## involving the integrals with more particle indices.   

new_target('F_prePPrint')

depend('FORM_AR_RSPNS_Q',
       'DEF_ME_INT_PP',
       'H_PP')

CLONE_OPERATOR({LABEL:'INT_PPr',
                TEMPLATE:'INT_PP'})

REPLACE({LABEL_RES:'F_prePPrint',
         LABEL_IN:'F_AR_rspns_q',
         OP_LIST:['H','H_PP']})

INVARIANT({LABEL_RES:'F_prePPrint',
           LABEL_IN:'F_prePPrint',
           OP_RES:'AR_rspns_q',
           OPERATORS:'H'})


new_target('F_PPrint')

depend('F_prePPrint')

DERIVATIVE({LABEL_RES:'F_PPrint',
           LABEL_IN:'F_prePPrint',
           OP_RES:'INT_PPr',
           OP_DERIV:'H_PP'})


#################################
# Setting up the formula for metric vector product 'AR_rspns_mu'
#################################

new_target('FORM_AR_RSPNS_MU')

depend('RSPNS_OP')

if (_lr_opt == 1):

    depend('F_OMG_C0',
           'E(MR)')
 
    DERIVATIVE({LABEL_RES:'F_AR_rspns_mu',
                LABEL_IN:'F_OMG_C0',
                OP_RES:'AR_rspns_mu',
                OP_DERIV:['T','C0'],
                OP_MULT:['R_q','R_mu']})

elif (_lr_opt == 2):

    depend('F_E_C0',
           'E(MR)')


    DERIVATIVE({LABEL_RES:'F_AR_rspns_c',
                LABEL_IN:'F_E_C0',
                OP_RES:'den12',
                OP_DERIV:'C0',
                OP_MULT:'R_mu'})

    EXPAND_OP_PRODUCT({LABEL:'F_AR_rspns_c',NEW:False,
                       OP_RES:'den12',
                       OPERATORS:['C0^+','H','R_q','C0'],
                       IDX_SV:[1,2,3,4],
                       CONNECT:[2,3]})

    EXPAND_OP_PRODUCT({LABEL:'F_AR_rspns_c',NEW:False,
                       OP_RES:'den12',
                       OPERATORS:['C0^+','H','T','R_q','C0'],
                       IDX_SV:[1,2,3,4,5]})

    SELECT_SPECIAL({LABEL_RES:'F_AR_rspns_c',
                    LABEL_IN:'F_AR_rspns_c',
                    TYPE:'nonzero',
                    MODE:'sum'})

    DERIVATIVE({LABEL_RES:'F_AR_rspns_mu',
                LABEL_IN:'F_AR_rspns_c',
                OP_RES:'AR_rspns_mu',
                OP_DERIV:'C0^+'})
    


EXPAND_OP_PRODUCT({LABEL:'F_AR_rspns_mu',
                  NEW:False,
                  OP_RES:'AR_rspns_mu',
                  OPERATORS:['AR_rspns_mu','E(MR)','R_mu','AR_rspns_mu'],
                  IDX_SV:[1,2,3,1],
                  FAC:-1.0})

##################################
#Setting up the formula for the metric vector products 'SR_rspns_q' and 'SR_rspns_mu' respectively
#################################
## The LR method involves a non unit metric matrix, it has expression with 'R_q' and 'T'. 
## For EOM it is simply an unit matrix



new_target('RSPNS_FORM')

depend('Dtr','F_T','Ttr','RSPNS_OP','C0')

if (_lr_opt == 1):

    _expand_product_basis={LABEL:'F_den12',
                           NEW:True, 
                           OP_RES:'den12'}


    _ops_contract={OPERATORS:['den12','C0^+','L','R_q','C0','den12'],
                   IDX_SV:[1,2,3,4,5,1]}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _expand_product_basis={LABEL:'F_den12',
                           NEW:False, 
                           OP_RES:'den12'}

    _ops_contract={OPERATORS:['den12','C0^+','L','R_q','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,1],
                   FAC:0.5}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','R_q','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,1],
                   FAC:-0.5}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','R_q','T','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:1.0/6}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','R_q','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:-0.33333333333333330}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','T','R_q','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,1],
                   FIX_VTX:True,
                   FAC:0.16666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','R_q','T','T','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:0.041666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','R_q','T','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:-0.125}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','T','R_q','T','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:0.125}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    _ops_contract={OPERATORS:['den12','C0^+','L','T','T','T','R_q','C0','den12'],
                   IDX_SV:[1,2,3,4,5,6,7,8,1],
                   FIX_VTX:True,
                   FAC:-0.041666666666666660}

    EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

    SELECT_SPECIAL({LABEL_RES:'F_den12',
                    LABEL_IN:'F_den12',
                    TYPE:'nonzero',
                    MODE:'sum'})

elif (_lr_opt == 2):

    EXPAND_OP_PRODUCT({LABEL:'F_den12',NEW:True,
                       OP_RES:'den12',
                       OPERATORS:['den12','C0^+','L','R_q','C0','den12'],
                       IDX_SV:[1,2,3,4,5,1]})



DERIVATIVE({LABEL_RES:'F_SR_rspns_q',
            LABEL_IN:'F_den12',
            OP_RES:'SR_rspns_q',
            OP_DERIV:'L'})

_expand_product_basis={LABEL:'F_SR_rspns_mu',
                       NEW:True, 
                       OP_RES:'SR_rspns_mu'}


_ops_contract={OPERATORS:['SR_rspns_mu','R_mu','SR_rspns_mu'],
               IDX_SV:[1,2,1]}

EXPAND_OP_PRODUCT(dict(_expand_product_basis.items()+_ops_contract.items()))

INVARIANT({LABEL_RES:'F_R_q',
           LABEL_IN:'F_T',
           OP_RES:'R_q',
           OPERATORS:'D'})

REPLACE({LABEL_RES:'F_R_q',LABEL_IN:'F_R_q',
          OP_LIST:['Ttr','R_prime_q']})



## Finish setting up formulas 
###########################
#Here starts the iterations to get excited states of all the different 
#spatial symmetry and spin multiplicity possible

_first_iter = True

for _icnt in range (0,_ncnt):

    _sym_arr = _sym[_icnt]
    _s2 =int(_mult[_icnt])

    _ms = 0
    if ((_s2 % 2) == 0):
        _ms = 1    
    if ((_ms == 0) and ((_s2 % 4) == 1)):
        _msc = 1
    elif ((_ms == 0) and ((_s2 % 4) == 3)):
        _msc = -1
    else:
        _msc = 0

    for _isym in range (0,_nsym):

        _no_root = int(_sym_arr[_isym])

        if( _no_root == 0):
            continue

        _isym_r = _multd2h[_isym,_isym_0-1]
       
        if (_s2 == _s2_0):
            s2_r = 1
        elif (abs(_s2-_s2_0) == 2):
            s2_r = 3
        else:
            quit_error('cannot handle this S2 difference')

        _ms_r = _ms - _ms_0
        _msc_r = 0
      
        if ((_ms_r == 0) and (s2_r == 1)):
                _msc_r = 1
        elif ((_ms_r == 0) and (s2_r == 3)):
                _msc_r = -1
      
        _extnsn = str(_isym+1) + '_' + str(_msc+1)
      
        _list_rspns = 'LIST_RSPNS_' + _extnsn

      
### Defining the lists for all the operators involved

        new_target(_list_rspns)
        depend('RSPNS_OP','DEF_ME_C0','DEF_ME_Dtrdag','H0','DEF_ME_T',
               'DIA_T','DIA_C0','DEF_ME_E(MR)','F_prePPrint')
        
        _op_list={'R_q':[_isym_r,_msc_r],
                  'R_mu':[_isym+1,_msc],
                  'R_prime_q':[_isym_r,0]}
        
        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':0,AB_SYM:_op_list[_op][1],MIN_REC:1,MAX_REC:_no_root})

#       DEF_ME_LIST({LIST:'ME_DIAG_t'+_extnsn,OPERATOR:'DAI_T',IRREP:_op_list[_op][0],
#                    '2MS':0,AB_SYM:_op_list[_op][1],MIN_REC:1,MAX_REC:_no_root})
#     
#       DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
#                    '2MS':0,AB_SYM:_op_list[_op][1],MIN_REC:1,MAX_REC:_no_root})

        _op_list={'AR_rspns_q':[_isym_r,_msc_r],
                  'AR_rspns_mu':[_isym+1,_msc],
                  'SR_rspns_q':[_isym_r,_msc_r],
                  'SR_rspns_mu':[_isym+1,_msc],
                  'INT_PPr':[_isym_r,_msc_r]}
        
        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':0,AB_SYM:_op_list[_op][1]})
        
        _op_list={'DIA_T':['ME_DIAG_t',_isym_r],
                  'DIA_C0':['ME_DIAG_c',_isym+1]}
        
        for _op in _op_list:
            DEF_ME_LIST({LIST:_op_list[_op][0]+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][1],
                        '2MS':0})
        
        DEF_ME_LIST({LIST:'ME_MINEN'+_extnsn,
                     OPERATOR:'E(MR)',
                     IRREP:1,
                     '2MS':0})
        
        ASSIGN_ME2OP({LIST:'ME_E(MR)',OPERATOR:'E(MR)'})
        
        _diag_cal_q = 'DIAG_CAL_q_' + _extnsn

### Setting up the preconditioner corresponding to the equations to solve 'R_q'

        new_target(_diag_cal_q)
        
        depend(_list_rspns)
        depend('EVAL_FREF','FOPT_Atr')
        
        PRECONDITIONER({LIST_PRC:'ME_DIAG_t'+_extnsn,
                        LIST_INP:'ME_FREF'})
        
        ASSIGN_ME2OP({LIST:'ME_Dtr',
                     OPERATOR:'Dtr'})
       
        EVALUATE({FORM:'FOPT_Atr'})
        
        EXTRACT_DIAG({LIST_RES:'ME_DIAG_t'+_extnsn,
                    LIST_IN:'ME_A',
                     MODE:'extend'})
        

        _diag_cal_mu = 'DIAG_CAL_mu_' + _extnsn

### Setting up the preconditioner corresponding to the equations to solve 'R_q'

        new_target(_diag_cal_mu)
      
        depend(_list_rspns)
       
        PRECONDITIONER({LIST_PRC:'ME_DIAG_c'+_extnsn,
                        LIST_INP:'H0',
                        MODE:'dia-H'})
        
        SCALE_COPY({LIST_RES:'ME_MINEN'+_extnsn,
                    LIST_INP:'ME_E(MR)',
                    FAC:-1.0})
        
        EXTRACT_DIAG({LIST_RES:'ME_DIAG_c'+_extnsn,
                     LIST_IN:'ME_MINEN'+_extnsn,
                     MODE:'ext_act'})
        
        
        _rspns_opt = 'RSPNS_OPT_' + _extnsn


### Optimising all the Formula.
        
        new_target(_rspns_opt)
        
        depend('RSPNS_FORM ',_list_rspns,'FORM_AR_RSPNS_Q','FORM_AR_RSPNS_MU',
               'DEF_ME_E(MR)','F_PPrint')
        
        OPTIMIZE({LABEL_OPT:'RSPNS_OPT'+_extnsn,
                  LABELS_IN:['F_AR_rspns_q','F_AR_rspns_mu','F_SR_rspns_q','F_SR_rspns_mu','F_R_q'],
                  INTERM:'F_PPrint'})
        

#### projecting out the elements of the ground state during the solution for the 
#### states with same symmetry as the ground state

        if (_first_iter):
            _prj_form = 'PRJ_FORM_'
            new_target(_prj_form)
            
            EXPAND_OP_PRODUCT({LABEL:'F_prj',
                               OP_RES:'R_mu',
                               OPERATORS:['R_mu','C0','C0^+','R_mu','R_mu'],
                               IDX_SV:[1,2,3,4,1],
                               AVOID:[1,4,3,5],
                               FAC:-1.0})
            
            OPTIMIZE({LABEL_OPT:'FOPT_prj',
                      LABELS_IN:'F_prj'})
        
   
            _first_iter = False

#### Finally solving the eigen value equation to get the excitation energies
        
        _solve_eqn = 'SOLVE_EQN_' + _extnsn

        new_target(_solve_eqn)
        depend('INPUT_INFO')
        depend(_rspns_opt)
        depend(_diag_cal_q)
        depend(_diag_cal_mu)
        depend(_prj_form)
        
        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_R_q'+_extnsn,'ME_R_mu'+_extnsn]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_t'+_extnsn,'ME_DIAG_c'+_extnsn] 
        _solve_evp_basis[OP_MVP]=['AR_rspns_q','AR_rspns_mu'] 
        _solve_evp_basis[OP_SVP]=['SR_rspns_q','SR_rspns_mu'] 
        _solve_evp_basis[FORM]='RSPNS_OPT'+_extnsn
        _solve_evp_basis[LIST_SPC]=['ME_R_prime_q'+_extnsn,'ME_Dtr','ME_Dtrdag']
        _solve_evp_basis[MODE]='TRF PRJ'
        _solve_evp_basis[FORM_SPC]='FOPT_prj'
        _solve_evp_basis[N_ROOTS]=_no_root
        _solve_evp_basis[CHOICE_OPT]=_choice

        PRINT({STRING: 'Doing calculation of irrep:    ' + str(_isym+1) + 
                       '  and of spin multiplicity:    ' + str(_s2)})

        PRINT({STRING: 'isym_r, msc_r:' + str(_isym_r) + ',  ' + str(_msc_r)})
        
        SOLVE_EVP(_solve_evp_basis)
 

        _get_mel_inf = 'GET_MEL_INF' + _extnsn 

        new_target(_get_mel_inf,True)
        
        depend(_solve_eqn)

        DEF_ME_LIST({LIST:'ME_SR_q'+_extnsn,
                     OPERATOR:'SR_rspns_q',
                     IRREP:_isym_r,
                     '2MS':0,
                     AB_SYM:_msc_r,
                     MIN_REC:1,MAX_REC:_no_root})


        OPTIMIZE({LABEL_OPT:'F_OPT_SR'+_extnsn,
                  LABELS_IN:'F_SR_rspns_q'})

        TRANSFORM({LIST_IN:'ME_R_q'+_extnsn,
                   LIST_OUT:'ME_SR_q'+_extnsn,
                   FORM:'F_OPT_SR'+_extnsn})
       
        ANALYZE_MEL({LISTS:['ME_R_mu'+_extnsn,'ME_R_q'+_extnsn],
                     LISTS_CV:['ME_R_mu'+_extnsn,'ME_SR_q'+_extnsn]})

        PRINT({STRING: 'Done calculation of irrep:    ' + str(_isym+1) + 
                       '  and of spin multiplicity:    ' + str(_s2)})

export_targets();
