#This is a python interface for GeCCo to solve the first order L(\omega) to 
#calculate static and dynamic polarizability using the linear response method.
#This works when any parameters related to CI coeffients are not relaxed

from gecco_interface import *
from get_response_data import _response_data
import math

_inp=GeCCo_Input()

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

    new_target('L(1)'+i_par)
    depend('DEF_ME_T')
    depend('DEF_ME_L')
    depend('DEF_ME_T(1)'+i_par)

    CLONE_OPERATOR({LABEL:'L(1)'+i_par,
                    TEMPLATE:'L'})

#Setting the order of the T(1) operator to be first

    SET_ORDER({LABEL:'L(1)'+i_par,
              ORDER:1,
              SPECIES:1})

    new_target('RSPNS(2)_OP'+i_par)

    depend('OMG')
    depend('T(1)'+i_par)

    CLONE_OPERATOR({LABEL:'Ltr(1)'+i_par,TEMPLATE:'L'})
    CLONE_OPERATOR({LABEL:'O(1)LL'+i_par,TEMPLATE:'OMG',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'O(1)RL'+i_par,TEMPLATE:'OMG',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'DIA_L(1)'+i_par,TEMPLATE:'L'})

    DEF_SCALAR({LABEL:'RED_LAG(1)LL'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RL'+i_par})
    DEF_SCALAR({LABEL:'preRED_LAG(1)RL_1'+i_par})
    DEF_SCALAR({LABEL:'preRED_LAG(1)RL_2'+i_par})


    new_target('F_O(1)LL'+i_par)

    depend('F_MRCC_LAG')
    depend('L(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LL'+i_par,
                LABEL_IN:'F_MRCC_LAG',
                OP_RES:'RED_LAG(1)LL'+i_par,
                OP_DERIV:'L',
                OP_MULT:'L(1)'+i_par})

    DERIVATIVE({LABEL_RES:'F_O(1)LL'+i_par,
                LABEL_IN:'F_RED_LAG(1)LL'+i_par,
                OP_RES:'O(1)LL'+i_par,
                OP_DERIV:'T'})

    PRINT_FORMULA({LABEL:'F_O(1)LL'+i_par})

    new_target('F_O(1)RL'+i_par)

    depend('L(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)
    depend('T(1)'+i_par)
    depend('IMPORT_V(1)')
    depend('F_MRCC_LAG')
    
    DERIVATIVE({LABEL_RES:'F_preRED_LAG(1)RL_1'+i_par,
                LABEL_IN:'F_MRCC_LAG',
                OP_RES:'preRED_LAG(1)RL_1'+i_par,
                OP_DERIV:'T',
                OP_MULT:'T(1)'+i_par}) 

    REPLACE({LABEL_RES:'F_preRED_LAG(1)RL_2'+i_par,
             LABEL_IN:'F_MRCC_LAG',
             OP_LIST:['H','V(1)']})
    
    INVARIANT({LABEL_RES:'F_preRED_LAG(1)RL_2'+i_par,
               LABEL_IN:'F_preRED_LAG(1)RL_2'+i_par,
               OP_RES:'preRED_LAG(1)RL_2'+i_par,
               OPERATORS:'H'})
    
    DEF_FORMULA({LABEL:'F_RED_LAG(1)RL'+i_par,
                FORMULA:'RED_LAG(1)RL'+i_par+'=preRED_LAG(1)RL_1'+i_par+'+preRED_LAG(1)RL_2'+i_par})    

    EXPAND({LABEL_RES:'F_RED_LAG(1)RL'+i_par,
            LABEL_IN:'F_RED_LAG(1)RL'+i_par,
            INTERM:['F_preRED_LAG(1)RL_1'+i_par,'F_preRED_LAG(1)RL_2'+i_par]}) 

    DERIVATIVE({LABEL_RES:'F_O(1)RL'+i_par,
                LABEL_IN:'F_RED_LAG(1)RL'+i_par,
                OP_RES:'O(1)RL'+i_par,
                OP_DERIV:'T'})

    PRINT_FORMULA({LABEL:'F_O(1)RL'+i_par})

    new_target('F_L(1)'+i_par)

    depend('F_L')

    INVARIANT({LABEL_RES:'F_L(1)'+i_par,
               LABEL_IN:'F_L',
               OP_RES:'L(1)'+i_par,
               OPERATORS:'H'})

    REPLACE({LABEL_RES:'F_L(1)'+i_par,LABEL_IN:'F_L(1)'+i_par,
              OP_LIST:['Ltr','Ltr(1)'+i_par]})
    
    new_target('DEF_ME_L(1)'+i_par)

    depend('L(1)'+i_par)

    _op_list={'L(1)'+i_par:'ME_L(1)'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    new_target('LIST_RSPNS(2)'+i_par) 

    depend('RSPNS(2)_OP'+i_par)

    _op_list={'O(1)LL'+i_par:'ME_O(1)LL'+i_par,
              'O(1)RL'+i_par:'ME_O(1)RL'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:1})

    _op_list={'Ltr(1)'+i_par:'ME_Ltr(1)'+i_par,
              'DIA_L(1)'+i_par:'ME_DIA_L(1)'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0})

    new_target('DIAG_L(1)'+i_par)

    depend('EVAL_FREF','FOPT_Atr','DIAG1SxxM00_T','EVAL_Atr')

    COPY_LIST({LIST_RES:'ME_DIA_L(1)'+i_par,LIST_INP:'DIAG1SxxM00_T',ADJOINT:True,FAC:1.0})

    new_target('OPT_LAMBDA(1)'+i_par)

    depend('DEF_ME_L(1)'+i_par)
    depend('LIST_RSPNS(2)'+i_par)
    depend('F_O(1)LL'+i_par)
    depend('F_O(1)RL'+i_par)
    depend('F_L(1)'+i_par)
    depend('DEF_ME_Ttr')
    depend('DIAG_L(1)'+i_par)

    OPTIMIZE({LABEL_OPT:'FOPT_L(1)'+i_par,
              LABELS_IN:'F_L(1)'+i_par})

    OPTIMIZE({LABEL_OPT:'FOPT_SOLVE_L(1)'+i_par,
              LABELS_IN:['F_O(1)LL'+i_par,'F_O(1)RL'+i_par]})

    new_target('SOLVE_LAMBDA(1)'+i_par,True)

    depend('OPT_LAMBDA(1)'+i_par)
    depend('DIAG_L(1)'+i_par)
    depend('DEF_ME_Dtrdag','DEF_ME_Ttr')

    SOLVE_LEQ({LIST_OPT:'ME_L(1)'+i_par,
               LIST_PRC:'ME_DIA_L(1)'+i_par,
               OP_MVP:'O(1)LL'+i_par,
               OP_SVP:'L(1)'+i_par,
               OP_RHS:'O(1)RL'+i_par,
               FORM:'FOPT_SOLVE_L(1)'+i_par,
               LIST_SPC:['ME_L(1)'+i_par,'ME_Ltr(1)'+i_par,'ME_Dtr','ME_Dtrdag'],
               FORM_SPC:'FOPT_L(1)'+i_par,
               MODE:'TRF',
               N_ROOTS:1})
