#This is a python interface for GeCCo to solve the first order L(\omega) and \bar{c}(\omega) to 
#calculate static and dynamic polarizability using the linear response method.

from gecco_interface import *
from get_response_data import _response_data
import math

_inp=GeCCo_Input()
_orb = Orb_Info()

# Get all the necessary options needed during the calculation

_s2 = _orb.get('imult')

_ms = _orb.get('ims')

_isym = _orb.get('lsym')

_root = _inp.get('method.MR.ciroot')

_option=_response_data['option']


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

    new_target('C0_bar(1)'+i_par)
    depend('DEF_ME_C0')
    depend('DEF_ME_C0_bar')
    depend('DEF_ME_C0(1)'+i_par)

    CLONE_OPERATOR({LABEL:'C0_bar(1)'+i_par,
                    TEMPLATE:'C0_bar'})

#Setting the order of the T(1) operator to be first

    SET_ORDER({LABEL:'C0_bar(1)'+i_par,
              ORDER:1,
              SPECIES:1})

    new_target('RSPNS(2)_OP'+i_par)

    depend('OMG')
    depend('T(1)'+i_par)
    depend('A_C0')

    CLONE_OPERATOR({LABEL:'Ltr(1)'+i_par,TEMPLATE:'L'})
    CLONE_OPERATOR({LABEL:'O(1)LL'+i_par,TEMPLATE:'OMG',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'O(1)RL'+i_par,TEMPLATE:'OMG',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'O(1)LCb'+i_par,TEMPLATE:'A_C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'PRJ_O(1)LCb'+i_par,TEMPLATE:'A_C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'O(1)RCb'+i_par,TEMPLATE:'A_C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'PRJ_O(1)RCb'+i_par,TEMPLATE:'A_C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'DIA_L(1)'+i_par,TEMPLATE:'L'})
    CLONE_OPERATOR({LABEL:'DIAG_C0_bar(1)'+i_par,TEMPLATE:'C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'C0_x'+i_par,TEMPLATE:'C0',ADJOINT:True})
    CLONE_OPERATOR({LABEL:'T(1)_x'+i_par,TEMPLATE:'T',ADJOINT:True})

    DEF_SCALAR({LABEL:'RED_LAG(1)LL'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LCb'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RL'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb'+i_par})

    DEF_SCALAR({LABEL:'RED_LAG(1)LCb_1'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LCb_2'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)LCb_3'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RL_1'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RL_2'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RL_3'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb_1'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb_2'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb_3'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb_4'+i_par})
    DEF_SCALAR({LABEL:'RED_LAG(1)RCb_5'+i_par})


    new_target('F_O(1)LL'+i_par)

    depend('F_MRCC_LAG_PROP')
    depend('L(1)'+i_par)
    depend('C0_bar(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)
    depend('E(MR)')

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LL'+i_par,
                LABEL_IN:'F_MRCC_LAG_PROP',
                OP_RES:'RED_LAG(1)LL'+i_par,
                OP_DERIV:['L','C0_bar'],
                OP_MULT:['L(1)'+i_par,'C0_bar(1)'+i_par]})

    DERIVATIVE({LABEL_RES:'F_O(1)LL'+i_par,
                LABEL_IN:'F_RED_LAG(1)LL'+i_par,
                OP_RES:'O(1)LL'+i_par,
                OP_DERIV:'T'})

#   PRINT_FORMULA({LABEL:'F_O(1)LL'+i_par})

    new_target('F_O(1)LCb'+i_par)

    depend('F_MRCC_LAG_PROP')
    depend('L(1)'+i_par)
    depend('C0_bar(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LCb_1'+i_par,
                LABEL_IN:'F_MRCC_LAG_PROP',
                OP_RES:'RED_LAG(1)LCb_1'+i_par,
                OP_DERIV:'L',
                OP_MULT:'L(1)'+i_par})

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)LCb_2'+i_par,
                LABEL_IN:'F_MRCC_LAG_PROP',
                OP_RES:'RED_LAG(1)LCb_2'+i_par,
                OP_DERIV:'C0_bar',
                OP_MULT:'C0_bar(1)'+i_par})

    INVARIANT({LABEL_RES:'F_RED_LAG(1)LCb_3'+i_par,
               LABEL_IN:'F_RED_LAG(1)LCb_1'+i_par,
               OP_RES:'RED_LAG(1)LCb_3'+i_par,
               OPERATORS:'C0_bar'})

    TRANSPS_FORMULA({LABEL_IN:'F_RED_LAG(1)LCb_3'+i_par,
                     LABEL_RES:'F_RED_LAG(1)LCb_3'+i_par,
                     OP_RES:'RED_LAG(1)LCb_3'+i_par,
                     INIT:False,MULTI:True})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)LCb_3'+i_par})

    _arg_formula='RED_LAG(1)LCb'+i_par+'=RED_LAG(1)LCb_1'+i_par+'+RED_LAG(1)LCb_2'+i_par+'+RED_LAG(1)LCb_3'+i_par

    DEF_FORMULA({LABEL:'F_RED_LAG(1)LCb'+i_par,
                 FORMULA:_arg_formula})

    EXPAND({LABEL_RES:'F_RED_LAG(1)LCb'+i_par,
            LABEL_IN:'F_RED_LAG(1)LCb'+i_par,
            INTERM:['F_RED_LAG(1)LCb_1'+i_par,'F_RED_LAG(1)LCb_2'+i_par,'F_RED_LAG(1)LCb_3'+i_par]})

# substracting E <\bar{Psi_0}|Psi_0> from F_LAG_2

    EXPAND_OP_PRODUCT({LABEL:'F_RED_LAG(1)LCb'+i_par,NEW:False,
                       OP_RES:'RED_LAG(1)LCb'+i_par,
                       OPERATORS:['E(MR)','C0_bar(1)'+i_par,'C0'],
                       IDX_SV:[1,2,3],
                       FAC:-1.0})

    DERIVATIVE({LABEL_RES:'F_O(1)LCb'+i_par,
                LABEL_IN:'F_RED_LAG(1)LCb'+i_par,
                OP_RES:'O(1)LCb'+i_par,
                OP_DERIV:'C0'})

#   PRINT_FORMULA({LABEL:'F_O(1)LCb'+i_par})

    new_target('F_PRJ_O(1)LCb'+i_par)

    depend('F_O(1)LCb'+i_par)

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_O(1)LCb'+i_par,
                       OP_RES:'PRJ_O(1)LCb'+i_par,
                       NEW:True,
                       OPERATORS:['PRJ_O(1)LCb'+i_par,'O(1)LCb'+i_par,'PRJ_O(1)LCb'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_prePRJ_O(1)LCb'+i_par,
                       OP_RES:'PRJ_O(1)LCb'+i_par,
                       NEW:False,
#                      OPERATORS:['PRJ_O(1)LCb'+i_par,'C0_bar','C0_bar^+','O(1)LCb'+i_par,'PRJ_O(1)LCb'+i_par],
                       OPERATORS:['PRJ_O(1)LCb'+i_par,'O(1)LCb'+i_par,'C0_bar^+','C0_bar','PRJ_O(1)LCb'+i_par],
                       IDX_SV:[1,2,3,4,1],
#                      AVOID:[1,4,3,5],
                       AVOID:[1,3,2,5],
                       FAC:-1.0})

#   PRINT_FORMULA({LABEL:'F_prePRJ_O(1)LCb'+i_par})

    EXPAND({LABEL_RES:'F_PRJ_O(1)LCb'+i_par,
            LABEL_IN:'F_prePRJ_O(1)LCb'+i_par,
            INTERM:'F_O(1)LCb'+i_par})

#   PRINT_FORMULA({LABEL:'F_PRJ_O(1)LCb'+i_par})

    new_target('F_O(1)RL'+i_par)

    depend('L(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)
    depend('T(1)'+i_par)
    depend('IMPORT_V(1)')
    depend('F_MRCC_LAG_PROP')
    depend('F_preRSPNS(2)')
    
    INVARIANT({LABEL_RES:'F_RED_LAG(1)RL_1'+i_par,
               LABEL_IN:'F_preRSPNS(2)',
               OP_RES:'RED_LAG(1)RL_1'+i_par,
               OPERATORS:'L(1)'+i_par})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RL_1'+i_par})

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)RL_2'+i_par,
                LABEL_IN:'F_MRCC_LAG_PROP',
                OP_RES:'RED_LAG(1)RL_2'+i_par,
                OP_DERIV:'T',
                OP_MULT:'T(1)'+i_par})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RL_2'+i_par})

    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)RL_3'+i_par,
                LABEL_IN:'F_MRCC_LAG_PROP',
                OP_RES:'RED_LAG(1)RL_3'+i_par,
                OP_DERIV:['C0','C0^+'],
                OP_MULT:['C0(1)'+i_par,'C0_x'+i_par]})

    REPLACE({LABEL_RES:'F_RED_LAG(1)RL_3'+i_par,
             LABEL_IN:'F_RED_LAG(1)RL_3'+i_par,
             OP_LIST:['C0_x'+i_par, 'C0(1)'+i_par+'^+']})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RL_3'+i_par})

    _arg_formula='RED_LAG(1)RL'+i_par+'=RED_LAG(1)RL_1'+i_par+'+RED_LAG(1)RL_2'+i_par+'+RED_LAG(1)RL_3'+i_par
    DEF_FORMULA({LABEL:'F_RED_LAG(1)RL'+i_par,
                FORMULA:_arg_formula})

    EXPAND({LABEL_RES:'F_RED_LAG(1)RL'+i_par,
            LABEL_IN:'F_RED_LAG(1)RL'+i_par,
            INTERM:['F_RED_LAG(1)RL_1'+i_par,'F_RED_LAG(1)RL_2'+i_par,'F_RED_LAG(1)RL_3'+i_par]})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RL'+i_par})

    DERIVATIVE({LABEL_RES:'F_O(1)RL'+i_par,
                LABEL_IN:'F_RED_LAG(1)RL'+i_par,
                OP_RES:'O(1)RL'+i_par,
                OP_DERIV:'T'})

#   PRINT_FORMULA({LABEL:'F_O(1)RL'+i_par})

    new_target('F_O(1)RCb'+i_par)

    depend('L(1)'+i_par)
    depend('RSPNS(2)_OP'+i_par)
    depend('T')
    depend('T(1)'+i_par)
    depend('IMPORT_V(1)')
    depend('F_MRCC_LAG_PROP')
    depend('F_preRSPNS(2)')
    depend('F_O(1)LCb'+i_par)

    INVARIANT({LABEL_RES:'F_RED_LAG(1)RCb_1'+i_par,
               LABEL_IN:'F_RED_LAG(1)LCb'+i_par,
               OP_RES:'RED_LAG(1)RCb_1'+i_par,
               OPERATORS:'E(MR)'})

    REPLACE({LABEL_RES:'F_RED_LAG(1)RCb_1'+i_par,
             LABEL_IN:'F_RED_LAG(1)RCb_1'+i_par,
             OP_LIST:['L(1)'+i_par,'L','L(1)'+i_par+'^+','L^+','C0_bar(1)'+i_par,'C0_bar']})


    DERIVATIVE({LABEL_RES:'F_RED_LAG(1)RCb_2'+i_par,
                LABEL_IN:'F_RED_LAG(1)RCb_1'+i_par,
                OP_RES:'RED_LAG(1)RCb_2'+i_par,
                OP_DERIV:['T','T^+'],
                OP_MULT:['T(1)'+i_par,'T(1)_x'+i_par]})

    REPLACE({LABEL_RES:'F_RED_LAG(1)RCb_2'+i_par,
             LABEL_IN:'F_RED_LAG(1)RCb_2'+i_par,
             OP_LIST:['T(1)_x'+i_par,'T(1)'+i_par+'^+']})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RCb_2'+i_par})

    REPLACE({LABEL_RES:'F_RED_LAG(1)RCb_3'+i_par,
             LABEL_IN:'F_RED_LAG(1)RCb_1'+i_par,
             OP_LIST:['H','V(1)']})

    INVARIANT({LABEL_RES:'F_RED_LAG(1)RCb_4'+i_par,
               LABEL_IN:'F_RED_LAG(1)RCb_3'+i_par,
               OP_RES:'RED_LAG(1)RCb_4'+i_par,
               OPERATORS:'H'})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RCb_4'+i_par})

    INVARIANT({LABEL_RES:'F_RED_LAG(1)RCb_5'+i_par,
               LABEL_IN:'F_RED_LAG(1)RCb_1'+i_par,
               OP_RES:'RED_LAG(1)RCb_5'+i_par,
               OPERATORS:'C0_bar'})

    REPLACE({LABEL_RES:'F_RED_LAG(1)RCb_5'+i_par,
             LABEL_IN:'F_RED_LAG(1)RCb_5'+i_par,
             OP_LIST:['C0^+','C0(1)'+i_par+'^+']})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RCb_5'+i_par})

    _arg_formula='RED_LAG(1)RCb'+i_par+'=RED_LAG(1)RCb_2'+i_par+'+RED_LAG(1)RCb_4'+i_par+'+RED_LAG(1)RCb_5'+i_par
    DEF_FORMULA({LABEL:'F_RED_LAG(1)RCb'+i_par,
                FORMULA:_arg_formula})

    EXPAND({LABEL_RES:'F_RED_LAG(1)RCb'+i_par,
            LABEL_IN:'F_RED_LAG(1)RCb'+i_par,
            INTERM:['F_RED_LAG(1)RCb_2'+i_par,'F_RED_LAG(1)RCb_4'+i_par,'F_RED_LAG(1)RCb_5'+i_par]})

#   PRINT_FORMULA({LABEL:'F_RED_LAG(1)RCb'+i_par})

    DERIVATIVE({LABEL_RES:'F_O(1)RCb'+i_par,
                LABEL_IN:'F_RED_LAG(1)RCb'+i_par,
                OP_RES:'O(1)RCb'+i_par,
                OP_DERIV:'C0'})

#   PRINT_FORMULA({LABEL:'F_O(1)RCb'+i_par})

    new_target('F_PRJ_O(1)RCb'+i_par)

    depend('F_O(1)RCb'+i_par)

    EXPAND_OP_PRODUCT({LABEL:'F_PRJ_O(1)RCb'+i_par,
                       OP_RES:'PRJ_O(1)RCb'+i_par,
                       NEW:True,
                       OPERATORS:['PRJ_O(1)RCb'+i_par,'O(1)RCb'+i_par,'PRJ_O(1)RCb'+i_par],
                       IDX_SV:[1,2,1],
                       FAC:1.0})

    EXPAND_OP_PRODUCT({LABEL:'F_PRJ_O(1)RCb'+i_par,
                       OP_RES:'PRJ_O(1)RCb'+i_par,
                       NEW:False,
#                      OPERATORS:['PRJ_O(1)LCb'+i_par,'C0_bar','C0_bar^+','O(1)LCb'+i_par,'PRJ_O(1)LCb'+i_par],
                       OPERATORS:['PRJ_O(1)RCb'+i_par,'O(1)RCb'+i_par,'C0_bar^+','C0_bar','PRJ_O(1)RCb'+i_par],
                       IDX_SV:[1,2,3,4,1],
#                      AVOID:[1,4,3,5],
                       AVOID:[1,3,2,5],
                       FAC:-1.0})

#   PRINT_FORMULA({LABEL:'F_prePRJ_O(1)RCb'+i_par})

#   EXPAND({LABEL_RES:'F_PRJ_O(1)RCb'+i_par,
#           LABEL_IN:'F_prePRJ_O(1)RCb'+i_par,
#           INTERM:'F_O(1)RCb'+i_par})

#   PRINT_FORMULA({LABEL:'F_PRJ_O(1)RCb'+i_par})

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

#creating ME list for C0_bar(1):
    new_target('DEF_ME_C0_bar(1)'+i_par)

    depend('C0_bar(1)'+i_par)

    _op_list={'C0_bar(1)'+i_par:'ME_C0_bar(1)'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_isym,
                     '2MS':0,
                     AB_SYM:_msc})

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

    _op_list={'PRJ_O(1)LCb'+i_par:'ME_PRJ_O(1)LCb'+i_par,
              'PRJ_O(1)RCb'+i_par:'ME_PRJ_O(1)RCb'+i_par,
              'O(1)LCb'+i_par:'ME_O(1)LCb'+i_par,
              'O(1)RCb'+i_par:'ME_O(1)RCb'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_isym,
                     '2MS':0,
                     AB_SYM:_msc})

    _op_list={'Ltr(1)'+i_par:'ME_Ltr(1)'+i_par,
              'DIA_L(1)'+i_par:'ME_DIA_L(1)'+i_par}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:1,
                     '2MS':0})

    _op_list={'DIAG_C0_bar(1)'+i_par:'ME_DIAG_C0_bar(1)'+i_par}
    
    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_isym,
                     '2MS':_ms})

    new_target('DIAG_L(1)'+i_par)

    depend('EVAL_FREF','FOPT_Atr','DIAG1SxxM00_T','EVAL_Atr')

    COPY_LIST({LIST_RES:'ME_DIA_L(1)'+i_par,LIST_INP:'DIAG1SxxM00_T',ADJOINT:True,FAC:1.0})

    new_target('DIAG_C0_bar(1)'+i_par)

    depend('DIAG_C0(1)'+i_par)

    COPY_LIST({LIST_RES:'ME_DIAG_C0_bar(1)'+i_par,LIST_INP:'ME_DIAG_C0(1)'+i_par,ADJOINT:True,FAC:1.0})

    new_target('OPT_LAMBDA(1)'+i_par)

    depend('H0','DEF_ME_E(MR)')
    depend('DEF_ME_L(1)'+i_par)
    depend('LIST_RSPNS(2)'+i_par)
    depend('F_O(1)LL'+i_par)
    depend('F_PRJ_O(1)LCb'+i_par)
    depend('F_O(1)RL'+i_par)
    depend('F_PRJ_O(1)RCb'+i_par)
    depend('F_L(1)'+i_par)
    depend('DEF_ME_Ttr')
    depend('DIAG_L(1)'+i_par)
    depend('DEF_ME_C0_bar(1)'+i_par)

    OPTIMIZE({LABEL_OPT:'FOPT_L(1)'+i_par,
              LABELS_IN:'F_L(1)'+i_par})

    OPTIMIZE({LABEL_OPT:'FOPT_O(1)RCb'+i_par,
              LABELS_IN:'F_O(1)RCb'+i_par})

    OPTIMIZE({LABEL_OPT:'FOPT_SOLVE_L(1)'+i_par,
              LABELS_IN:['F_O(1)LL'+i_par,'F_O(1)RL'+i_par,'F_PRJ_O(1)LCb'+i_par,'F_PRJ_O(1)RCb'+i_par]})
#             LABELS_IN:['F_O(1)LL'+i_par,'F_O(1)RL'+i_par,'F_O(1)LCb'+i_par,'F_O(1)RCb'+i_par]})

    new_target('SOLVE_LAMBDA(1)'+i_par,True)

    depend('OPT_LAMBDA(1)'+i_par)
    depend('DIAG_L(1)'+i_par)
    depend('DIAG_C0_bar(1)'+i_par)
    depend('DEF_ME_Dtrdag','DEF_ME_Ttr')

    EVALUATE({FORM:'FOPT_O(1)RCb'+i_par})

    SOLVE_LEQ({LIST_OPT:['ME_L(1)'+i_par,'ME_C0_bar(1)'+i_par],
               LIST_PRC:['ME_DIA_L(1)'+i_par,'ME_DIAG_C0_bar(1)'+i_par],
               OP_MVP:['O(1)LL'+i_par,'PRJ_O(1)LCb'+i_par],
#              OP_MVP:['O(1)LL'+i_par,'O(1)LCb'+i_par],
               OP_SVP:['L(1)'+i_par,'C0_bar(1)'+i_par],
               OP_RHS:['O(1)RL'+i_par,'PRJ_O(1)RCb'+i_par],
#              OP_RHS:['O(1)RL'+i_par,'O(1)RCb'+i_par],
               FORM:'FOPT_SOLVE_L(1)'+i_par,
               LIST_SPC:['ME_L(1)'+i_par,'ME_Ltr(1)'+i_par,'ME_Dtr','ME_Dtrdag'],
               FORM_SPC:'FOPT_L(1)'+i_par,
               MODE:'TRF DIA',
               N_ROOTS:1})

    PRINT_MEL({LIST:'ME_C0_bar(1)'+i_par})

