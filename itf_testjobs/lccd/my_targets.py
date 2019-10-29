from python_interface.gecco_interface import *
import python_interface.gecco_modules.string_to_form as stf

# Operators ------------------------------------------------------
new_target('CCD_OPS')
DEF_SCALAR({LABEL:'LCCD'})
DEF_SCALAR({LABEL:'ECCD'})
DEF_EXCITATION({LABEL:'T2',MIN_RANK:2,MAX_RANK:2})
CLONE_OPERATOR({LABEL:'O2',TEMPLATE:'T2'})
DEF_OP_FROM_OCC({LABEL:'INT_PP',DESCR:"PP,PP"})


# Formula --------------------------------------------------------
new_target('CCD_FORM')
depend('H')
depend('CCD_OPS')
form_l1 = stf.Formula("CCD_FORM:LCCD=<H>")
form_l1.append("<[H,T2]>")
form_l1.append("<T2^+*H>")
form_l1.append("<T2^+*[H,T2]>")

form_l1.set_rule()

PRINT_FORMULA({LABEL:'CCD_FORM',MODE:'SHORT'})

EXPAND_OP_PRODUCT({LABEL:'FORM_PP',OP_RES:'INT_PP',
                   OPERATORS:['H'],
                   IDX_SV:[1],
                   LABEL_DESCR:["1,,PP,PP"]})

# Residual
new_target('CCD_RES')
depend('CCD_FORM')

#FACTOR_OUT({LABEL_IN:'CCD_FORM',
#            LABEL_RES:'CCD_FORM',
#            INTERM:'H'})
DERIVATIVE({LABEL_RES:'CCD_RES',LABEL_IN:'CCD_FORM',
            OP_RES:'O2',OP_DERIV:'T2^+'})
PRINT_FORMULA({LABEL:'CCD_RES',MODE:'SHORT'})

# Energy
new_target('CCD_EN')
depend('CCD_FORM')
form_e1 = stf.Formula("CCD_EN:ECCD=<H>")
form_e1.append("<[H,T2]>")
form_e1.set_rule()
PRINT_FORMULA({LABEL:'CCD_EN',MODE:'SHORT'})


# ME lists ------------------------------------------------------
new_target('CCD_LISTS')
depend('CCD_OPS')
DEF_ME_LIST({LIST:'ME_LCCD',OPERATOR:'LCCD',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_ECCD',OPERATOR:'ECCD',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_T2',OPERATOR:'T2',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_O2',OPERATOR:'O2',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_PP',OPERATOR:'INT_PP',IRREP:1,'2MS':0,AB_SYM:+1})


# Optimise formula -----------------------------------------------
new_target('CCD_OPT')
depend('CCD_RES','CCD_EN','CCD_LISTS','H0')
OPTIMIZE({LABEL_OPT:'CCD_OPT',LABELS_IN:['CCD_RES','CCD_EN']})

TRANSLATE_ITF({LABEL: 'CCD_OPT',
               OUTPUT:'icmrcc_lccd.itfaa',
               TITLE: 'icmrcc_lccd.formulae',
               MULTI:False,
               ITIN:True})
PRINT_FORMULA({LABEL:'CCD_OPT',MODE:'SHORT'})


# Set up preconditioner -----------------------------------------
new_target('DIAG')
depend('H0','CCD_OPS')
CLONE_OPERATOR({LABEL:'D2',TEMPLATE:'T2'})
DEF_ME_LIST({LIST:'DIAG',OPERATOR:'D2',IRREP:1,'2MS':0,AB_SYM:+1})
PRECONDITIONER({LIST_PRC:'DIAG',LIST_INP:'H0'})


# Solve equations -----------------------------------------------
new_target('CCD_METHOD')
required()
depend('DIAG')
depend('CCD_OPT')

SOLVE_NLEQ({LIST_OPT:'ME_T2',LIST_RESID:'ME_O2',LIST_PRC:'DIAG',
            MODE:'DIA',LIST_E:'ME_ECCD',FORM:'CCD_OPT'})

export_targets();
