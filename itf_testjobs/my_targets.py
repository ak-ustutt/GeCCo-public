from python_interface.gecco_interface import *
import python_interface.gecco_modules.string_to_form as stf

# Operators ------------------------------------------------------
new_target('CCD_OPS')
DEF_SCALAR({LABEL:'LCCD'})
DEF_SCALAR({LABEL:'ECCD'})
DEF_EXCITATION({LABEL:'T2',MIN_RANK:2,MAX_RANK:2})
DEF_EXCITATION({LABEL:'T1',MIN_RANK:1,MAX_RANK:1})
CLONE_OPERATOR({LABEL:'O2',TEMPLATE:'T2'})
CLONE_OPERATOR({LABEL:'O1',TEMPLATE:'T1'})


# Formula --------------------------------------------------------
new_target('CCD_FORM')
depend('H')
depend('CCD_OPS')
form_l1 = stf.Formula("CCD_FORM:LCCD=<H>")
form_l1.append("<[H,T1]>")
form_l1.append("<[H,T2]>")

form_l1.append("<T1^+*H>")
form_l1.append("<T2^+*H>")

form_l1.append("<T1^+*[H,T1]>")
form_l1.append("<T1^+*[H,T2]>")
form_l1.append("<T1^+*(1/2)*[[H,T1],T1]>")
form_l1.append("<T1^+[[H,T2],T1]>")
form_l1.append("<T1^+*(1/6)*[[[H,T1],T1],T1]>")

form_l1.append("<T2^+*[H,T1]>")
form_l1.append("<T2^+*[H,T2]>")
form_l1.append("<T2^+*(1/2)*[[H,T2],T2]>")
# Above works

# Trouble with one term...
#form_l1.append("<T2^+*(1/2)*[[H,T1],T1]>")

#form_l1.append("<T2^+*(1/2)*[[H,T2],T1]>")
#form_l1.append("<T2^+*(1/6)*[[[H,T1],T1],T1]>")
#form_l1.append("<T2^+*(1/6)*[[[H,T2],T1],T1]>")

form_l1.set_rule()

# new terms
# [[H,T1],T1]
EXPAND_OP_PRODUCT({LABEL:'CCD_FORM',NEW:False,OP_RES:'LCCD',FAC:1.0,
                   OPERATORS:['T2^+','H','T1','T1'],
                   IDX_SV   :[1,2,3],
                   LABEL_DESCR:["1,2,HH,","1,3,,P","2,3,H,","1,4,,P","2,4,H,"]})

# This one doesn't work
#EXPAND_OP_PRODUCT({LABEL:'CCD_FORM',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T1','T1'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,2,H,P","1,3,H,","2,3,,P","1,4,,P","2,4,H,"]})

# 3ext
#EXPAND_OP_PRODUCT({LABEL:'CCD_FORM',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T1','T1'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,2,,PP","1,3,H,","2,3,,P","1,4,H,","2,4,,P"]})

# [[H,T2],T1]
# 3 ext
#EXPAND_OP_PRODUCT({LABEL:'CCD_FORM',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T2','T1'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,3,HH,P","2,3,,P","1,4,,P","2,4,H,"]})

#EXPAND_OP_PRODUCT({LABEL:'CCD_FORM',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T2','T1'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,3,H,PP","2,3,H,","1,4,H,","2,4,,P"]})




# 3 external integrals
EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:True,OP_RES:'LCCD',FAC:1.0,
                   OPERATORS:['T1^+','H','T1','T1'],
                   IDX_SV   :[1,2,3,4],
                   LABEL_DESCR:["1,2,,P","1,3,H,","2,3,,P","2,4,H,P"]})

#### Test term involving K and J in T2+ [H,T2]
###EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:1.0,
###                   OPERATORS:['T2^+','H','T2'],
###                   IDX_SV   :[1,2,3],
###                   LABEL_DESCR:["1,2,H,P","1,3,H,P","2,3,H,P"]})

#### Test term involving 3 internal integral in T1+ [H,T2]
###EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:1.0,
###                   OPERATORS:['T1^+','H','T2'],
###                   IDX_SV   :[1,2,3],
###                   LABEL_DESCR:["1,2,H,","1,3,,P","2,3,HH,P"]})

# Four external
#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T1','T1'],
#                   IDX_SV   :[1,2,3,4],
#                   LABEL_DESCR:["1,2,,PP","1,3,H,","2,3,,P","1,4,H,","2,4,,P"]})

#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:0.5,
#                   OPERATORS:['T2^+','H','T2','T1'],
#                   IDX_SV   :[1,2,3,4],
#                   LABEL_DESCR:["1,2,,P","1,3,HH,","2,3,,PP","1,4,,P","2,4,H,"]})

#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:0.5,
#                   OPERATORS:['T2^+','H','T2','T1'],
#                   IDX_SV   :[1,2,3,4],
#                   LABEL_DESCR:["1,2,,P","1,3,H,P","2,3,H,P","1,4,H,","2,4,,P"]})

#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:0.5,
#                   OPERATORS:['T2^+','H','T2','T1'],
#                   IDX_SV   :[1,2,3,4],
#                   LABEL_DESCR:["1,2,,P","1,3,HH,P","2,3,,P","2,4,H,P"]})

#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:False,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T1','T1','T1'],
#                   IDX_SV   :[1,2,3,4],
#                   LABEL_DESCR:["1,2,,P","1,3,H,","2,3,,P","1,4,H,","2,4,,P","1,5,,P","2,5,H,"]})



PRINT_FORMULA({LABEL:'CCD_FORM',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'Ext3',MODE:'SHORT'})

#FACTOR_OUT({LABEL_RES:'CCD_FORM',LABEL_IN:'CCD_FORM',INTERM:'Ext3'})
PRINT_FORMULA({LABEL:'CCD_FORM',MODE:'SHORT'})


# Residual
new_target('CCD_RES')
depend('CCD_FORM')

DERIVATIVE({LABEL_RES:'CCD_RES2',LABEL_IN:'CCD_FORM',
            OP_RES:'O2',OP_DERIV:'T2^+'})
DERIVATIVE({LABEL_RES:'CCD_RES1',LABEL_IN:'CCD_FORM',
            OP_RES:'O1',OP_DERIV:'T1^+'})

PRINT_FORMULA({LABEL:'CCD_RES1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'CCD_RES2',MODE:'SHORT'})
#TEX_FORMULA({LABEL:'CCD_RES1',OUTPUT:'texform1.tex'})
#TEX_FORMULA({LABEL:'CCD_RES2',OUTPUT:'texform2.tex'})

# Energy
new_target('CCD_EN')
depend('CCD_FORM')
form_e1 = stf.Formula("CCD_EN:ECCD=<H>")
form_e1.append("<[H,T2]>")
form_e1.append("<[H,T1]>")
#form_e1.append("<(1/2)*[[H,T1],T1]>")
form_e1.set_rule()
PRINT_FORMULA({LABEL:'CCD_EN',MODE:'SHORT'})


# ME lists ------------------------------------------------------
new_target('CCD_LISTS')
depend('CCD_OPS')
DEF_ME_LIST({LIST:'ME_LCCD',OPERATOR:'LCCD',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_ECCD',OPERATOR:'ECCD',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_T1',OPERATOR:'T1',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_T2',OPERATOR:'T2',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_O1',OPERATOR:'O1',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_O2',OPERATOR:'O2',IRREP:1,'2MS':0,AB_SYM:+1})


# Optimise formula -----------------------------------------------
new_target('CCD_OPT')
depend('CCD_RES','CCD_EN','CCD_LISTS','H0')
OPTIMIZE({LABEL_OPT:'CCD_OPT',LABELS_IN:['CCD_RES1','CCD_RES2','CCD_EN']})


TRANSLATE_ITF({LABEL: 'CCD_OPT',
               OUTPUT:'icmrcc_ccsd.itfaa',
               TITLE: 'icmrcc_ccsd.formulae',
               MULTI:False})
PRINT_FORMULA({LABEL:'CCD_OPT',MODE:'SHORT'})


# Set up preconditioner -----------------------------------------
new_target('DIAG')
depend('H0','CCD_OPS')
CLONE_OPERATOR({LABEL:'D2',TEMPLATE:'T2'})
DEF_ME_LIST({LIST:'DIAG2',OPERATOR:'D2',IRREP:1,'2MS':0,AB_SYM:+1})
PRECONDITIONER({LIST_PRC:'DIAG2',LIST_INP:'H0'})

CLONE_OPERATOR({LABEL:'D1',TEMPLATE:'T1'})
DEF_ME_LIST({LIST:'DIAG1',OPERATOR:'D1',IRREP:1,'2MS':0,AB_SYM:+1})
PRECONDITIONER({LIST_PRC:'DIAG1',LIST_INP:'H0'})


# Solve equations -----------------------------------------------
new_target('CCD_METHOD')
required()
depend('DIAG')
depend('CCD_OPT')

SOLVE_NLEQ({LIST_OPT:['ME_T1','ME_T2'],
            LIST_RESID:['ME_O1','ME_O2'],
            LIST_PRC:['DIAG1','DIAG2'],
            MODE:'DIA DIA',LIST_E:'ME_ECCD',FORM:'CCD_OPT'})

export_targets();

# 3 ext integrals
#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:True,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T1^+','H','T2'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,2,,P","1,3,H,","2,3,H,PP"]})

#EXPAND_OP_PRODUCT({LABEL:'Ext3',NEW:True,OP_RES:'LCCD',FAC:1.0,
#                   OPERATORS:['T2^+','H','T1'],
#                   IDX_SV   :[1,2,3],
#                   LABEL_DESCR:["1,2,H,PP","1,3,H,","2,3,,P"]})
