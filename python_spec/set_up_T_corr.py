# Prepare everything for a further (T)-correction
# Based on the experimental file of Matthias for the 11/15/14 fix for (T)
#
# Probably just a provisional solution
#
# Yuri, Nov 2014
#

import sys,os
sys.path=sys.path+[os.getenv("GECCO_DIR")+"/python_interface"]


from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

msc = 1
if (orb.get('ims') != 0):
    msc = 0

solve_mrcc = 'SOLVE_MRCC'
if (inp.is_keyword_set('method.R12')):
    solve_mrcc += '_F12' 

densmix = inp.get('method.MR.densmix')

# create the Tfix list, to be used by the (T) correction
new_target('save_Tfix')
depend(solve_mrcc)
DEF_ME_LIST({LIST:'ME_Tfix',
             OPERATOR:'T',
             '2MS':0,
             IRREP:1,
             AB_SYM:msc})
SCALE_COPY({LIST_RES:'ME_Tfix',
            LIST_INP:'ME_T',
            FAC:1.0})


# define product of metric and cluster operator: S*t
new_target('FORM_ST')
depend('L','DEF_ME_T','DEF_ME_OMG','DEF_ME_C0','NORM','save_Tfix')
EXPAND_OP_PRODUCT({LABEL:'F_LST',
                   OP_RES:'NORM',
                   OPERATORS:['C0^+','L','T','C0'],
                   IDX_SV:   [     1,  2,  3,   4]})

# Make sure to use the proper (averaged) density
if (densmix>0):
    depend('F_DENS0')
    FACTOR_OUT({LABEL_RES:'F_LST',
                LABEL_IN:'F_LST',
                INTERM:'F_DENS0'})

DERIVATIVE({LABEL_RES:'F_ST',
            LABEL_IN:'F_LST',
            OP_RES:'OMG',
            OP_DERIV:'L'})
#PRINT_FORMULA({LABEL:'F_ST'})
OPTIMIZE({LABEL_OPT:'FOPT_ST',
          LABELS_IN:'F_ST'})


# truncate transformation formula such that only T2' is computed in the orthogonal basis
new_target('FOPT_T_extract')
depend ('F_T','DEF_ME_Dtr','DEF_ME_T','DEF_ME_Ttr')
SELECT_SPECIAL({LABEL_RES:'F_T_extract',
                LABEL_IN:'F_T',
                TYPE:'rank',
                MODE:'22',
                OPERATORS:['T','Ttr']})
#PRINT_FORMULA({LABEL:'F_T_extract'})
OPTIMIZE({LABEL_OPT:'FOPT_T_extract',
          LABELS_IN:'F_T_extract'})


# carry out the projection t_new = X P2' X^T S t_old = X P2' X^(-1) t_old
new_target('EVAL_extract')
depend ('FORM_ST','DIAG1SxxM00_T','EVAL_D','DEF_ME_Dtrdag','FOPT_T_extract','EVAL_Atr')
# evaluate metric vector product
EVALUATE({FORM:'FOPT_ST'})
# apply sign correction
SCALE_COPY({LIST_RES:'ME_T',
            LIST_INP:'ME_OMG',
            LIST_SHAPE:'ME_OMG'})
# transform to orthogonalized basis (and extract wanted rank in this step)
ASSIGN_ME2OP({LIST:'ME_Dtrdag',OPERATOR:'Dtr'})
ASSIGN_ME2OP({LIST:'ME_Ttr',OPERATOR:'T'})
ASSIGN_ME2OP({LIST:'ME_T',OPERATOR:'Ttr'})
EVALUATE({FORM:'FOPT_T_extract'})
# transform back to original basis
ASSIGN_ME2OP({LIST:'ME_Dtr',OPERATOR:'Dtr'})
ASSIGN_ME2OP({LIST:'ME_T',OPERATOR:'T'})
ASSIGN_ME2OP({LIST:'ME_Ttr',OPERATOR:'Ttr'})
EVALUATE({FORM:'FOPT_T_extract'})
PRINT_MEL({LIST:'ME_T',
           FORMAT:'NORM F20.12',
           COMMENT:'>>> Norm of extracted T operator :'})

# create the T2fix list, to be used by the (T) correction
new_target('save_T2fix',True)
depend('EVAL_extract')
DEF_ME_LIST({LIST:'ME_T2fix',
             OPERATOR:'T',
             '2MS':0,
             IRREP:1,
             AB_SYM:msc})
SCALE_COPY({LIST_RES:'ME_T2fix',
            LIST_INP:'ME_T',
            FAC:1.0})


export_targets();
