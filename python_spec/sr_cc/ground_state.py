#
# Set targets for general ground state CC (now using the python set up)
# 
#

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

import python_interface.gecco_modules.default_keywords as dk

# get input and orbital information
#inp = GeCCo_Input()
orb = Orb_Info()

msc = 1
if (orb.get('ims') != 0):
    msc = 0

# read main info about method
minexc = keywords.get('method.CC.minexc')
minexc = int(minexc) if minexc is not None else 1
maxexc = int(keywords.get('method.CC.maxexc'))
maxexc = int(maxexc) if maxexc is not None else 2
truncate = keywords.get('method.CC.truncate')
if truncate is None:
    truncate = "no"

# define the method label
mlabel = 'CC'
mchars = ['S','D','T','Q','P','H']
if minexc < 2 and maxexc <= 6:
    for ex in range(1,maxexc+1):
        if ex>= minexc:
            mlabel += mchars[ex-1]
elif minexc==1:
    mlabel = 'CC(X='+str(maxexc)+')'
else:
    mlabel = 'CC(strange)'
    

# define the operators
new_target('CC_OPS')
DEF_SCALAR({LABEL:'LCC'})  # Lagrange energy functional
DEF_SCALAR({LABEL:'ECC'})  # projective energy

# operator, diagonal, residual and LAMBDA
DEF_EXCITATION({LABEL:'T',MIN_RANK:minexc,MAX_RANK:maxexc})
CLONE_OPERATOR({LABEL:'DIAG',TEMPLATE:'T'})
CLONE_OPERATOR({LABEL:'OMG',TEMPLATE:'T'})
CLONE_OPERATOR({LABEL:'LAM',TEMPLATE:'T',ADJOINT:True})

DEF_OP_FROM_OCC({LABEL:'CC1DEN',JOIN:2,DESCR:',H;H,|,H;P,|,P;H,|,P;P,'})

# define the Lagrangian
new_target('CC_LAGRANGIAN')
depend('H')
depend('CC_OPS')

form_lag = stf.Formula("F_CC_LAG:LCC=<H>")

if truncate == 'no':
    PRINT({STRING:'Setting up the CC Lagrange functional'})
    form_lag.append("<[H,T]+(1/2)*[[H,T],T]>")
    form_lag.append("<LAM*H>")
    form_lag.append("<LAM*[H,T]>")
    form_lag.append("<LAM*(1/2)*[[H,T],T]>")
    form_lag.append("<LAM*(1/6)*[[[H,T],T],T]>")
    form_lag.append("<LAM*(1/24)*[[[[H,T],T],T],T]>")

    form_lag.set_rule()
else:
    quit_error('truncate not yet implemented, sorry')

# define CC equations:
new_target('CC_EQS')
depend('CC_LAGRANGIAN')

# energy
INVARIANT({LABEL_RES:'F_CC_EN',LABEL_IN:'F_CC_LAG',OP_RES:'ECC',OPERATORS:'LAM'})
# residual
DERIVATIVE({LABEL_RES:'F_CC_OMG',LABEL_IN:'F_CC_LAG',OP_RES:'OMG',OP_DERIV:'LAM'})

# set up the lists:
new_target('CC_LISTS')
depend('CC_OPS')
DEF_ME_LIST({LIST:'ME_LCC',OPERATOR:'LCC',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_ECC',OPERATOR:'ECC',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_T',OPERATOR:'T',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_DIAG',OPERATOR:'DIAG',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_OMG',OPERATOR:'OMG',IRREP:1,'2MS':0,AB_SYM:+1})

# do the computation
new_target('SOLVE_CC')
required()
depend('CC_LISTS','CC_EQS','H0')

# set up a diagonal for preconditioning
# the generic list for H is called H0 
PRECONDITIONER({LIST_PRC:'ME_DIAG',LIST_INP:'H0'})

PRINT({STRING:'Optimizing the CC equations for numerical evaluation.'})
OPTIMIZE({LABEL_OPT:'FOPT_CC',LABELS_IN:['F_CC_EN','F_CC_OMG']})

PRINT({STRING:'Now solving the CC equations:'})
PRINT_MEL_INFO({LIST:'ME_T'})

SOLVE_NLEQ({LIST_OPT:'ME_T',
            LIST_RESID:'ME_OMG',
            LIST_PRC:'ME_DIAG',
            MODE:'DIA',LIST_E:'ME_ECC',FORM:'FOPT_CC'})

PRINT_MEL({LIST:'ME_ECC',COMMENT:mlabel+' energy',FORMAT:'SCAL F24.14'})
PUSH_RESULT({LIST:'ME_ECC',COMMENT:mlabel, FORMAT:"SCAL F20.14"})

## Lambda equations:
new_target('CC_LAMBDA_EQS')
depend('CC_LAGRANGIAN')

CLONE_OPERATOR({LABEL:'LAM_A',TEMPLATE:'LAM'})
CLONE_OPERATOR({LABEL:'ETA',TEMPLATE:'LAM'})
DERIVATIVE({LABEL_RES:'F_CC_LAM',LABEL_IN:'F_CC_LAG',OP_RES:'LAM_A',OP_DERIV:'T'})
# split into RHS and Jacobian trafo:
LEQ_SPLIT({LABEL_RHS:'F_CC_ETA',LABEL_TRF:'F_CC_LAM_A',LABEL_RAW:'F_CC_LAM',
    OP_X:'LAM',OP_RHS:'ETA',OP_TRF:'LAM_A'})

# make also the 1-particle density
DERIVATIVE({LABEL_RES:'F_CC1DEN',LABEL_IN:'F_CC_LAG',OP_RES:'CC1DEN',OP_DERIV:'H'})

# additional lists
new_target('CC_LAMBDA_LISTS')
depend('CC_LAMBDA_EQS')

DEF_ME_LIST({LIST:'ME_LAM',OPERATOR:'LAM',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_ETA',OPERATOR:'ETA',IRREP:1,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_LAM_A',OPERATOR:'LAM_A',IRREP:1,'2MS':0,AB_SYM:+1})

DEF_ME_LIST({LIST:'ME_CC1DEN',OPERATOR:'CC1DEN',IRREP:1,'2MS':0,AB_SYM:+1})


# solve the equations
new_target('SOLVE_LAMBDA')
depend('SOLVE_CC','CC_LAMBDA_LISTS','CC_LAMBDA_EQS')

PRINT({STRING:'Optimizing the LAMBDA equatons for numerical evaluation.'})
OPTIMIZE({LABEL_OPT:'FOPT_CC_LAM',LABELS_IN:['F_CC_ETA','F_CC_LAM_A']})
OPTIMIZE({LABEL_OPT:'FOPT_CC1DEN',LABELS_IN:['F_CC1DEN']})

PRINT({STRING:'Now solving the CC LAMBDA equations:'})

SOLVE_LEQ({LIST_OPT:'ME_LAM',
           LIST_PRC:'ME_DIAG',
           OP_MVP:'LAM_A',
           OP_SVP:'LAM',
           OP_RHS:'ETA',
           FORM:'FOPT_CC_LAM',
           MODE:'DIA',
           N_ROOTS:1})

EVALUATE({FORM:'FOPT_CC1DEN'})

PRINT({STRING:''})
PRINT({STRING:'Evaluating unrelaxed ground state properties (electronic part):'})
EVAL_PROP({RANK:1,DENS:'ME_CC1DEN'})

export_targets()

