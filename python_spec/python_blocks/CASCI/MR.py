
from python_interface.gecco_interface import *


from python_interface.gecco_modules.NoticeUtil import *

oldref = False
word = keywords.get('method.MR.oldref')
if word is not None:
    if word=='T':
        oldref = True

spinadapt = keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0

# parse MR_P for backward compatibility
ciroot=keywords.get('method.MR_P.ciroot')
ciroot = int(ciroot) if ciroot is not None else 1
if ciroot==1:
    ciroot=keywords.get('method.MR.ciroot')
    ciroot = int(ciroot) if ciroot is not None else 1

maxroots = keywords.get('method.MR_P.maxroot')
maxroots = int(maxroots) if maxroots is not None else ciroot
if maxroots==ciroot:
    maxroots = keywords.get('method.MR.maxroot')
    maxroots = int(maxroots) if maxroots is not None else ciroot

###################################################################
# ---- solve for reference state
###################################################################




new_target('MakeRefState')
if spinadapt > 0:
    depend('SpinProjSpec')
depend('RefState-Operators')
depend('MAKE_D0')


heading('Preparing equations ...')
# Formula for Eigenvalue-Problem 
#E0=<|C0^+*H*C0|>
#dE0/dC0^+ :       H_C0=H*C0   <= E0*C0=H*C0
EXPAND_OP_PRODUCT({
        LABEL:'FORM_EREF',
        NEW:True,
        OP_RES:'E0',
        OPERATORS:['C0^+','H','C0'],
        IDX_SV:[1,2,3]})

DERIVATIVE({
        LABEL_RES:'FORM_H_C0',
        LABEL_IN:'FORM_EREF',
        OP_RES:'H_C0',
        OP_DERIV:'C0^+'})



OPTIMIZE({
        LABEL_OPT:'FOPT_H_C0',
        LABELS_IN:'FORM_H_C0'})

if not oldref:
    comment('Solving equations ...')

    SOLVE_map={
            LIST_OPT:'ME_C0',
            LIST_PRC:'ME_D0',
            OP_MVP:'H_C0',
            OP_SVP:'C0',
            FORM:'FOPT_H_C0',
            MODE:'DIA',
            N_ROOTS:maxroots,
            TARG_ROOT:ciroot,
    }

    if (spinadapt != 0): #and refproj = 0  
        SOLVE_map[MODE]='SPP' #'DIA' will be overwritten
        SOLVE_map[LIST_SPC]='ME_C0_sp'
        SOLVE_map[FORM_SPC]='FOPT_C0_sp'
    
    elif False: # ref_proj !=0
        SOLVE_map[MODE]='PRJ' #'DIA' will be overwritten
        SOLVE_map[FORM_SPC]='FOPT_C0_prj' #not yet defined
    
    SOLVE_EVP(SOLVE_map)

    # get a copy of the reference function
    SCALE_COPY({LIST_RES:"ME_C00",
                LIST_INP:"ME_C0",
                FAC:1.0})
else:
    PRINT({STRING:'Restart active, assuming files for C0 and C00 in place ...'})

# evaluate reference energy from C00 (needed for restart)
REPLACE({LABEL_IN:'FORM_EREF',LABEL_RES:'FORM_EREF',OP_LIST:['C0','C00','C0^+','C00^+']})

OPTIMIZE({
        LABEL_OPT:'FOPT_EREF',
        LABELS_IN:'FORM_EREF'})

EVALUATE({
        FORM:'FOPT_EREF'})

PRINT_MEL({LIST:'ME_E0',COMMENT:'Reference energy',FORMAT:'SCAL F24.14'})
PUSH_RESULT({LIST:'ME_E0',COMMENT:"Reference", FORMAT:"SCAL F24.14"})


new_target("FOPT_GAM0")
depend('MakeRefState')
depend('RefState-Operators')
# Formula for Densities
#mark:C0 change
EXPAND_OP_PRODUCT({
        LABEL:'FORM_GAM0',
        NEW:True,
        OP_RES:'GAM0',
        OPERATORS:['GAM0','C0^+','GAM0','GAM0','C0','GAM0'],
        IDX_SV:[1,2,1,1,3,1]})

debug_FORM('FORM_GAM0')

EXPAND_OP_PRODUCT({
        LABEL:'FORM_GAM00',
        NEW:True,
        OP_RES:'GAM00',
        OPERATORS:['GAM00','C00^+','GAM00','GAM00','C00','GAM00'],
        IDX_SV:[1,2,1,1,3,1]})




OPTIMIZE({
        LABEL_OPT:'FOPT_GAM0',
        LABELS_IN:'FORM_GAM0'})
OPTIMIZE({
        LABEL_OPT:'FOPT_GAM00',
        LABELS_IN:'FORM_GAM00'})
PRINT({
        STRING:'Setup density matrices'})

new_target('GAM0_CALC')
depend("FOPT_GAM0")
EVALUATE({
        FORM:'FOPT_GAM0'})

debug_MEL('GAM0_LST')


new_target('GAM00_CALC')
depend("FOPT_GAM0")
EVALUATE({
        FORM:'FOPT_GAM00'})

debug_MEL('GAM00_LST')


# YAA tmp: Small formula for testing the density
# Just for testing.
# Why does the density needs AB_sym=msc for a
# correct behaviour?
#
# EXPAND_OP_PRODUCT({
#         LABEL:'FORM_GAM0_cpy',
#         NEW:True,
#         OP_RES:'GAM0_cpy',
#         OPERATORS:['GAM0_cpy','C0^+','GAM0_cpy','GAM0_cpy','C0','GAM0_cpy'],
#         IDX_SV:[1,2,1,1,3,1]})
# debug_FORM('FORM_GAM0_cpy')
# OPTIMIZE({
#         LABEL_OPT:'FOPT_GAM0_cpy',
#         LABELS_IN:'FORM_GAM0_cpy'})
# EVALUATE({
#         FORM:'FOPT_GAM0_cpy'})
# debug_MEL('GAM0_LST')
# debug_MEL('GAM0_LST_cpy')


new_target("EVAL_E0")
depend("MakeRefState")
# moved evaluation of E0 into target MakeRefState
