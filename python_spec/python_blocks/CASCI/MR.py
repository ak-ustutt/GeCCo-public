
from python_interface.gecco_interface import *


from python_interface.gecco_modules.NoticeUtil import *



spinadapt = keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0

ciroot=keywords.get('method.MR_P.ciroot')
ciroot = int(ciroot) if ciroot is not None else 1

maxroots = keywords.get('method.MR_P.maxroot')
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





OPTIMIZE({
        LABEL_OPT:'FOPT_GAM0',
        LABELS_IN:'FORM_GAM0'})
PRINT({
        STRING:'Setup densities up to 3rd order'})

new_target('GAM0_CALC')#for comparison to cut of the reference state calculation 
depend("FOPT_GAM0")
EVALUATE({
        FORM:'FOPT_GAM0'})



new_target("EVAL_E0")
depend("MakeRefState")
OPTIMIZE({
        LABEL_OPT:'FOPT_EREF',
        LABELS_IN:'FORM_EREF'})

EVALUATE({
        FORM:'FOPT_EREF'})


