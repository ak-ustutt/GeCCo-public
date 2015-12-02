
from gecco_interface import *


from gecco_modules.NoticeUtil import *



spinadapt=0
if keywords.is_keyword_set('calculate.routes.spinadapt'):
    spinadapt=int(keywords.get('calculate.routes.spinadapt'))

ciroot=1
if keywords.is_keyword_set('method.MR_P.ciroot')
    ciroot=int(keywords.get('method.MR_P.ciroot'))
print "ciroot: "+str(ciroot)


maxroots=ciroot
if keywords.is_keyword_set('method.MR_P.maxroot')
    maxroots=int(keywords.get('method.MR_P.maxroot'))
print "maxroot: "+str(maxroot)
###################################################################
# ---- solve for reference state
###################################################################
new_target('MakeRefState')
if spinadapt > 0:
    depend('SpinProjSpec')
depend('RefState-Operators')
depend('H0')


heading('Preparing equations ...')
# Formula for Eigenvalue-Problem 
#E0=<|C0^+*H*C0|>
#dE0/dC0^+
#H_C0=H*C0         <= E0*C0=H*C0
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
comment('Prepare diagonal ...')
PRECONDITIONER({
        LIST_PRC:'D0_LST',
        LIST_INP:'H0',
        MODE:'dia-H'})


comment('Solving equations ...')

SOLVE_map={
        LIST_OPT:'C0_LST',
        LIST_PRC:'D0_LST',
        OP_MVP:'H_C0',
        OP_SVP:'C0',
        FORM:'FOPT_H_C0',
        MODE:'DIA',
        N_ROOTS:maxroots
        TARG_ROOT:ciroot
}

if (spinadapt != 0): #and refproj = 0  
    SOLVE_map[MODE]='SPP' #'DIA' will be overwritten
    SOLVE_map[LIST_SPC]='C0_sp_LST'
    SOLVE_map[FORM_SPC]='FOPT_C0_sp'

elif False: # ref_proj !=0
    SOLVE_map[MODE]='PRJ' #'DIA' will be overwritten
    SOLVE_map[FORM_SPC]='FOPT_C0_prj' #not yet defined

SOLVE_EVP(SOLVE_map)






new_target('GAM0_CALC')#for comparison to cut of the reference state calculation 
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
EVALUATE({
        FORM:'FOPT_GAM0'})








