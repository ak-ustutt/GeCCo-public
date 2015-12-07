from gecco_interface import *
from python_blocks.Arnes_helfer import *

###################################################################
#Setting up A_TRF
###################################################################


new_target('MAKE_A_TRF')
depend('MakeOrthBasis')
depend('T-Operators')
depend('MakeRefState')
depend('H0')

comment('Setting up A_TRF')

DEF_SCALAR({
        LABEL:'A_TRF_SCAL'})
DEF_OP_FROM_OCC({
        LABEL:'A_TRF',
        JOIN:3,
        DESCR:',;V,V;,|,V;VV,VV;V,|,;VV,VV;,|,V;,;V,|,VV;V,V;VV,|,V;V,V;V,|,VV;,;VV,'})
DEF_ME_LIST({
        LIST:'A_TRF_LST',
        OPERATOR:'A_TRF',
        IRREP:1,
        '2MS':0})







#The actual Formula for A_TRF
#A_TRF=C0^+*t'^+*X^+*[H,t'*X]*C0
#procedure:
#SCAL=<C0^+*L*[H,T]*C0>
#Insert transformation formulas for L and T
#derive with respect to L' and T'

EXPAND_OP_PRODUCT({
        LABEL:'FORM_A_TRF',
        OP_RES:'A_TRF_SCAL',
        OPERATORS:['C0^+','L_NONTRF','H','T2_ca','C0'],
        IDX_SV:[1,2,3,4,5],
        NEW:True})
EXPAND_OP_PRODUCT({LABEL:'FORM_A_TRF',OP_RES:'A_TRF_SCAL',
        OPERATORS:['C0^+','L_NONTRF','T2_ca','H','C0'],
        IDX_SV:[1,2,3,4,5],
        FIX_VTX:True,
        FAC:-1,
        NEW:False})
#eliminate all terms that are ruled out by the commutation(show up with +1 and -1 as prefactor)
SUM_TERMS({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF'})

FACTOR_OUT({
        LABEL_RES:'FORM_A_TRF',
        LABEL_IN:'FORM_A_TRF',
        INTERM:'FORM_GAM0'})

EXPAND({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF',
        INTERM:'FORM_L_TRF'})
EXPAND({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF',
        INTERM:'FORM_T2_orth'})
SELECT_SPECIAL({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF', 
        OPERATORS:['T2_orth','L_TRF'], 
        TYPE:'SAME'})

#Double derivative to T' and T'^+ to get A_TRF=<|t^+[H,t]|>
DERIVATIVE({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF_INTERM',
        OP_RES:'Ov',
        OP_DERIV:'L_TRF'})
DERIVATIVE({
        LABEL_IN:'FORM_A_TRF_INTERM',
        LABEL_RES:'FORM_A_TRF_FINAL',
        OP_RES:'A_TRF',
        OP_DERIV:'T2_orth'})



debug_FORM('FORM_A_TRF_FINAL')


OPTIMIZE({
        LABEL_OPT:'FOPT_A_TRF_FINAL',
        LABELS_IN:'FORM_A_TRF_FINAL'})
EVALUATE({
        FORM:'FOPT_A_TRF_FINAL'})

debug_MEL('A_TRF_LST')
