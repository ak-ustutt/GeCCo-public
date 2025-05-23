from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

i_am='diag_precon'

minexc= keywords.get('method.MR.minexc')
minexc = int(minexc) if minexc is not None else 1

maxexc= keywords.get('method.MR.maxexc')
maxexc = int(maxexc) if maxexc is not None else 1

useT1=(minexc==1)

# perturbative triples are handled differently, keep this for non-iterative cases
if maxexc>2:
    word = keywords.get('method.MR.triples')
    if word is None:
      triples=3
    elif word == "B" or word == "3":
      triples=3
    elif word == "E" or word == "4":
      triples=4
    elif word == "F" or word == "5":
      triples=5
    else:
      quit_error('triples must be one of B,E,F,3,4,5; found: '+word)
    quit_error('full triples: diag_precon requires adaptation')
      
#------------------------------------------------------------------
#Setting up A_TRF
#-------------------------------------------------------------------

new_target('MAKE_A_TRF')
depend('MakeOrthBasis')

depend('DEF_T')

depend('DEF_LAM')

depend('DEF_O')

depend('MakeRefState')
depend('H0')

comment('Setting up A_TRF')

_A_TRF_shape=',;V,V;,|'\
              ',V;VV,VV;V,|'\
              ',;VV,VV;,|'\
              ',V;,;V,|'\
              ',VV;V,V;VV,|'\
              ',V;V,V;V,|'\
              ',VV;,;VV,'

DEF_SCALAR({
        LABEL:'A_TRF_SCAL'})
DEF_OP_FROM_OCC({
        LABEL:'A_TRF',
        JOIN:3,
        DESCR:_A_TRF_shape})
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
_hamil_='H'
# there is a option to try more here: FOCK_EFF alone is not good, needs shifts (I guess)
#_hamil_='FOCK_EFF'
#depend('MAKE_FOCK_EFF')

EXPAND_OP_PRODUCT({
        LABEL:'FORM_A_TRF',
        OP_RES:'A_TRF_SCAL',
        OPERATORS:['C0^+','LAM2g',_hamil_,'T2g','C0'],
        IDX_SV:[1,2,3,4,5],
        NEW:True})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_A_TRF',OP_RES:'A_TRF_SCAL',
        OPERATORS:['C0^+','LAM2g','T2g',_hamil_,'C0'],
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
        INTERM:'FORM_LAM_TRF'})
EXPAND({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF',
        INTERM:'FORM_T2_orth'})
SELECT_SPECIAL({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF',
        OPERATORS:['T2_orth','LAM_TRF'],
        TYPE:'SAME'})

#Double derivative to T' and T'^+ to get A_TRF=<|t^+[H,t]|>
DERIVATIVE({
        LABEL_IN:'FORM_A_TRF',
        LABEL_RES:'FORM_A_TRF_INTERM',
        OP_RES:'Ov',
        OP_DERIV:'LAM_TRF'})
DERIVATIVE({
        LABEL_IN:'FORM_A_TRF_INTERM',
        LABEL_RES:'FORM_A_TRF_FINAL',
        OP_RES:'A_TRF',
        OP_DERIV:'T2_orth'})



debug_FORM('FORM_A_TRF_FINAL')


OPTIMIZE({
        LABEL_OPT:'FOPT_A_TRF_FINAL',
        LABELS_IN:'FORM_A_TRF_FINAL'})
#PRINT_FORMULA({LABEL:'FOPT_A_TRF_FINAL', MODE:'SHORT'})
#PRINT_FORMULA({LABEL:'FORM_A_TRF_FINAL', MODE:'SHORT'})


EVALUATE({
        FORM:'FOPT_A_TRF_FINAL'})


debug_MEL('A_TRF_LST')










#-------------------------------------------------------------------
#building the actual preconditioner
#-------------------------------------------------------------------
new_target('BUILD_PRECON')
depend('MAKE_A_TRF')
depend('EVAL_F_EFF_INACT')



CLONE_OPERATOR({
        LABEL:'PRECON',
        TEMPLATE:'T2g'})
DEF_ME_LIST({
        LIST:'ME_PRECON2g',
        OPERATOR:'PRECON',
        IRREP:1,
        '2MS':0,
        AB_SYM:0}) # AB not well defined for orthog. basis!
#        AB_SYM:+1})

# Uncomment to get a print out of the active preconditioner
#CLONE_OPERATOR({
#        LABEL:'TEST_PRECON',
#        TEMPLATE:'T2g'})
#DEF_ME_LIST({
#        LIST:'TEST',
#        OPERATOR:'TEST_PRECON',
#        IRREP:1,
#        '2MS':0,
#        AB_SYM:0})
#DEF_ME_LIST({
#        LIST:'TEST_DIFF',
#        OPERATOR:'TEST_PRECON',
#        IRREP:1,
#        '2MS':0,
#        AB_SYM:0})

#extract preconditioner
PRECONDITIONER({
        LIST_PRC:'ME_PRECON2g',
        LIST_INP:'FOCK_EFF_INACT_LST'})
#PRECONDITIONER({
#        LIST_PRC:'TEST',
#        LIST_INP:'FOCK_EFF_INACT_LST'})
#
#PRINT_MEL({
#        LIST:'ME_PRECON2g',
#        COMMENT:'hello: inactive'})

EXTRACT_DIAG({
        LIST_RES:'ME_PRECON2g',
        LIST_IN:'A_TRF_LST',
        MODE:'extend'})
SCALE_COPY({
        LIST_RES:'ME_PRECON2g',
        LIST_INP:'ME_PRECON2g',
        FAC:0.2,
        MODE:'prc_thresh'})

#PRINT_MEL({
#        LIST:'ME_PRECON2g',
#        COMMENT:'hello: inactive + active'})

#PRINT_MEL({
#        LIST:'TEST',
#        COMMENT:'hello: test1'})
#ADD({LIST_SUM:'TEST_DIFF',
#     LISTS:['TEST', 'ME_PRECON2g'],
#     FAC:[-1,1]})
#
#PRINT_MEL({
#        LIST:'TEST_DIFF',
#        COMMENT:'hello: test2'})

PRINT_MEL_INFO({
        LIST:'ME_PRECON2g'})

#PRINT_MEL({
#        LIST:'A_TRF_LST',
#        COMMENT:'hello: active'})

debug_MEL('ME_PRECON2g')

#----------------------------------------------------------
#Preconditioner for T1
#----------------------------------------------------------

if (useT1):
  CLONE_OPERATOR({
        LABEL:'PRECON1',
        TEMPLATE:'T1'})
  DEF_ME_LIST({
        LIST:'ME_PRECON1',
        OPERATOR:'PRECON1',
        IRREP:1,
        '2MS':0,
        AB_SYM:0})

  #extract preconditioner
  PRECONDITIONER({
        LIST_PRC:'ME_PRECON1',
        LIST_INP:'FOCK_EFF_INACT_LST'})

  EXTRACT_DIAG({
        LIST_RES:'ME_PRECON1',
        LIST_IN:'A_TRF_LST',
        MODE:'extend'})
  SCALE_COPY({
        LIST_RES:'ME_PRECON1',
        LIST_INP:'ME_PRECON1',
        FAC:0.2,
        MODE:'prc_thresh'})

  debug_MEL('ME_PRECON1')


