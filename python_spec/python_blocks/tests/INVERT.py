from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import * 

spinadapt = keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0
  
_s2 = orbitals.get('imult')

_ms = orbitals.get('ims')

if ((_ms == 0) and ((_s2-1 % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and ((_s2+1 % 4) == 0)):
    _msc = -1
else:
    _msc = 0

project= keywords.get('method.MR.project')
print("project:",project)
project = int(project) if project is not None else 2
    
new_target('TEST_INVERT',True)
depend('MakeRefState')
depend('GAM0_CALC')


#Note: this testing routine does some basic testing of INVERT. there are many options that are not hit here

if project in(2,3,4) :
    _Tv_shape='V,H|P,V|P,H|PP,VV|PV,HV|VV,HH|PV,VV|VV,VH'
    _Ov_shape=',;V,H|,V;P,|,;P,H|,VV;PP,|,V;PV,H|,;VV,HH|,VV;PV,|,V;VV,H'
    _GAM_S_shape=',V;VV,VV;V,|,V;VV,V;,|,;V,VV;V,|,;V,V;,|,;VV,VV;,|,VV;V,V;VV,|,VV;V,;V,|,V;,V;VV,|,V;,;V,|,V;V,V;V,|,V;V,;,|,;,V;V,|,;,;,|,VV;,;VV,'
    _X_TRM_shape='VV,VV;V,V|VV,V;,V|V,VV;V,|V,V;,|VV,VV;,|V,V;VV,VV|V,;V,VV|,V;VV,V|,;V,V|V,V;V,V|V,;,V|,V;V,|,;,|,;VV,VV'
    useT1=True
elif project == 1 :
  _Tv_shape='V,H|P,V|P,H|PP,VV|PV,HV|VV,HH|PV,VV|VV,VH' # need P,H to generate ,;,;,
  _Ov_shape=',;V,H|,V;P,|,;P,H|,VV;PP,|,V;PV,H|,;VV,HH|,VV;PV,|,V;VV,H'
  _GAM_S_shape=',V;VV,VV;V,|,;V,V;,|,;VV,VV;,|,VV;V,V;VV,|,V;,;V,|,V;V,V;V,|,VV;,;VV,|,;,;,' # <-- this I meant above
  _X_TRM_shape='VV,VV;V,V|V,V;,|VV,VV;,|V,V;VV,VV|,;V,V|V,V;V,V|,;VV,VV|,;,'
  useT1=False
else:
    raise Exception("unknown projector for Test")

#       just setting up some helper entities

DEF_SCALAR({
        LABEL:'SSCAL'})

DEF_OP_FROM_OCC({
        LABEL:'1v',
        DESCR:'V,V|VV,VV'})
SET_HERMITIAN({
        LABEL:'1v',
        CA_SYMMETRY:+1})



DEF_ME_LIST({
        LIST:'ME_1v',
        OPERATOR:'1v',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1,
        DIAG_TYPE:1})
ADD_UNITY({
        LIST:'ME_1v',
        FAC:1.0,
        INIT:True})




DEF_OP_FROM_OCC({
        LABEL:'Tv',
        DESCR:_Tv_shape})

DEF_OP_FROM_OCC({
        LABEL:'Ov',
        JOIN:2,
        DESCR:_Ov_shape})

DEF_OP_FROM_OCC({
        LABEL:'GAM_S',
        JOIN:3,
        DESCR:_GAM_S_shape})

if spinadapt >= 2 : # and gno >= 0
    S2_val=0
else:
    S2_val=-1  # the default

DEF_ME_LIST({
        LIST:'ME_GAM_S',
        OPERATOR:'GAM_S',
        IRREP:1,
        '2MS':0,
        AB_SYM:_msc,
        S2:S2_val})


EXPAND_OP_PRODUCT({
        LABEL:'FORM_SMAT',
        NEW:True,
        OP_RES:'SSCAL',
        OPERATORS:['C0^+','Tv^+','Tv','C0'],
        IDX_SV:[1,2,3,4]})

FACTOR_OUT({
        LABEL_RES:'FORM_SMAT',
        LABEL_IN:'FORM_SMAT',
        INTERM:'FORM_GAM0'})




INSERT({
        LABEL_RES:'FORM_SMAT',
        LABEL_IN:'FORM_SMAT',
        OP_RES:'SSCAL',
        OP_INS:'1v',
        OP_INCL:['Tv^+','Tv']})

DERIVATIVE({
        LABEL_RES:'FORM_TGAM0',
        LABEL_IN:'FORM_SMAT',
        OP_RES:'Ov',
        OP_DERIV:'Tv^+'})

DERIVATIVE({
        LABEL_RES:'FORM_GAM_S',
        LABEL_IN:'FORM_TGAM0',
        OP_RES:'GAM_S',
        OP_DERIV:'Tv'})

OPTIMIZE({
        LABEL_OPT:'FOPT_GAM_S',
        LABELS_IN:'FORM_GAM_S'})

EVALUATE({
        FORM:'FOPT_GAM_S'})




debug_FORM('FORM_GAM_S', only_this=True)


CLONE_OPERATOR({
        LABEL:'ISQ_GAM',
        TEMPLATE:'GAM_S'})

DEF_ME_LIST({
        LIST:'ME_GAM_S_ISQ',
        OPERATOR:'ISQ_GAM',
        IRREP:1,
        '2MS':0,
        'S2':0,
        AB_SYM:0,
        CA_SYM:0})


INVERT({
        LIST_INV:'ME_GAM_S',
        LIST:'ME_GAM_S_ISQ',
        MODE:'invsqrt'}) 



debug_MEL('ME_GAM_S', only_this=True)

debug_MEL('ME_GAM_S_ISQ', only_this=True)
