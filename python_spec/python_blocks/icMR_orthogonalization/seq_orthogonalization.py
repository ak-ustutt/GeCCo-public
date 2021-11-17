
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import * 

i_am='seq_orthogonalization'

spinadapt = keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0

minexc= keywords.get('method.MR.minexc')
minexc = int(minexc) if minexc is not None else 1

maxexc= keywords.get('method.MR.maxexc')
maxexc = int(maxexc) if maxexc is not None else 1
# perturbative correction requested?
word = keywords.get('method.MR.pertCorr')
if word is None:
    pertCorr =  False
else:
    if word == "F":
        pertCorr = False
    elif word == "T":
        pertCorr = True
    else:
        quit_error('pertCorr must be T or F, found: '+word)


if (minexc==1):
  _Tv_shape='V,H|P,V|P,H|PP,VV|PV,HV|VV,HH|PV,VV|VV,VH'
  _Ov_shape=',;V,H|,V;P,|,;P,H|,VV;PP,|,V;PV,H|,;VV,HH|,VV;PV,|,V;VV,H'
  _GAM_S_shape=',V;VV,VV;V,|,V;VV,V;,|,;V,VV;V,|,;V,V;,|,;VV,VV;,|,VV;V,V;VV,|,VV;V,;V,|,V;,V;VV,|,V;,;V,|,V;V,V;V,|,V;V,;,|,;,V;V,|,;,;,|,VV;,;VV,'
  _X_TRM_shape='VV,VV;V,V|VV,V;,V|V,VV;V,|V,V;,|VV,VV;,|V,V;VV,VV|V,;V,VV|,V;VV,V|,;V,V|V,V;V,V|V,;,V|,V;V,|,;,|,;VV,VV'
  useT1=True
else:
  _Tv_shape='V,H|P,V|P,H|PP,VV|PV,HV|VV,HH|PV,VV|VV,VH' # need P,H to generate ,;,;,
  _Ov_shape=',;V,H|,V;P,|,;P,H|,VV;PP,|,V;PV,H|,;VV,HH|,VV;PV,|,V;VV,H'
  _GAM_S_shape=',V;VV,VV;V,|,;V,V;,|,;VV,VV;,|,VV;V,V;VV,|,V;,;V,|,V;V,V;V,|,VV;,;VV,|,;,;,' # <-- this I meant above
  _X_TRM_shape='VV,VV;V,V|V,V;,|VV,VV;,|V,V;VV,VV|,;V,V|V,V;V,V|,;VV,VV|,;,'
  useT1=False

  # TODO: extend for triples etc.
if (maxexc>2 or pertCorr):
    _Tv_shape+='|PPP,VVV|VVV,HHH'
    _Ov_shape+='|,VVV;PPP,|,;VVV,HHH'
    _GAM_S_shape+='|,;VVV,VVV;,|,VVV;,;VVV,'
    _X_TRM_shape+='|VVV,VVV;,|,;VVV,VVV'

  
_s2 = orbitals.get('imult')

_ms = orbitals.get('ims')

if ((_ms == 0) and ((_s2-1 % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and ((_s2+1 % 4) == 0)):
    _msc = -1
else:
    _msc = 0
###################################################################
###################################################################
# ... set up transformation matrices for orthog. basis
###################################################################
###################################################################
new_target('MakeOrthBasis')
depend('DEF_T')
depend('DEF_O')
depend('DEF_LAM')

depend('MakeRefState')
depend('GAM0_CALC')


DEF_SCALAR({
        LABEL:'SSCAL'})

DEF_OP_FROM_OCC({
        LABEL:'1v_WE',
        DESCR:'V,V|VV,VV'})
SET_HERMITIAN({
        LABEL:'1v_WE',
        CA_SYMMETRY:+1})

DEF_ME_LIST({
        LIST:'1vLST',
        OPERATOR:'1v_WE',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1,
        DIAG_TYPE:1})
ADD_UNITY({
        LIST:'1vLST',
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

if spinadapt >= 2 : # and gno >=0 
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





CLONE_OPERATOR({
        LABEL:'T2_orth',
        TEMPLATE:'T2g'})
DEF_ME_LIST({
        LIST:'ME_T2_orth',
        OPERATOR:'T2_orth',
        IRREP:1,
        '2MS':0,
        AB_SYM:0})



# Formula for overlap matrix
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

debug_FORM('FORM_SMAT')

INSERT({
        LABEL_RES:'FORM_SMAT',
        LABEL_IN:'FORM_SMAT',
        OP_RES:'SSCAL',
        OP_INS:'1v_WE',
        OP_INCL:['Tv^+','Tv']})

# Now the double derivative, eliminating T^+ T
comment('Derivative 1')
DERIVATIVE({
        LABEL_RES:'FORM_TGAM0',
        LABEL_IN:'FORM_SMAT',
        OP_RES:'Ov',
        OP_DERIV:'Tv^+'})

comment('Derivative 2')

DERIVATIVE({
        LABEL_RES:'FORM_GAM_S',
        LABEL_IN:'FORM_TGAM0',
        OP_RES:'GAM_S',
        OP_DERIV:'Tv'})

debug_FORM('FORM_GAM_S')







# Evaluate the overlap densities
comment('Evaluate effective densities for overlap ...')
OPTIMIZE({
        LABEL_OPT:'FOPT_GAM_S',
        LABELS_IN:['FORM_GAM_S',"FORM_GAM0"]})

EVALUATE({
        FORM:'FOPT_GAM_S'})

debug_MEL('ME_GAM_S')
 
if False : # (gno >= 0 and l_iccc) 
    if (spinadapt >= 2): 
        SPIN_PROJECT({ 
                LIST:'ME_GAM_S',
                S2:0})
    EVALUATE({
            FORM:'',  # was  'FOPT_D_GNO'
            INIT:False
            })    


# Singular value decomposition
# does this overwrite ME_GAM_S ? yes it does! ME_GAM_S is now the projector
INVERT({
        LIST_INV:'ME_GAM_S',
        LIST:'ME_GAM_S_ISQ',
        MODE:'invsqrt'}) 

debug_MEL('ME_GAM_S_ISQ')





# Reordering of S^-0.5' to Transformation matrix X (X_Tr(ansformation)M(atrix))
DEF_OP_FROM_OCC({
        LABEL:'X_TRM',
        JOIN:2,
        DESCR:_X_TRM_shape})

CLONE_OPERATOR({
        LABEL:'X_TRM_DAG',
        TEMPLATE:'X_TRM'})
DEF_ME_LIST({
        LIST:'ME_X_TRM',
        OPERATOR:'X_TRM',
        IRREP:1,
        '2MS':0,
        S2:0,
        AB_SYM:0,
        CA_SYM:0})


DEF_ME_LIST({LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM_DAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:0})


REORDER_MEL({
        LIST_RES:'ME_X_TRM',
        LIST_IN:'ME_GAM_S_ISQ',
        FROMTO: 13})




#IMPORTANT: transposes on the input list. 
REORDER_MEL({
        LIST_RES:'ME_X_TRM_DAG',
        LIST_IN:'ME_GAM_S_ISQ',
        FROMTO: 13,
        ADJOINT:True})


debug_MEL('ME_X_TRM')



CLONE_OPERATOR({
        LABEL:'P_PROJ',
        TEMPLATE:'X_TRM'})

DEF_ME_LIST({
        LIST:'ME_P_PROJ',
        OPERATOR:'P_PROJ',
        IRREP:1,
        '2MS':0,
        S2:0,
        AB_SYM:0,
        CA_SYM:0})
REORDER_MEL({
        LIST_RES:'ME_P_PROJ',
        LIST_IN:'ME_GAM_S',
        FROMTO:13
})


# Transformation Formula for Excitation operators (OMEGA and T2_orth are 
# transformed by the same Formula with different definitions of X_TRM)  
# T2_orth:T_Tr(ans)f(ormed) 
#t=t'*X
EXPAND_OP_PRODUCT({
        LABEL:'FORM_T2_orth',
        NEW:True,
        OP_RES:'T2g',
        OPERATORS:['T2g','X_TRM','T2_orth','X_TRM','T2g'],
        IDX_SV:[1,2,3,2,1],
# No self contractions of the density matrix, no open lines from the wrong vertex
        AVOID:[2,4,1,4,2,5],
        })

#delete all terms where T2_orth has active external lines. external lines should be modified by X_trm

SELECT_LINE({
        LABEL_IN:'FORM_T2_orth',
        LABEL_RES:'FORM_T2_orth',
        OP_RES:'T2g',
        OP_INCL:'T2_orth',
        IGAST:3,
        MODE:'no_ext'})

debug_FORM('FORM_T2_orth')

OPTIMIZE({
        LABEL_OPT:'FOPT_T2_orth',
        LABELS_IN:'FORM_T2_orth'})


#Building the adjoint of FORM_T2_orth to transform T^+ # only needed for Building of A_trf

CLONE_OPERATOR({
        LABEL:'LAM_TRF',
        TEMPLATE:'LAM2g'})


EXPAND_OP_PRODUCT({
        LABEL:'FORM_LAM_TRF',
        NEW:True,
        OP_RES:'LAM2g',
        OPERATORS:['LAM2g','X_TRM^+','LAM_TRF','X_TRM^+','LAM2g'],
        IDX_SV:[1,2,3,2,1],
# No self contractions of the density matrix, no open lines from the wrong vertex
        AVOID:[2,4,
               1,4,
               2,5],
        })

SELECT_LINE({
        LABEL_IN:'FORM_LAM_TRF',
        LABEL_RES:'FORM_LAM_TRF',
        OP_RES:'LAM2g',
        OP_INCL:'LAM_TRF',
        IGAST:3,         #virtual
        MODE:'no_ext'})

debug_FORM('FORM_LAM_TRF')


#Formula to transform T1 only
#

if (useT1):
  CLONE_OPERATOR({
        LABEL:'T1_orth',
        TEMPLATE:'T1'})
  DEF_ME_LIST({
        LIST:'ME_T1_orth',
        OPERATOR:'T1_orth',
        IRREP:1,
        '2MS':0,
        AB_SYM:0})

  EXPAND_OP_PRODUCT({
        LABEL:'FORM_T1_orth',
        NEW:True,
        OP_RES:'T1',
        OPERATORS:['T1','X_TRM','T1_orth','X_TRM','T1'],
        IDX_SV:[1,2,3,2,1],
        AVOID:[2,4,
               1,4,
               2,5],
        })

  #delete all terms where T1_orth has active external lines. external lines should be modified by X_trm
  SELECT_LINE({
        LABEL_IN:'FORM_T1_orth',
        LABEL_RES:'FORM_T1_orth',
        OP_RES:'T1',
        OP_INCL:'T1_orth',
        IGAST:3,
        MODE:'no_ext'})

  OPTIMIZE({
        LABEL_OPT:'FOPT_T1_orth',
        LABELS_IN:'FORM_T1_orth'})
# end of block for (useT1)
  debug_FORM('FORM_T1_orth')


OPTIMIZE({
        LABEL_OPT:'FOPT_GES',
        LABELS_IN:['FORM_T2_orth']})


