"""Solve/evaluate the icMRCC (T) equations

History:

Andreas Nov 2021: Creation of initial stub based on some earlier experiments

"""
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *
#test MRCCSD+(T)
#from python_spec.python_blocks.solve.icMRCCSDsolve.py import MRCC_LAG_LST
verbosity=00

run_checks=False

i_am="icMRCCSDpTsolve.py"

_inp = GeCCo_Input()
# get some info:
_orb = Orb_Info()
_ninact_o = _orb.get('nactt_hpv',1)
_nact_o = _orb.get('nactt_hpv',3)
#_nact_o = _orb.get('nactorb')
_nact_e = _orb.get('nactel')

_spinadapt = _inp.get('calculate.routes.spinadapt')
if _spinadapt == None:
  _spinadapt = 3
  if _s2 == 1:
    _spinadapt = 0

_ms = 0
if _spinadapt:
   _s2 = 0
else:
   _s2 = -1
_msc = +1

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

pTH0 = "Dyall"
word = keywords.get('method.MR.H0')
if word is None:
    pTH0 = "Dyall"
elif word.upper() == "DYALL":
    pTH0 = "Dyall"
elif word.upper() == "FOCK":
    pTH0 = "Fock"
else:
    quit_error('H0 must be either "Dyall" or "Fock"')

shift = 0.0
word = keywords.get('method.MR.shift')
if word is None:
    shfit = 0.0
else:
    try:
        shift = float(word)
    except:
        quit_error('Undigestible input for shift: '+word)


new_target('SOLVE_MRCC_PT')
depend('SOLVE_MRCC')
depend('MakeOrthBasisPT') ## compute an extended set of lists for orth --> X_TRM_PT, X_TRM_PT_DAG
heading('Computing the (T) correction')

#DBG:
if run_checks:
    PRINT({STRING:'X_TRM:'})
    ANALYZE_MEL({LISTS:'ME_X_TRM',LISTS_CV:'ME_X_TRM'})
    PRINT({STRING:'X_TRM_DAG:'})
    ANALYZE_MEL({LISTS:'ME_X_TRM_DAG',LISTS_CV:'ME_X_TRM_DAG'})
    PRINT({STRING:'X_TRM_PT:'})
    ANALYZE_MEL({LISTS:'ME_X_TRM_PT',LISTS_CV:'ME_X_TRM_PT'})
    PRINT({STRING:'X_TRM_PT_DAG:'})
    ANALYZE_MEL({LISTS:'ME_X_TRM_PT_DAG',LISTS_CV:'ME_X_TRM_PT_DAG'})
    PRINT({STRING:'ME_P_PROJ_PT:'})
    ANALYZE_MEL({LISTS:'ME_P_PROJ_PT',LISTS_CV:'ME_P_PROJ_PT'})
    debug_MEL('ME_P_PROJ_PT')#,only_this=True)
    debug_MEL('ME_X_TRM_PT')#,only_this=True)


# for "FOCK" we have to compute <F> for current density (note: the density defining F_EFF is still the CASCI density!)
if (pTH0 == "Fock"):
   OPTIMIZE({
        LABELS_IN:['FORM_F_EFF_EXP'],
        LABEL_OPT:'FOPT_F_EFF_EXP'})
   EVALUATE({
        FORM:'FOPT_F_EFF_EXP'})


# set the various occupation classes
# might be done more elegantly, here we just enumerate everything:
Tblocks={}
if _ninact_o > 1:
   Tblocks['PPP-HHH']={}
   Tblocks['PPP-HHH']['occ']='PPP,HHH'
   Tblocks['PPP-HHH']['res']=',;PPP,HHH'
   Tblocks['PPP-HHH']['nact']=0
   Tblocks['PPP-HHH']['shiftfac']=0.

if _ninact_o > 0:
   Tblocks['PPP-HHV']={}
   Tblocks['PPP-HHV']['occ']='PPP,HHV'
   Tblocks['PPP-HHV']['res']=',V;PPP,HH'
   Tblocks['PPP-HHV']['nact']=1
   Tblocks['PPP-HHV']['shiftfac']=1.

if _ninact_o > 1:
   Tblocks['PPV-HHH']={}
   Tblocks['PPV-HHH']['occ']='PPV,HHH'
   Tblocks['PPV-HHH']['res']=',;PPV,HHH'
   Tblocks['PPV-HHH']['nact']=1
   Tblocks['PPV-HHH']['shiftfac']=1.

if _ninact_o > 0:
   Tblocks['PPP-HVV']={}
   Tblocks['PPP-HVV']['occ']='PPP,HVV'
   Tblocks['PPP-HVV']['res']=',VV;PPP,H'
   Tblocks['PPP-HVV']['nact']=2
   Tblocks['PPP-HVV']['shiftfac']=1.

if _ninact_o > 1:
   Tblocks['PVV-HHH']={}
   Tblocks['PVV-HHH']['occ']='PVV,HHH'
   Tblocks['PVV-HHH']['res']=',;PVV,HHH'
   Tblocks['PVV-HHH']['nact']=2
   Tblocks['PVV-HHH']['shiftfac']=1.

if _ninact_o > 0:
   Tblocks['PPV-HHV']={}
   Tblocks['PPV-HHV']['occ']='PPV,HHV|PP,HH'
   Tblocks['PPV-HHV']['res']=',V;PPV,HH|,;PP,HH'
   Tblocks['PPV-HHV']['nact']=2
   Tblocks['PPV-HHV']['shiftfac']=2.

if _nact_e > 2:
   Tblocks['PPP-VVV']={}
   Tblocks['PPP-VVV']['occ']='PPP,VVV'
   Tblocks['PPP-VVV']['res']=',VVV;PPP,'
   Tblocks['PPP-VVV']['nact']=3
   Tblocks['PPP-VVV']['shiftfac']=1.
   Tblocks['PPP-VVV']['A_extra']=',VVV;,;VVV,'
   Tblocks['PPP-VVV']['A_offdia']=''

if ((2*_nact_o) - _nact_e) > 2 and _ninact_o > 1:
   Tblocks['VVV-HHH']={}
   Tblocks['VVV-HHH']['occ']='VVV,HHH'
   Tblocks['VVV-HHH']['res']=',;VVV,HHH'
   Tblocks['VVV-HHH']['nact']=3
   Tblocks['VVV-HHH']['shiftfac']=1.
   Tblocks['VVV-HHH']['A_extra']=',;VVV,VVV;,'
   Tblocks['VVV-HHH']['A_offdia']=''

if _ninact_o > 0: # and _nact_e > 1:
   Tblocks['PPV-HVV']={}
   Tblocks['PPV-HVV']['occ']='PPV,HVV|PP,HV'
   Tblocks['PPV-HVV']['res']=',VV;PPV,H|,V;PP,H'
   Tblocks['PPV-HVV']['nact']=3
   Tblocks['PPV-HVV']['shiftfac']=2.

if _ninact_o > 0:  
   Tblocks['PVV-HHV']={}
   Tblocks['PVV-HHV']['occ']='PVV,HHV|PV,HH'
   Tblocks['PVV-HHV']['res']=',V;PVV,HH|,;PV,HH'
   Tblocks['PVV-HHV']['nact']=3
   Tblocks['PVV-HHV']['shiftfac']=2.

# Should be redundant for CAS(4,4) according to Hanauer 2012

if (triples>3):
    if _ninact_o > 0: #and _nact_e > 3:
       Tblocks['PVV-HVV']={}
       Tblocks['PVV-HVV']['occ']='PVV,HVV|PV,HV|P,H'
       Tblocks['PVV-HVV']['res']=',VV;PVV,H|,V;PV,H|,;P,H'
       Tblocks['PVV-HVV']['nact']=4
       Tblocks['PVV-HVV']['shiftfac']=2.
       Tblocks['PVV-HVV']['A_extra']=',VV;VV,VV;VV,|,V;V,V;V,'
       Tblocks['PVV-HVV']['A_offdia']='|,VV;VV,V;V,|,V;V,VV;VV,'

    Tblocks['PPV-VVV']={}
    Tblocks['PPV-VVV']['occ']='PPV,VVV|PP,VV'
    Tblocks['PPV-VVV']['res']=',VVV;PPV,|,VV;PP,'
    Tblocks['PPV-VVV']['nact']=4
    Tblocks['PPV-VVV']['shiftfac']=2.
    Tblocks['PPV-VVV']['A_extra']=',VVV;V,V;VVV,|,VV;,;VV,'
    Tblocks['PPV-VVV']['A_offdia']='|,VVV;V,;VV,|,VV;,V;VVV,'

    if _ninact_o > 0:
       Tblocks['VVV-HHV']={}
       Tblocks['VVV-HHV']['occ']='VVV,HHV|VV,HH'
       Tblocks['VVV-HHV']['res']=',V;VVV,HH|,;VV,HH'
       Tblocks['VVV-HHV']['nact']=4
       Tblocks['VVV-HHV']['shiftfac']=2.
       Tblocks['VVV-HHV']['A_extra']=',V;VVV,VVV;V,|,;VV,VV;,'
       Tblocks['VVV-HHV']['A_offdia']='|,V;VVV,VV;,|,;VV,VVV;V,'

if (triples>4):
    if _ninact_o > 0:
       Tblocks['VVV-HVV']={}
       Tblocks['VVV-HVV']['occ']='VVV,HVV|VV,HV|V,H'
       Tblocks['VVV-HVV']['res']=',VV;VVV,H|,V;VV,H|,;V,H'
       Tblocks['VVV-HVV']['nact']=5
       Tblocks['VVV-HVV']['shiftfac']=2.
       Tblocks['VVV-HVV']['A_extra']=',VV;VVV,VVV;VV,|,V;VV,VV;V,|,;V,V;,'
       Tblocks['VVV-HVV']['A_offdia']='|,VV;VVV,VV;V,|VVV;V,V;,|,V;VV,VVV;VV,|,V;VV,V;,|,;V,VVV;VV,|,;V,VV;V,'

    Tblocks['PVV-VVV']={}
    Tblocks['PVV-VVV']['occ']='PVV,VVV|PV,VV|P,V'
    Tblocks['PVV-VVV']['res']=',VVV;PVV,|,VV;PV,|,V;P,'
    Tblocks['PVV-VVV']['nact']=5
    Tblocks['PVV-VVV']['shiftfac']=2.
    Tblocks['PVV-VVV']['A_extra']=',VVV;VV,VV;VVV,|,VV;V,V;VV,|,V;,;V,'
    Tblocks['PVV-VVV']['A_offdia']='|,VVV;VV,V;VV,|,VVV;VV,;V,|,VV;V,VV;VVV,|,VV;V,;V,|,V;,VV;VVV,|,V;,V;VV,'

# for later: we need a dummy operator for virtual space only:
DEF_OP_FROM_OCC({LABEL:'HVV_PT',DESCR:'V,V|VV,VV'})

    
#  DEF_FORMULA_TESTSTUFF
#  TOTAL SUM
DEF_SCALAR({LABEL:'EPTtotal'})
DEF_ME_LIST({LIST:'ME-EPTtotal',OPERATOR:'EPTtotal',IRREP:1,'2MS':0,'ABSYM':0})

DEF_SCALAR({LABEL:'Etotal'})
DEF_ME_LIST({LIST:'ME-Etotal',OPERATOR:'Etotal',IRREP:1,'2MS':0,'ABSYM':0})

DEF_SCALAR({LABEL:'Etotal_B'})
DEF_ME_LIST({LIST:'ME-Etotal_B',OPERATOR:'Etotal_B',IRREP:1,'2MS':0,'ABSYM':0})
DEF_SCALAR({LABEL:'Etotal_E'})
DEF_ME_LIST({LIST:'ME-Etotal_E',OPERATOR:'Etotal_E',IRREP:1,'2MS':0,'ABSYM':0})
DEF_SCALAR({LABEL:'Etotal_F'})
DEF_ME_LIST({LIST:'ME-Etotal_F',OPERATOR:'Etotal_F',IRREP:1,'2MS':0,'ABSYM':0})

# EPT4 EP5 scalar and ME-list test
EPT4str ='EPT4total=' # fills the EPT4 string for the energy summation later on 
DEF_SCALAR({LABEL:'EPT4total'})
DEF_ME_LIST({LIST:'ME-EPT4total',OPERATOR:'EPT4total',IRREP:1,'2MS':0,'ABSYM':0})

EPT5str ='EPT5total=' # fills the EPT5 string for the summation later 
DEF_SCALAR({LABEL:'EPT5total'}) 
DEF_ME_LIST({LIST:'ME-EPT5total',OPERATOR:'EPT5total',IRREP:1,'2MS':0,'ABSYM':0})
# Cleaner output.
#ETRIP_Bstr # shoiws energy up to {3} (3 active indices)
ETRIP4_Bstr ='ETRIP4_B='
DEF_SCALAR({LABEL:'ETRIP4_B'})
DEF_ME_LIST({LIST:'ME-ETRIP4_B',OPERATOR:'ETRIP4_B',IRREP:1,'2MS':0,'ABSYM':0})
ETRIP5_Bstr ='ETRIP5_B='
DEF_SCALAR({LABEL:'ETRIP5_B'})
DEF_ME_LIST({LIST:'ME-ETRIP5_B',OPERATOR:'ETRIP5_B',IRREP:1,'2MS':0,'ABSYM':0})
#ETRIP_Estr # shows energy up to {4} (4 active indices)
ETRIP4_Estr ='ETRIP4_E='
DEF_SCALAR({LABEL:'ETRIP4_E'})
DEF_ME_LIST({LIST:'ME-ETRIP4_E',OPERATOR:'ETRIP4_E',IRREP:1,'2MS':0,'ABSYM':0})
ETRIP5_Estr ='ETRIP5_E='
DEF_SCALAR({LABEL:'ETRIP5_E'})
DEF_ME_LIST({LIST:'ME-ETRIP5_E',OPERATOR:'ETRIP5_E',IRREP:1,'2MS':0,'ABSYM':0})
#ETRIP_Fstr # shows energy up to {5} once it works
ETRIP4_Fstr ='ETRIP4_F='
DEF_SCALAR({LABEL:'ETRIP4_F'})
DEF_ME_LIST({LIST:'ME-ETRIP4_F',OPERATOR:'ETRIP4_F',IRREP:1,'2MS':0,'ABSYM':0})
ETRIP5_Fstr ='ETRIP5_F='
DEF_SCALAR({LABEL:'ETRIP5_F'})
DEF_ME_LIST({LIST:'ME-ETRIP5_F',OPERATOR:'ETRIP5_F',IRREP:1,'2MS':0,'ABSYM':0})

# for DBG
if run_checks:
    DEF_SCALAR({LABEL:'OVL-TEST'})
    DEF_ME_LIST({LIST:'ME-OVL-TEST',OPERATOR:'OVL-TEST',IRREP:1,'2MS':0,'ABSYM':0})


# generate some operators to address the individual T3 blocks:
for _Tb in Tblocks:

   DEF_OP_FROM_OCC({LABEL:'T3-'+_Tb,DESCR:Tblocks[_Tb]['occ']})
   DEF_ME_LIST({LIST:'ME-T3-'+_Tb,OPERATOR:'T3-'+_Tb,IRREP:1,'2MS':_ms,'S2':_s2,'ABSYM':_msc})
   # maybe: use the T3g list instead? will require a bit of change to the formulae below, but I think it is doable ...

   # transformed T3 ... will probably need this
   CLONE_OPERATOR({LABEL:'T3tr-'+_Tb,TEMPLATE:'T3-'+_Tb})
   DEF_ME_LIST({LIST:'ME-T3tr-'+_Tb,OPERATOR:'T3tr-'+_Tb,IRREP:1,'2MS':_ms,'ABSYM':_msc})

   # create own diagonal
   CLONE_OPERATOR({LABEL:'DIA-'+_Tb,TEMPLATE:'T3-'+_Tb})
   DEF_ME_LIST({LIST:'ME-DIA-'+_Tb,OPERATOR:'DIA-'+_Tb,IRREP:1,'2MS':_ms,'ABSYM':_msc})

   DEF_OP_FROM_OCC({LABEL:'O3-'+_Tb,JOIN:2,DESCR:Tblocks[_Tb]['res']})

   if (verbosity>=100):
       PRINT({STRING:'-->'+Tblocks[_Tb]['res']})
       PRINT_OP_OCC({LABEL:'O3-'+_Tb})
   CLONE_OPERATOR({LABEL:'PTrhs-'+_Tb,TEMPLATE:'O3-'+_Tb})
   DEF_ME_LIST({LIST:'ME-PTrhs-'+_Tb,OPERATOR:'PTrhs-'+_Tb,IRREP:1,'2MS':_ms,'S2':_s2,'ABSYM':_msc})
   CLONE_OPERATOR({LABEL:'PTtrf-'+_Tb,TEMPLATE:'O3-'+_Tb})
   DEF_ME_LIST({LIST:'ME-PTtrf-'+_Tb,OPERATOR:'PTtrf-'+_Tb,IRREP:1,'2MS':_ms,'S2':_s2,'ABSYM':_msc})

   DEF_SCALAR({LABEL:'EPT-'+_Tb})
   DEF_ME_LIST({LIST:'ME-EPT-'+_Tb,OPERATOR:'EPT-'+_Tb,IRREP:1,'2MS':0,'ABSYM':0})
   DEF_SCALAR({LABEL:'EPT4-'+_Tb})
   DEF_ME_LIST({LIST:'ME-EPT4-'+_Tb,OPERATOR:'EPT4-'+_Tb,IRREP:1,'2MS':0,'ABSYM':0})
   DEF_SCALAR({LABEL:'EPT5-'+_Tb})
   DEF_ME_LIST({LIST:'ME-EPT5-'+_Tb,OPERATOR:'EPT5-'+_Tb,IRREP:1,'2MS':0,'ABSYM':0})
   DEF_SCALAR({LABEL:'EPTLG-'+_Tb})
   DEF_ME_LIST({LIST:'ME-EPTLG-'+_Tb,OPERATOR:'EPTLG-'+_Tb,IRREP:1,'2MS':0,'ABSYM':0})
   

# insert density in all equations:
FACTOR_OUT({LABEL_RES:'FORM_MRCC_PT_LAG',LABEL_IN:'FORM_MRCC_PT_LAG',INTERM:'FORM_GAM0'})
FACTOR_OUT({LABEL_RES:'FORM_MRCC_PT_E4',LABEL_IN:'FORM_MRCC_PT_E4',INTERM:'FORM_GAM0'})
FACTOR_OUT({LABEL_RES:'FORM_MRCC_PT_E5',LABEL_IN:'FORM_MRCC_PT_E5',INTERM:'FORM_GAM0'})


PRINT({STRING:'Computing additional denominators'})
# expressions for each block
for _Tb in Tblocks:
   # we replace L3 and T3 by the relevant subblock
   REPLACE({LABEL_RES:'F_PT_LG_'+_Tb,LABEL_IN:'FORM_MRCC_PT_LAG',
            OP_LIST:['T3g','T3-'+_Tb,'LAM3g','T3-'+_Tb+'^+']})
   # and remove the rest ... 
   INVARIANT({LABEL_RES:'F_PT_LG_'+_Tb,LABEL_IN:'F_PT_LG_'+_Tb,
              OP_RES:'EPTLG-'+_Tb,OPERATORS:['T3g','LAM3g']})
   # turn this into an residual equation
   DERIVATIVE({LABEL_IN:'F_PT_LG_'+_Tb,LABEL_RES:'F_PTeq_'+_Tb,OP_RES:'O3-'+_Tb,OP_DERIV:'T3-'+_Tb+'^+'})
   # finally: split into RHS and TRF (matrix-vector trafo) part
   LEQ_SPLIT({LABEL_RAW:'F_PTeq_'+_Tb,
              LABEL_RHS:'F_PTrhs_'+_Tb,LABEL_TRF:'F_PTtrf_'+_Tb,
              OP_RHS:'PTrhs-'+_Tb,OP_TRF:'PTtrf-'+_Tb,OP_X:'T3-'+_Tb})
   # 4th order
   REPLACE({LABEL_RES:'F_PT_E4_'+_Tb,LABEL_IN:'FORM_MRCC_PT_E4',
            OP_LIST:['T3g','T3-'+_Tb,'T3g^+','T3-'+_Tb+'^+']})
   INVARIANT({LABEL_RES:'F_PT_E4_'+_Tb,LABEL_IN:'F_PT_E4_'+_Tb,
              OP_RES:'EPT4-'+_Tb,OPERATORS:'T3g'})
   # 5th order
   REPLACE({LABEL_RES:'F_PT_E5_'+_Tb,LABEL_IN:'FORM_MRCC_PT_E5',
            OP_LIST:['T3g','T3-'+_Tb,'T3g^+','T3-'+_Tb+'^+']})
   INVARIANT({LABEL_RES:'F_PT_E5_'+_Tb,LABEL_IN:'F_PT_E5_'+_Tb,
              OP_RES:'EPT5-'+_Tb,OPERATORS:'T3g'})
   # Lagrange correction term ... not req.'d if we solve the equations
#   REPLACE({LABEL_RES:'F_PT_LG_'+_Tb,LABEL_IN:'FORM_MRCC_PT_LAG',
#            OP_LIST:['T3g','T3-'+_Tb,'T3g^+','T3-'+_Tb+'^+']})
#   INVARIANT({LABEL_RES:'F_PT_LG_'+_Tb,LABEL_IN:'F_PT_LG_'+_Tb,
#              OP_RES:'EPTLG-'+_Tb,OPERATORS:['T3g','T3g^+']})
   # we also need some transformations for the T3g block
   # ...
   debug_FORM('F_PTrhs_'+_Tb)#,only_this=True)
   debug_FORM('F_PTtrf_'+_Tb)#,only_this=True)
   debug_FORM('F_PT_E4_'+_Tb)#,only_this=True)
   debug_FORM('F_PT_E5_'+_Tb)#,only_this=True)
   debug_FORM('F_PT_LG_'+_Tb)#,only_this=True)

   # ... and the optimized formulas
   OPTIMIZE({LABEL_OPT:'FOPT_PTeq-'+_Tb,
             LABELS_IN:['F_PTrhs_'+_Tb,'F_PTtrf_'+_Tb]})
   OPTIMIZE({LABEL_OPT:'FOPT_PTE-'+_Tb,
             LABELS_IN:['F_PT_E4_'+_Tb,'F_PT_E5_'+_Tb,'F_PT_LG_'+_Tb]})
   
   # Sums all terms and adds a + in the end
   EPT4str = EPT4str + 'EPT4-' +_Tb + '+' 
   EPT5str = EPT5str + 'EPT5-' +_Tb + '+'
   if Tblocks[_Tb]["nact"]<4: #generate separate strings for better output for {3} etc..
       ETRIP4_Bstr =ETRIP4_Bstr + 'EPT4-' +_Tb + '+'
       ETRIP5_Bstr =ETRIP5_Bstr + 'EPT5-' +_Tb + '+'
   if Tblocks[_Tb]["nact"]<5:
       ETRIP4_Estr =ETRIP4_Estr + 'EPT4-' +_Tb + '+'
       ETRIP5_Estr =ETRIP5_Estr + 'EPT5-' +_Tb + '+'
   if Tblocks[_Tb]["nact"]<6:
       ETRIP4_Fstr =ETRIP4_Fstr + 'EPT4-' +_Tb + '+'
       ETRIP5_Fstr =ETRIP5_Fstr + 'EPT5-' +_Tb + '+'
   if Tblocks[_Tb]["nact"]>0:
      # and a transformation (that keeps T3 orth to T2 and T1)
      EXPAND_OP_PRODUCT({LABEL:'F_T3tr-'+_Tb,NEW:True,
                      OP_RES:'T3-'+_Tb,
                      OPERATORS:['T3-'+_Tb,'X_TRM_PT','T3tr-'+_Tb,'X_TRM_PT','T3-'+_Tb],
                      IDX_SV:   [1,2,3,2,1],
                      AVOID:[2,4,1,4,2,5]})
      SELECT_LINE({
        LABEL_IN:'F_T3tr-'+_Tb,
        LABEL_RES:'F_T3tr-'+_Tb,
        OP_RES:'T3-'+_Tb,
        OP_INCL:'T3tr-'+_Tb,
        IGAST:3,
        MODE:'no_ext'})
      debug_FORM('F_T3tr-'+_Tb)#,only_this=True)
      # we can remove some parts of the transformation:
      # for the transformation of the residual into the orthog. basis (where X_TRM is replaced by X_TRM_DAG)
      # we only need the contributions of the low-rank residuals to the largest rank one
      # for the transformation of the amplitude update, only the contributions from the largest rank will be
      # kept and distributed to the lower ranks to get the orthogonalization tails
      SELECT_SPECIAL({LABEL_RES:'F_T3tr-'+_Tb,LABEL_IN:'F_T3tr-'+_Tb,
             TYPE:'rank',MODE:'33',OPERATORS:['T3-'+_Tb,'T3tr-'+_Tb]})
      debug_FORM('F_T3tr-'+_Tb)#,only_this=True)
      OPTIMIZE({LABEL_OPT:'FOPT_T3tr-'+_Tb,
             LABELS_IN:'F_T3tr-'+_Tb})

      # DBG -- test orthogonality
      if Tblocks[_Tb]["nact"]>3 and run_checks:
          EXPAND_OP_PRODUCT({LABEL:'FORM_ORTH-'+_Tb,OP_RES:'OVL-TEST',
                OPERATORS:['C0^+','T1^+','T3-'+_Tb,'C0'],
                IDX_SV:[1,2,3,4],
                NEW:True})
          EXPAND_OP_PRODUCT({LABEL:'FORM_ORTH-'+_Tb,OP_RES:'OVL-TEST',
                OPERATORS:['C0^+','T2g^+','T3-'+_Tb,'C0'],
                IDX_SV:[1,2,3,4],
                NEW:False})
          debug_FORM('FORM_ORTH-'+_Tb,only_this=True)
          OPTIMIZE({LABEL_OPT:'FOPT_ORTH-'+_Tb,LABELS_IN:'FORM_ORTH-'+_Tb})

          # create formula to project T3 into the contravariant basis:
          CLONE_OPERATOR({LABEL:'ST3-'+_Tb,TEMPLATE:'O3-'+_Tb})
          DEF_ME_LIST({LIST:'ME-ST3-'+_Tb,OPERATOR:'ST3-'+_Tb,IRREP:1,'2MS':_ms,'S2':_s2,'ABSYM':_msc})
          EXPAND_OP_PRODUCT({LABEL:'F_ST3-0-'+_Tb,OP_RES:'OVL-TEST',  # OVL-TEST is only a dummy intermediate here
                OPERATORS:['C0^+','T3-'+_Tb+'^+','T3-'+_Tb,'C0'],
                IDX_SV:[1,2,3,4],
                NEW:True})
          FACTOR_OUT({LABEL_RES:'F_ST3-0-'+_Tb,LABEL_IN:'F_ST3-0-'+_Tb,INTERM:'FORM_GAM0'})
          DERIVATIVE({LABEL_IN:'F_ST3-0-'+_Tb,LABEL_RES:'F_ST3-'+_Tb,OP_RES:'ST3-'+_Tb,OP_DERIV:'T3-'+_Tb+'^+'})
          OPTIMIZE({LABEL_OPT:'FOPT-ST3-'+_Tb,LABELS_IN:'F_ST3-'+_Tb})


      # for some T blocks, we need to compute extra blocks of the A matrix
      if "A_extra" in Tblocks[_Tb].keys():
          
          _A_T_shape = Tblocks[_Tb]["A_extra"]
          _A_TX_shape = _A_T_shape + Tblocks[_Tb]["A_offdia"]
          DEF_SCALAR({LABEL:'A_TRF_SCAL-'+_Tb})
          DEF_OP_FROM_OCC({LABEL:'A_TRF-'+_Tb,JOIN:3,DESCR:_A_T_shape})
          DEF_OP_FROM_OCC({LABEL:'A_INT-'+_Tb,JOIN:3,DESCR:_A_TX_shape})
          #CLONE_OPERATOR({LABEL:'A_INT-'+_Tb,TEMPLATE:'A_TRF-'+_Tb})
          DEF_ME_LIST({LIST:'ME_A_TRF-'+_Tb,OPERATOR:'A_TRF-'+_Tb,IRREP:1,'2MS':0})
          DEF_ME_LIST({LIST:'ME_A_INT-'+_Tb,OPERATOR:'A_INT-'+_Tb,IRREP:1,'2MS':0})

          # define only the terms involving the active part of H ... all other terms would be eliminated below anyway
          EXPAND_OP_PRODUCT({LABEL:'FORM_A_TRF-'+_Tb,OP_RES:'A_TRF_SCAL-'+_Tb,
            OPERATORS:['C0^+','T3-'+_Tb+'^+','HVV_PT','T3-'+_Tb,'C0'],
            IDX_SV:[1,2,3,4,5],
            NEW:True})
          EXPAND_OP_PRODUCT({LABEL:'FORM_A_TRF-'+_Tb,OP_RES:'A_TRF_SCAL-'+_Tb,
            OPERATORS:['C0^+','T3-'+_Tb+'^+','T3-'+_Tb,'HVV_PT','C0'],
            IDX_SV:[1,2,3,4,5],
            FIX_VTX:True,FAC:-1,NEW:False})
          REPLACE({LABEL_RES:'FORM_A_TRF-'+_Tb,LABEL_IN:'FORM_A_TRF-'+_Tb,OP_LIST:['HVV_PT','H']})

          SUM_TERMS({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_TRF-'+_Tb})

          FACTOR_OUT({LABEL_RES:'FORM_A_TRF-'+_Tb,LABEL_IN:'FORM_A_TRF-'+_Tb,INTERM:'FORM_GAM0'})
          debug_FORM('FORM_A_TRF-'+_Tb)#,only_this=True)

          # insert the unit operator in V space: (is already defined when metric is set up for MRCCSD)
          # we need this in order to define a proper intermediate after taking the double derivative
          INSERT({LABEL_RES:'FORM_A_TRF-'+_Tb,LABEL_IN:'FORM_A_TRF-'+_Tb,
              OP_RES:'A_TRF_SCAL-'+_Tb,OP_INS:'1v_WE',OP_INCL:['T3-'+_Tb+'^+','T3-'+_Tb]})

          # generate intermediate without trafo matrix:
          DERIVATIVE({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_INT0-'+_Tb,OP_RES:'O3-'+_Tb,OP_DERIV:'T3-'+_Tb+'^+'})
          DERIVATIVE({LABEL_IN:'FORM_A_INT0-'+_Tb,LABEL_RES:'FORM_A_INT-'+_Tb,OP_RES:'A_INT-'+_Tb,OP_DERIV:'T3-'+_Tb})
          REORDER_FORMULA({LABEL_IN:'FORM_A_INT-'+_Tb,LABEL_RES:'FORM_A_INT-'+_Tb})
          debug_FORM('FORM_A_INT-'+_Tb)#,only_this=True)

          FACTOR_OUT({LABEL_RES:'FORM_A_TRF-'+_Tb,LABEL_IN:'FORM_A_TRF-'+_Tb,INTERM:'FORM_A_INT-'+_Tb})
          debug_FORM('FORM_A_TRF-'+_Tb)#,only_this=True)

          EXPAND({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_TRF-'+_Tb,INTERM:'F_T3tr-'+_Tb+'^+'})
          EXPAND({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_TRF-'+_Tb,INTERM:'F_T3tr-'+_Tb})
          debug_FORM('FORM_A_TRF-'+_Tb)#,only_this=True)
          # Probably unnecessary:
          SELECT_SPECIAL({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_TRF-'+_Tb,OPERATORS:['T3tr-'+_Tb,'T3tr-'+_Tb+'^+'],TYPE:'SAME'})
          debug_FORM('FORM_A_TRF-'+_Tb)#,only_this=True)
          
          DERIVATIVE({LABEL_IN:'FORM_A_TRF-'+_Tb,LABEL_RES:'FORM_A_TRF_INT-'+_Tb,OP_RES:'O3-'+_Tb,OP_DERIV:'T3tr-'+_Tb+'^+'})
          DERIVATIVE({LABEL_IN:'FORM_A_TRF_INT-'+_Tb,LABEL_RES:'FORM_A_TRF_FIN-'+_Tb,OP_RES:'A_TRF-'+_Tb,OP_DERIV:'T3tr-'+_Tb})
          debug_FORM('FORM_A_TRF_FIN-'+_Tb)#,only_this=True)

          OPTIMIZE({LABEL_OPT:'FOPT_A_TRF-'+_Tb,LABELS_IN:['FORM_A_INT-'+_Tb,'FORM_A_TRF_FIN-'+_Tb]})
          #OPTIMIZE({LABEL_OPT:'FOPT_A_TRF-'+_Tb,LABELS_IN:['FORM_A_TRF_FIN-'+_Tb]})

          EVALUATE({FORM:'FOPT_A_TRF-'+_Tb})
          debug_MEL('ME_A_TRF-'+_Tb)#,only_this=True)
          


# TESTING DEF_FORMULA

# removes the last + 

EPT4str=EPT4str[:-1]
EPT5str=EPT5str[:-1]
#print(EPT4str)
### better output for {3} etc..
ETRIP4_Bstr =ETRIP4_Bstr[:-1]
ETRIP5_Bstr =ETRIP5_Bstr[:-1]

ETRIP4_Estr =ETRIP4_Estr[:-1] 
ETRIP5_Estr =ETRIP5_Estr[:-1] 

ETRIP4_Fstr =ETRIP4_Fstr[:-1] 
ETRIP5_Fstr =ETRIP5_Fstr[:-1] 
# optimizing summation formula
DEF_FORMULA({LABEL:'FORM_EPT4tot',FORMULA:EPT4str})
OPTIMIZE({LABEL_OPT:'FOPT_EPT4tot',LABELS_IN:'FORM_EPT4tot' })

DEF_FORMULA({LABEL:'FORM_EPT5tot',FORMULA:EPT5str})
OPTIMIZE({LABEL_OPT:'FOPT_EPT5tot',LABELS_IN:'FORM_EPT5tot' })

DEF_FORMULA({LABEL:'FORM_EPTtot',FORMULA:'EPTtotal=EPT4total+EPT5total'})
OPTIMIZE({LABEL_OPT:'FOPT_EPTtot',LABELS_IN:'FORM_EPTtot' })

#DEF_SCALAR({LABEL:'Etotal'})
#DEF_ME_LIST({LIST:'ME-Etotal',OPERATOR:'Etotal',IRREP:1,'2MS':0,'ABSYM':0})
### better output for {3} -> could be more elegant though..
DEF_FORMULA({LABEL:'FORM_ETRIP4_B',FORMULA:ETRIP4_Bstr})
OPTIMIZE({LABEL_OPT:'FOPT_ETRIP4_B',LABELS_IN:'FORM_ETRIP4_B' })
DEF_FORMULA({LABEL:'FORM_ETRIP5_B',FORMULA:ETRIP5_Bstr})
OPTIMIZE({LABEL_OPT:'FOPT_ETRIP5_B',LABELS_IN:'FORM_ETRIP5_B' })
DEF_FORMULA({LABEL:'FORM_Etot_B',FORMULA:'Etotal_B=ETRIP4_B+ETRIP5_B'})
OPTIMIZE({LABEL_OPT:'FOPT_Etot_B',LABELS_IN:'FORM_Etot_B' })
if triples > 3:

    DEF_FORMULA({LABEL:'FORM_ETRIP4_E',FORMULA:ETRIP4_Estr})
    OPTIMIZE({LABEL_OPT:'FOPT_ETRIP4_E',LABELS_IN:'FORM_ETRIP4_E' })
    DEF_FORMULA({LABEL:'FORM_ETRIP5_E',FORMULA:ETRIP5_Estr})
    OPTIMIZE({LABEL_OPT:'FOPT_ETRIP5_E',LABELS_IN:'FORM_ETRIP5_E' })
    DEF_FORMULA({LABEL:'FORM_Etot_E',FORMULA:'Etotal_E=ETRIP4_E+ETRIP5_E'})
    OPTIMIZE({LABEL_OPT:'FOPT_Etot_E',LABELS_IN:'FORM_Etot_E' })

if triples > 4:

    DEF_FORMULA({LABEL:'FORM_ETRIP4_F',FORMULA:ETRIP4_Fstr})
    OPTIMIZE({LABEL_OPT:'FOPT_ETRIP4_F',LABELS_IN:'FORM_ETRIP4_F' })
    DEF_FORMULA({LABEL:'FORM_ETRIP5_F',FORMULA:ETRIP5_Fstr})
    OPTIMIZE({LABEL_OPT:'FOPT_ETRIP5_F',LABELS_IN:'FORM_ETRIP5_F' })
    DEF_FORMULA({LABEL:'FORM_Etot_F',FORMULA:'Etotal_F=ETRIP4_F+ETRIP5_F'})
    OPTIMIZE({LABEL_OPT:'FOPT_Etot_F',LABELS_IN:'FORM_Etot_F' })
# Trying to sum MRCCSD and (T)

DEF_FORMULA({LABEL:'FORM_Etot',FORMULA:'Etotal=EPTtotal+MRCC_LAG'})
OPTIMIZE({LABEL_OPT:'FOPT_Etot',LABELS_IN:'FORM_Etot' })

# Printing Test terms
debug_FORM('FORM_EPT4tot')
debug_FORM('FORM_EPT5tot')
debug_FORM('FORM_EPTtot')

# loop over blocks of T3
for _Tb in Tblocks:
   # set up preconditioner
   PRECONDITIONER({LIST_PRC:'ME-DIA-'+_Tb,LIST_INP:'FOCK_EFF_INACT_LST',MODE:'dia-F'})

   if Tblocks[_Tb]['nact']>0:
       if 'A_extra' not in Tblocks[_Tb].keys():
           # standard list:
           EXTRACT_DIAG({LIST_RES:'ME-DIA-'+_Tb,LIST_IN:'A_TRF_LST',MODE:'extend'})
       else:
           # special list:
           EXTRACT_DIAG({LIST_RES:'ME-DIA-'+_Tb,LIST_IN:'ME_A_TRF-'+_Tb,MODE:'extend'})
           
   # solve linear equations
   PRINT({STRING:'Now solving (T) equations for block '+Tblocks[_Tb]['occ']})
   if Tblocks[_Tb]['nact']>0:
      # apply shift if requested:
      if shift>0.0 and pTH0=="Fock":
          shift_block = shift*Tblocks[_Tb]['shiftfac']
          # have to subtract as ( F - <F> ) is evaluated 
          PRINT({STRING:'Applied a shift of {:4.2f} to this block'.format(shift_block)})
          MODIFY_BLOCKS({LIST:'FOCK_EFF_EXP_LST',FAC:-shift_block,DESCR:',',MODE:'SHIFT'})

      SOLVE_LEQ({LIST_OPT:'ME-T3-'+_Tb,MODE:'TRF',
              OP_MVP:'PTtrf-'+_Tb,OP_RHS:'PTrhs-'+_Tb,
              OP_SVP:'T3-'+_Tb,N_ROOTS:1,
              LIST_PRC:'ME-DIA-'+_Tb,
              FORM:'FOPT_PTeq-'+_Tb,
              LIST_SPC:['ME-T3-'+_Tb,'ME-T3tr-'+_Tb,'ME_X_TRM_PT','ME_X_TRM_PT_DAG'],
              FORM_SPC:'FOPT_T3tr-'+_Tb})

      # undo shift:
      if shift>0.0 and pTH0=="Fock":
          MODIFY_BLOCKS({LIST:'FOCK_EFF_EXP_LST',FAC:shift_block,DESCR:',',MODE:'SHIFT'})

      # DBG
      if Tblocks[_Tb]["nact"]>3 and run_checks:
          ANALYZE_MEL({LISTS:'ME-T3-'+_Tb,LISTS_CV:'ME-T3-'+_Tb})
          EVALUATE({FORM:'FOPT-ST3-'+_Tb})
          PRINT({STRING:'Biorthogonal norm:'})
          ANALYZE_MEL({LISTS:'ME-T3-'+_Tb,LISTS_CV:'ME-ST3-'+_Tb})
          PRINT({STRING:'Vector in contrav. basis:'})
          ANALYZE_MEL({LISTS:'ME-ST3-'+_Tb,LISTS_CV:'ME-ST3-'+_Tb})
          ASSIGN_ME2OP({LIST:'ME_X_TRM_PT_DAG',OPERATOR:'X_TRM_PT'}) # reassign trafo matrix
          ASSIGN_ME2OP({LIST:'ME-T3tr-'+_Tb,OPERATOR:'T3-'+_Tb}) # reassign result matrix
          ASSIGN_ME2OP({LIST:'ME-ST3-'+_Tb,OPERATOR:'T3tr-'+_Tb}) # reassign input matrix
          EVALUATE({FORM:'FOPT_T3tr-'+_Tb})
          PRINT({STRING:'Vector in orth. basis:'})
          ANALYZE_MEL({LISTS:'ME-T3tr-'+_Tb,LISTS_CV:'ME-T3tr-'+_Tb})
          # and back trafo
          ASSIGN_ME2OP({LIST:'ME_X_TRM_PT',OPERATOR:'X_TRM_PT'}) # reassign trafo matrix
          ASSIGN_ME2OP({LIST:'ME-T3-'+_Tb,OPERATOR:'T3-'+_Tb}) # reassign result matrix
          ASSIGN_ME2OP({LIST:'ME-T3tr-'+_Tb,OPERATOR:'T3tr-'+_Tb}) # reassign input matrix
          EVALUATE({FORM:'FOPT_T3tr-'+_Tb})
          PRINT({STRING:'Regenerated T3 vector:'})
          ANALYZE_MEL({LISTS:'ME-T3-'+_Tb,LISTS_CV:'ME-T3-'+_Tb})


          # test orthogonality
          EVALUATE({FORM:'FOPT_ORTH-'+_Tb})
          PRINT_MEL({LIST:'ME-OVL-TEST',
              COMMENT:'<(T1+T2)T3> = ',
              FORMAT:'SCAL F24.14'})


   else:
      # simpler case here (consider this non-iteratively)
      SOLVE_LEQ({LIST_OPT:'ME-T3-'+_Tb,MODE:'DIA',
              OP_MVP:'PTtrf-'+_Tb,OP_RHS:'PTrhs-'+_Tb,
              OP_SVP:'T3-'+_Tb,N_ROOTS:1,
              LIST_PRC:'ME-DIA-'+_Tb,
              FORM:'FOPT_PTeq-'+_Tb})

   # compute E4, E5, ELG
   EVALUATE({FORM:'FOPT_PTE-'+_Tb})

   PRINT_MEL({
       LIST:'ME-EPT4-'+_Tb,
       COMMENT:'EPT4 ('+_Tb+'):',
       FORMAT:'SCAL F24.14'})

   PRINT_MEL({
       LIST:'ME-EPT5-'+_Tb,
       COMMENT:'EPT5 ('+_Tb+'):',
       FORMAT:'SCAL F24.14'})

#MRCC-ENERGY for ease of use
PRINT_MEL({
    LIST:'MRCC_LAG_LST',
    COMMENT:'MRCCSD energy                           :',
    FORMAT:'SCAL F24.14'})

# DEF_FORMULA TEST   
EVALUATE({FORM:'FOPT_EPT4tot'})
PRINT_MEL({
    LIST:'ME-EPT4total',
    COMMENT:'Total 4th order triples correction      :',
    FORMAT:'SCAL F24.14'})

EVALUATE({FORM:'FOPT_EPT5tot'})
PRINT_MEL({
    LIST:'ME-EPT5total',
    COMMENT:'Total 5th order triples correction      :',
    FORMAT:'SCAL F24.14'})

# get elements from ME-PT, otherwise the formula will result in 0
if triples > 3:
    EVALUATE({FORM:'FOPT_ETRIP4_B'})
    EVALUATE({FORM:'FOPT_ETRIP5_B'})
    EVALUATE({FORM:'FOPT_Etot_B'})
    PRINT_MEL({
        LIST:'ME-Etotal_B',
        COMMENT:'Total icMRCCSD(T) {3} triples correction:',
        FORMAT:'SCAL F24.14'})
    EVALUATE({FORM:'FOPT_ETRIP4_E'})
    EVALUATE({FORM:'FOPT_ETRIP5_E'})
    EVALUATE({FORM:'FOPT_Etot_E'})
    PRINT_MEL({
        LIST:'ME-Etotal_E',
        COMMENT:'Total icMRCCSD(T) {4} triples correction:',
        FORMAT:'SCAL F24.14'})


if triples > 4:

    EVALUATE({FORM:'FOPT_ETRIP4_F'})
    EVALUATE({FORM:'FOPT_ETRIP5_F'})
    EVALUATE({FORM:'FOPT_Etot_F'})
    PRINT_MEL({
        LIST:'ME-Etotal_F',
        COMMENT:'Total icMRCCSD(T) {5} triples correction:',
        FORMAT:'SCAL F24.14'})

EVALUATE({FORM:'FOPT_EPTtot'})
PRINT_MEL({
    LIST:'ME-EPTtotal',
    COMMENT:'Total icMRCCSD(T) triples correction    :',
    FORMAT:'SCAL F24.14'})
      
EVALUATE({FORM:'FOPT_Etot'})
PRINT_MEL({
    LIST:'ME-Etotal',
    COMMENT:'Total icMRCCSD(T) energy                :',
    FORMAT:'SCAL F24.14'})
#EVALUATE({FORM:'FOPT_EPT4tot_test'})
#PRINT_MEL({
#    LIST:'ME-EPT4total_test',
#    COMMENT:'TEST Total 4. order triples correction:',
#    FORMAT:'SCAL F24.14'})
