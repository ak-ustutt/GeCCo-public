"""Solve/evaluate the icMRCC (T) equations

History:

Andreas Nov 2021: Creation of initial stub based on some earlier experiments

"""
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

verbosity=100

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


new_target('SOLVE_MRCC_PT')
depend('SOLVE_MRCC')
heading('Computing the (T) correction')

# set the various occupation classes
# might be done more elegantly, here we just enumerate everything:
Tblocks={}
Tblocks['PPP-HHH']={}
Tblocks['PPP-HHH']['occ']='PPP,HHH'
Tblocks['PPP-HHH']['res']=',;PPP,HHH'
Tblocks['PPP-HHH']['nact']=0

Tblocks['PPP-HHV']={}
Tblocks['PPP-HHV']['occ']='PPP,HHV'
Tblocks['PPP-HHV']['res']=',V;PPP,HH'
Tblocks['PPP-HHV']['nact']=1

Tblocks['PPV-HHH']={}
Tblocks['PPV-HHH']['occ']='PPV,HHH'
Tblocks['PPV-HHH']['res']=',;PPV,HHH'
Tblocks['PPV-HHH']['nact']=1

Tblocks['PPP-HVV']={}
Tblocks['PPP-HVV']['occ']='PPP,HVV'
Tblocks['PPP-HVV']['res']=',VV;PPP,H'
Tblocks['PPP-HVV']['nact']=2

Tblocks['PVV-HHH']={}
Tblocks['PVV-HHH']['occ']='PVV,HHH'
Tblocks['PVV-HHH']['res']=',;PVV,HHH'
Tblocks['PVV-HHH']['nact']=2

Tblocks['PPV-HHV']={}
Tblocks['PPV-HHV']['occ']='PPV,HHV|PP,HH'
Tblocks['PPV-HHV']['res']=',V;PPV,HH|,;PP,HH'
Tblocks['PPV-HHV']['nact']=2

# for cas(2,2) we can skip these ...
#Tblocks['PPP-VVV']={}
#Tblocks['PPP-VVV']['occ']='PPP,VVV'
#Tblocks['PPP-VVV']['res']=',VVV;PPP,'
#Tblocks['PPP-VVV']['nact']=3

#Tblocks['VVV-HHH']={}
#Tblocks['VVV-HHH']['occ']='VVV,HHH'
#Tblocks['VVV-HHH']['res']=',;VVV,HHH'
#Tblocks['VVV-HHH']['nact']=3

#Tblocks['PPV-HVV']={}
#Tblocks['PPV-HVV']['occ']='PPV,HVV|PP,HV'
#Tblocks['PPV-HVV']['res']=',VV;PPV,H|,V;PP,H'
#Tblocks['PPV-HVV']['nact']=3

#Tblocks['PVV-HHV']={}
#Tblocks['PVV-HHV']['occ']='PVV,HHV|PV,HH'
#Tblocks['PVV-HHV']['res']=',V;PVV,HH|,;PV,HH'
#Tblocks['PVV-HHV']['nact']=3

# TODO: also consider here the cases PVV-HVV VVV-HHV PPV-VVV and maybe PVV-VVV and VVV-VVH

# generat some operators to address the individual T3 blocks:
for _Tb in Tblocks:

   DEF_OP_FROM_OCC({LABEL:'T3-'+_Tb,DESCR:Tblocks[_Tb]['occ']})
   DEF_ME_LIST({LIST:'ME-T3-'+_Tb,OPERATOR:'T3-'+_Tb,IRREP:1,'2MS':_ms,'S2':_s2,'ABSYM':_msc})
   # maybe: use the T3g list instead? will require a bit of change to the formulae below, but I think it is doable ...

   # transformed T3 ... will probably need this
   CLONE_OPERATOR({LABEL:'T3tr-'+_Tb,TEMPLATE:'T3-'+_Tb})
   DEF_ME_LIST({LIST:'ME-T3tr-'+_Tb,OPERATOR:'T3tr-'+_Tb,IRREP:1,'2MS':_ms,'ABSYM':_msc})

   # create own diagonal?
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
   debug_FORM('F_PTrhs_'+_Tb,only_this=True)
   debug_FORM('F_PTtrf_'+_Tb,only_this=True)
   debug_FORM('F_PT_E4_'+_Tb,only_this=True)
   debug_FORM('F_PT_E5_'+_Tb,only_this=True)
   debug_FORM('F_PT_LG_'+_Tb,only_this=True)

   # ... and the optimized formulas
   OPTIMIZE({LABEL_OPT:'FOPT_PTeq-'+_Tb,
             LABELS_IN:['F_PTrhs_'+_Tb,'F_PTtrf_'+_Tb]})
   OPTIMIZE({LABEL_OPT:'FOPT_PTE-'+_Tb,
             LABELS_IN:['F_PT_E4_'+_Tb,'F_PT_E5_'+_Tb,'F_PT_LG_'+_Tb]})

   if Tblocks[_Tb]["nact"]>0:
      # and a transformation (that keeps T3 orth to T2 and T1)
      EXPAND_OP_PRODUCT({LABEL:'F_T3tr-'+_Tb,NEW:True,
                      OP_RES:'T3-'+_Tb,
                      OPERATORS:['T3-'+_Tb,'X_TRM','T3tr-'+_Tb,'X_TRM','T3-'+_Tb],
                      IDX_SV:   [1,2,3,2,1],
                      AVOID:[2,4,1,4,2,5]})
      SELECT_LINE({
        LABEL_IN:'F_T3tr-'+_Tb,
        LABEL_RES:'F_T3tr-'+_Tb,
        OP_RES:'T3-'+_Tb,
        OP_INCL:'T3tr-'+_Tb,
        IGAST:3,
        MODE:'no_ext'})
      debug_FORM('F_T3tr-'+_Tb,only_this=True)
      SELECT_SPECIAL({LABEL_RES:'F_T3tr-'+_Tb,LABEL_IN:'F_T3tr-'+_Tb,
             TYPE:'rank',MODE:'33',OPERATORS:['T3-'+_Tb,'T3tr-'+_Tb]})
      debug_FORM('F_T3tr-'+_Tb,only_this=True)
      OPTIMIZE({LABEL_OPT:'FOPT_T3tr-'+_Tb,
             LABELS_IN:'F_T3tr-'+_Tb})

# loop over blocks of T3
for _Tb in Tblocks:
   # set up preconditioner
   PRECONDITIONER({LIST_PRC:'ME-DIA-'+_Tb,LIST_INP:'FOCK_EFF_INACT_LST',MODE:'dia-F'})

   if Tblocks[_Tb]['nact']>0:
      EXTRACT_DIAG({LIST_RES:'ME-DIA-'+_Tb,LIST_IN:'A_TRF_LST',MODE:'extend'})

   # solve linear equations
   PRINT({STRING:'Now solving (T) equations for block '+Tblocks[_Tb]['occ']})
   if Tblocks[_Tb]['nact']>0:
      SOLVE_LEQ({LIST_OPT:'ME-T3-'+_Tb,MODE:'TRF',
              OP_MVP:'PTtrf-'+_Tb,OP_RHS:'PTrhs-'+_Tb,
              OP_SVP:'T3-'+_Tb,N_ROOTS:1,
              LIST_PRC:'ME-DIA-'+_Tb,
              FORM:'FOPT_PTeq-'+_Tb,
              LIST_SPC:['ME-T3-'+_Tb,'ME-T3tr-'+_Tb,'ME_X_TRM','ME_X_TRM_DAG'],
              FORM_SPC:'FOPT_T3tr-'+_Tb})
   else:
      # simpler case here (consider this non-iteratively)
      SOLVE_LEQ({LIST_OPT:'ME-T3-'+_Tb,MODE:'DIA',
              OP_MVP:'PTtrf-'+_Tb,OP_RHS:'PTrhs-'+_Tb,
              OP_SVP:'T3-'+_Tb,N_ROOTS:1,
              LIST_PRC:'ME-DIA-'+_Tb,
              FORM:'FOPT_PTeq-'+_Tb})

   # compute E4, E5, ELG
   EVALUATE({FORM:'FOPT_PTE-'+_Tb})

      
