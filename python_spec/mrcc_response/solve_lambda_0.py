from python_interface.gecco_interface import *
from python_spec.mrcc_response.get_response_data import _response_data

_inp = GeCCo_Input()
_orb = Orb_Info()

# Get all the necessary options needed during the calculation

_s2 = _orb.get('imult')

_ms = _orb.get('ims')

_isym = _orb.get('lsym')

_root = _inp.get('method.MR.ciroot')

_spinadapt = _inp.get('calculate.routes.spinadapt')
# unfortunately, GeCCo_Input does not provide the defaults
if _spinadapt == None:
  _spinadapt = 3 
  if _s2 == 1:
    _spinadapt = 0

_option=_response_data['option']

_eig_zero=_response_data['eig_zero']

_triplet=_response_data['triplet']

if ((_ms == 0) and (((_s2-1) % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and (((_s2+1) % 4) == 0)):
    _msc = -1
else:
    _msc = 0

if (_root == None):
    _root = 1
elif (int(_root) > 2):
    quit_error('Solution of Lambda equation not possible for ciroot greater than 2')
_root = int(_root)

#set the '2ms' value for the operators associated with the ci coefficients
if (_triplet or _ms!=0):
    _ms_c0= -1.0*_ms
else:
    _ms_c0 = _ms
    
# Get the restart option to skip the calculation of lower order response properties
_restart=_response_data['restart']


#Defining C0_bar separately from other operators used here. This is done as CO_bar will be needed outside.

new_target('C0_bar')
depend('C0')

CLONE_OPERATOR({LABEL:'C0_bar',TEMPLATE:'C0',ADJOINT:True})

new_target('LMBD_OP')

#defining the necessary operators

depend('OMG','L','Ltr','C0','D')

DEF_SCALAR({LABEL:'RED_LAG'})
DEF_SCALAR({LABEL:'RED_LAG_1'})
DEF_SCALAR({LABEL:'RED_LAG_2'})
DEF_SCALAR({LABEL:'RED_LAG_3'})
DEF_SCALAR({LABEL:'RED_LAG_4'})
DEF_SCALAR({LABEL:'C_NORM'})
CLONE_OPERATOR({LABEL:'OMG_L',TEMPLATE:'OMG',ADJOINT:True})
CLONE_OPERATOR({LABEL:'OMG_SL',TEMPLATE:'OMG',ADJOINT:True})
CLONE_OPERATOR({LABEL:'OMG_C0_bar',TEMPLATE:'C0',ADJOINT:True})
CLONE_OPERATOR({LABEL:'OMG_SC0',TEMPLATE:'C0',ADJOINT:True})
CLONE_OPERATOR({LABEL:'DIA_L',TEMPLATE:'L'})
#CLONE_OPERATOR({LABEL:'DIA_C0_bar',TEMPLATE:'C0',ADJOINT:True})
CLONE_OPERATOR({LABEL:'DIA_C0_bar',TEMPLATE:'C0',ADJOINT:False})
CLONE_OPERATOR({LABEL:'D_1',TEMPLATE:'D'})

new_target('F_LMBD')

depend('C0_bar','LMBD_OP','F_MRCC_LAG','F_OMG','F_OMG_C0','"E(MR)"','F_D','F_DENS0')

# first constructing <Psi_0|Lambda \bar{H}|Psi_0> from MRCC_LAG

DERIVATIVE({LABEL_RES:'F_LAG_1',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'RED_LAG_1',
            OP_DERIV:'L',
            OP_MULT:'L'})

# constructing <\bar{Psi_0}|\bar{H}|Psi_0> from MRCC_LAG

INVARIANT({LABEL_RES:'F_TEMP',
           LABEL_IN:'F_MRCC_LAG',
           OP_RES:'RED_LAG_4',
           OPERATORS:'L'})

DERIVATIVE({LABEL_RES:'F_LAG_2',
            LABEL_IN:'F_TEMP',
            OP_RES:'RED_LAG_2',
            OP_DERIV:'C0^+',
            OP_MULT:'C0_bar'}) 

if (_eig_zero):
# substracting E <\bar{Psi_0}|Psi_0> from F_LAG_2

    EXPAND_OP_PRODUCT({LABEL:'F_LAG_2',NEW:False,
                       OP_RES:'RED_LAG_2',
                       OPERATORS:['E(MR)','C0_bar','C0'],
                       IDX_SV:[1,2,3],
                       FAC:-1.0})

# adding F_LAG_1 and F_LAG_2 to get F_LAG

DEF_FORMULA({LABEL:'F_LAG',
             FORMULA:'RED_LAG=RED_LAG_1+RED_LAG_2'})

EXPAND({LABEL_RES:'F_LAG',
        LABEL_IN:'F_LAG',
        INTERM:['F_LAG_1','F_LAG_2']})


### adding the extra part in the lagrangian transposing the formula F_LAG_1

INVARIANT({LABEL_RES:'F_LAG_3',
           LABEL_IN:'F_LAG_1',
           OP_RES:'RED_LAG_3',
           OPERATORS:'Ltr'})

TRANSPS_FORMULA({LABEL_IN:'F_LAG_3',LABEL_RES:'F_LAG_3',
                 OP_RES:'RED_LAG_3',INIT:False,MULTI:True})


## adding F_LAG_3 to F_LAG for option=1

if (_option == 1):
    DEF_FORMULA({LABEL:'F_LAG_ext',
                 FORMULA:'RED_LAG=RED_LAG+RED_LAG_3'})
 
    EXPAND({LABEL_RES:'F_LAG_ext',
            LABEL_IN:'F_LAG_ext',
            INTERM:['F_LAG','F_LAG_3']})
elif (_option == 2):
    DEF_FORMULA({LABEL:'F_LAG_ext',
                 FORMULA:'RED_LAG=RED_LAG'})

    EXPAND({LABEL_RES:'F_LAG_ext',
            LABEL_IN:'F_LAG_ext',
            INTERM:['F_LAG']})

# Constructing A_{\lambda t}+A_{\bar{c} t} => F_L_1 from F_LAG

if (not _eig_zero):
    EXPAND_OP_PRODUCT({LABEL:'F_LAG',NEW:False,
                       OP_RES:'RED_LAG',
                       OPERATORS:['C0^+','E(MR)','L','T','C0'],
                       IDX_SV:[1,2,3,4,5],
                       FAC:1.0})

DERIVATIVE({LABEL_RES:'F_L_1',
            LABEL_IN:'F_LAG',
            OP_RES:'OMG_L',
            OP_DERIV:'T'})

# Constructing A_{\lambda \bar{c}}+A_{\bar{c} c} => F_L_2 from F_LAG_ext

DERIVATIVE({LABEL_RES:'F_L_2',
            LABEL_IN:'F_LAG_ext',
            OP_RES:'OMG_C0_bar',
            OP_DERIV:'C0'})

# obtaining formulas for metric vector products for L and C0_bar

EXPAND_OP_PRODUCT({LABEL:'F_LAG_4',NEW:True,
                   OP_RES:'RED_LAG',
                   OPERATORS:['RED_LAG','C0^+','L','T','C0','RED_LAG'],
                   IDX_SV:[1,2,3,4,5,1]})

#EXPAND_OP_PRODUCT({LABEL:'F_LAG_4',NEW:False,
#                   OP_RES:'RED_LAG',
#                   OPERATORS:['RED_LAG','C0^+','L','"E(MR)"','T','C0','RED_LAG'], ##### CAREFUL !!! Added extra term
#                   IDX_SV:[1,2,3,4,5,6,1], FAC:-1.0})    #### CAREFUL !!! Added extra term

DERIVATIVE({LABEL_RES:'F_SL',
            LABEL_IN:'F_LAG_4',
            OP_RES:'OMG_SL',
            OP_DERIV:'T'})

EXPAND_OP_PRODUCT({LABEL:'F_OMG_SC0',NEW:True,
                   OP_RES:'OMG_SC0',
                   OPERATORS:['OMG_SC0','C0_bar','OMG_SC0'],
                   IDX_SV:[1,2,1]})

#EXPAND_OP_PRODUCT({LABEL:'F_OMG_SC0',NEW:False,
#                   OP_RES:'OMG_SC0',
#                   OPERATORS:['OMG_SC0','C0_bar','"E(MR)"','OMG_SC0'],  ##### CAREFUL !!! Added extra term
#                   IDX_SV:[1,2,3,1], FAC:-1.0})  ##### CAREFUL !!! Added extra term

#PRINT_FORMULA({LABEL:'F_DENS0'})

REPLACE({LABEL_RES:'F_DENS(1)',
         LABEL_IN:'F_DENS0',
         OP_LIST:['C0^+','C0_bar','C0','C0_bar^+']})

#PRINT_FORMULA({LABEL:'F_DENS(1)'})

#for spin projection
new_target('FOPT_C0_bar_sp')
depend('DEF_ME_C0_bar','DEF_ME_S+','DEF_ME_S-','DEF_ME_Sz')

CLONE_OPERATOR({LABEL:'C0_bar_sp',TEMPLATE:'C0_bar'})
_op_list={'C0_bar_sp':'ME_C0_bar_sp'}
for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:_isym,
                 '2MS':_ms_c0,
                 AB_SYM:_msc})
# FIXME: in future we should set s2 of C0_bar_sp here (check side effects first)
#        currently, the imult form orb_info is read
    # <C0bar| 1/2*(S+S- + S-S+) + Sz^2
    # <C0bar| {S+S-} + 1/2*(S+S-)_c + 1/2*(S-S+)_c + Sz^2
    # I admit I still don't know, why the upper does not work directly
    EXPAND_OP_PRODUCT({LABEL:'F_C0_bar_sp',NEW:True,OP_RES:'C0_bar_sp',
           FAC:1.0,OPERATORS:['C0_bar_sp','C0_bar','S+','S-','C0_bar_sp'],
                   IDX_SV:[1,2,3,4,1],AVOID:[3,4]})
    EXPAND_OP_PRODUCT({LABEL:'F_C0_bar_sp',NEW:False,OP_RES:'C0_bar_sp',
           FAC:0.5,OPERATORS:['C0_bar_sp','C0_bar','S+','S-','C0_bar_sp'],
                   IDX_SV:[1,2,3,4,1],CONNECT:[3,4]})
    EXPAND_OP_PRODUCT({LABEL:'F_C0_bar_sp',NEW:False,OP_RES:'C0_bar_sp',
           FAC:0.5,OPERATORS:['C0_bar_sp','C0_bar','S-','S+','C0_bar_sp'],
                   IDX_SV:[1,2,3,4,1],CONNECT:[3,4]})
    EXPAND_OP_PRODUCT({LABEL:'F_C0_bar_sp',NEW:False,OP_RES:'C0_bar_sp',
           FAC:1.0,OPERATORS:['C0_bar_sp','C0_bar','Sz','Sz','C0_bar_sp'],
                   IDX_SV:[1,2,3,4,1],FIX_VTX:True})
    OPTIMIZE({LABEL_OPT:'FOPT_C0_bar_sp',
              LABELS_IN:'F_C0_bar_sp'})
    

#some intermediates for fast calculations

new_target('F_prePPlint')

depend('F_LMBD',
       'DEF_ME_INT_PP',
       'H_PP')

REPLACE({LABEL_RES:'F_prePPlint',
         LABEL_IN:'F_L_1',
         OP_LIST:['H','H_PP']})

INVARIANT({LABEL_RES:'F_prePPlint',
           LABEL_IN:'F_prePPlint',
           OP_RES:'OMG_L',
           OPERATORS:'H'})

new_target('F_PPlint')

depend('F_prePPlint')

CLONE_OPERATOR({LABEL:'INT_PPl',
                TEMPLATE:'INT_PP',ADJOINT:True})

DERIVATIVE({LABEL_RES:'F_PPlint',
           LABEL_IN:'F_prePPlint',
           OP_RES:'INT_PPl',
           OP_DERIV:'H_PP'})


#creating ME list for L:

new_target('DEF_ME_L')
depend('L')

_s=-1
if _spinadapt > 1:
    _s=0

_op_list={'L':'ME_L'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,S2:_s,
                 AB_SYM:1,
                 MIN_REC:1,MAX_REC:_root})

#creating ME list for C0_bar:
new_target('DEF_ME_C0_bar')
depend('C0_bar')

_op_list={'C0_bar':'ME_C0_bar'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:_isym,
                 '2MS':_ms_c0,
                 AB_SYM:_msc,
                 MIN_REC:1,MAX_REC:_root})
# FIXME: in future we should set s2 of C0_bar here (check side effects first)
#        currently, the imult form orb_info is read

# creating ME list for the new operators:

new_target('LIST_LMBD')

depend('F_LMBD','F_L','DEF_ME_E(MR)')
depend('F_PPlint')

_op_list={'OMG_L':'ME_OMG_L',
          'OMG_SL':'ME_OMG_SL'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,S2:_s,
                 AB_SYM:1})

_op_list={'INT_PPl':'ME_INT_PPl'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:1})

_op_list={'OMG_C0_bar':'ME_OMG_C0_bar',
          'OMG_SC0':'ME_OMG_SC0'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:_isym,
                 '2MS':_ms_c0,
                 AB_SYM:_msc})

_op_list={'Ltr':'ME_Ltr'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,
                 MIN_REC:1,MAX_REC:_root})

_op_list={'DIA_L':'ME_DIA_L'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

_op_list={'DIA_C0_bar':'ME_DIA_C0_bar'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:_isym,
                 '2MS':_ms})
                 #'2MS':_ms_c0})

DEF_ME_LIST({LIST:'ME_MINEN',
             OPERATOR:'E(MR)',
             IRREP:1,
             '2MS':0})

ASSIGN_ME2OP({LIST:'ME_E(MR)',OPERATOR:'E(MR)'})

# Preconditioner for L

new_target('DIAG_L')

depend('EVAL_FREF','FOPT_Atr','DIAG1SxxM00_T','EVAL_Atr')

COPY_LIST({LIST_RES:'ME_DIA_L',LIST_INP:'DIAG1SxxM00_T',ADJOINT:True,FAC:1.0})

# Preconditioner for \bar{c}

new_target('DIAG_C0_bar')

depend('H0','LIST_LMBD')

PRECONDITIONER({LIST_PRC:'ME_DIA_C0_bar',
                LIST_INP:'H0',
                MODE:'dia-H'})

SCALE_COPY({LIST_RES:'ME_MINEN',
            LIST_INP:'ME_E(MR)',
            FAC:-1.0})

EXTRACT_DIAG({LIST_RES:'ME_DIA_C0_bar',
              LIST_IN:'ME_MINEN',
              MODE:'ext_act'})

PRINT_MEL({LIST:'ME_E(MR)',COMMENT:'ENERGY ICMRCC:'})

# Optimizing all the five formulas together for solve_evp

new_target('LMBD_OPT')

depend('DEF_ME_L','DEF_ME_C0_bar','DEF_ME_INT_D')
depend('LIST_LMBD','DEF_ME_C0','H0','DEF_ME_T','DIAG_L','DIAG_C0_bar','DEF_ME_Ttr','F_L','F_INT_D')

OPTIMIZE({LABEL_OPT:'LMBD_OPT',
          LABELS_IN:['F_L_1','F_L_2','F_SL','F_L','F_OMG_SC0'],
          INTERM:['F_PPlint','F_INT_D']})

OPTIMIZE({LABEL_OPT:'FOPT_L',
          LABELS_IN:'F_L'})

# Soving the equation using SOLVE_EVP

new_target('SOLVE_LMBD')

depend('LMBD_OPT','DEF_ME_Dtrdag','DEF_ME_Ttr')

if (_root == 1):
    SET_MEL({LIST:'ME_L',IDX_LIST:1,VAL_LIST:0.0})
    COPY_LIST({LIST_RES:'ME_C0_bar',LIST_INP:'ME_C0',ADJOINT:True})
elif (_root == 2):
    SET_STATE({LISTS:'ME_C0',ISTATE:1})
    SET_STATE({LISTS:'ME_C0_bar',ISTATE:1})
    COPY_LIST({LIST_RES:'ME_C0_bar',LIST_INP:'ME_C0',ADJOINT:True})
    SET_STATE({LISTS:'ME_L',ISTATE:1})
    SET_MEL({LIST:'ME_L',IDX_LIST:1,VAL_LIST:0.0})
    SET_STATE({LISTS:'ME_C0',ISTATE:2})
    SET_STATE({LISTS:'ME_C0_bar',ISTATE:2})
    COPY_LIST({LIST_RES:'ME_C0_bar',LIST_INP:'ME_C0',ADJOINT:True})
    SET_STATE({LISTS:'ME_L',ISTATE:2})
    SET_MEL({LIST:'ME_L',IDX_LIST:1,VAL_LIST:0.0})

_solve_evp_basis={}
_solve_evp_basis[LIST_OPT]=['ME_L','ME_C0_bar']
_solve_evp_basis[LIST_PRC]=['ME_DIA_L','ME_DIA_C0_bar']
_solve_evp_basis[OP_MVP]=['OMG_L','OMG_C0_bar']
_solve_evp_basis[OP_SVP]=['OMG_SL','OMG_SC0']
_solve_evp_basis[FORM]='LMBD_OPT'
_solve_evp_basis[INIT]=True
_solve_evp_basis[N_ROOTS]=_root
if _spinadapt == 0:
  _solve_evp_basis[MODE]='TRF DIA'
  _solve_evp_basis[LIST_SPC]=['ME_Ltr','ME_Dtr','ME_Dtrdag']
else:
  depend('FOPT_C0_bar_sp')
  _solve_evp_basis[MODE]='TRF SPP'
  _solve_evp_basis[LIST_SPC]=['ME_Ltr','ME_Dtr','ME_Dtrdag','ME_C0_bar_sp']
  _solve_evp_basis[FORM_SPC]=['FOPT_C0_bar_sp']

SOLVE_EVP(_solve_evp_basis)

## Now renormalizing L and C0_bar vectors so that C0_bar remains normalized respect to C0

if (_restart<2):
    new_target('NORMALIZE',True)
else:
    new_target('NORMALIZE')

depend('SOLVE_LMBD')

EXPAND_OP_PRODUCT({LABEL:'F_NORM_C',NEW:True,
                   OP_RES:'C_NORM',
                   OPERATORS:['C0_bar','C0'],
                   IDX_SV:[1,2],
                   FAC:1.0})

DEF_ME_LIST({LIST:'ME_C_NORM',
             OPERATOR:'C_NORM',
             IRREP:1,
             '2MS':0})
OPTIMIZE({LABEL_OPT:'FOPT_NORM_C',
          LABELS_IN:'F_NORM_C'})
EVALUATE({FORM:'FOPT_NORM_C'})

SCALE({LIST_RES:'ME_C0_bar',LIST_INP:'ME_C0_bar',
       LIST_SCAL:'ME_C_NORM',FAC:1.0,INV:True})

SCALE({LIST_RES:'ME_L',LIST_INP:'ME_L',
       LIST_SCAL:'ME_C_NORM',FAC:1.0,INV:True})
