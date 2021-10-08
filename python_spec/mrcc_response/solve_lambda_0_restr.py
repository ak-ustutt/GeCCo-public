from python_interface.gecco_interface import *
from python_spec.mrcc_response.get_response_data import _response_data

_inp = GeCCo_Input()
_orb = Orb_Info()

_s2 = _orb.get('imult')

_ms = _orb.get('ims')

_isym = _orb.get('lsym')

_triplet=_response_data['triplet']

# unfortunately, GeCCo_Input does not provide the defaults
_spinadapt = 3 #_inp.get('calculate.routes.spinadapt')
if _s2 == 1:
  _spinadapt = 0

if ((_ms == 0) and ((_s2-1 % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and ((_s2+1 % 4) == 0)):
    _msc = -1
else:
    _msc = 0

# Get the restart option to skip the calculation of lower order response properties
_restart=_inp.get('calculate.properties.restart')
if(_restart == None):
    _restart=1


new_target('LMBD_OP')

depend('OMG','L','Ltr')

CLONE_OPERATOR({LABEL:'L_res',TEMPLATE:'OMG',ADJOINT:True})
CLONE_OPERATOR({LABEL:'L_res_rhs',TEMPLATE:'OMG',ADJOINT:True})
CLONE_OPERATOR({LABEL:'L_res_trf',TEMPLATE:'OMG',ADJOINT:True})
CLONE_OPERATOR({LABEL:'DIA_L',TEMPLATE:'L'})
DEF_SCALAR({LABEL:'RED_LAG'})

new_target('F_LMBD')

depend('LMBD_OP','T','F_T','F_MRCC_LAG','D','Dtr')

DERIVATIVE({LABEL_RES:'F_lres',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'L_res',
            OP_DERIV:'T'})

LEQ_SPLIT({LABEL_RAW:'F_lres',
            LABEL_RHS:'F_lres_rhs',
            LABEL_TRF:'F_lres_trf',
            OP_X:'L',
            OP_RHS:'L_res_rhs',
            OP_TRF:'L_res_trf'})

#PRINT_FORMULA({LABEL:'F_lres_rhs'})
#PRINT_FORMULA({LABEL:'F_lres_trf'})

#creating ME list for L:

new_target('DEF_ME_L')
depend('L')

_s = -1
if _spinadapt > 1:
  _s = 0

_op_list={'L':'ME_L'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,S2:_s,
                 AB_SYM:1})

new_target('LIST_LMBD')

depend('F_LMBD','F_L')

_op_list={'L_res_rhs':'ME_Lresr',
          'L_res_trf':'ME_Lrest'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,S2:_s,
                 AB_SYM:1})

_op_list={'Ltr':'ME_Ltr',
          'DIA_L':'ME_DIA_L'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0})

new_target('DIAG_L')

depend('EVAL_FREF','FOPT_Atr','DIAG1SxxM00_T','EVAL_Atr')

COPY_LIST({LIST_RES:'ME_DIA_L',LIST_INP:'DIAG1SxxM00_T',ADJOINT:True,FAC:1.0})

new_target('LMBD_OPT')

depend('DEF_ME_L','DEF_ME_INT_D')
depend('LIST_LMBD','H0','DEF_ME_T','DIAG_L','DEF_ME_Ttr','F_L','F_INT_D')

#PRINT_FORMULA({LABEL:'F_L'})
OPTIMIZE({LABEL_OPT:'LMBD_OPT',
          LABELS_IN:['F_lres_rhs','F_lres_trf'],
          INTERM:['F_INT_D']})

OPTIMIZE({LABEL_OPT:'FOPT_L',
          LABELS_IN:'F_L'})

if (_restart<2):
    new_target('SOLVE_LMBD',True)
else:
    new_target('SOLVE_LMBD')

depend('LMBD_OPT','DEF_ME_Dtrdag','DEF_ME_Ttr')

SOLVE_LEQ({LIST_OPT:'ME_L',
           LIST_PRC:'ME_DIA_L',
           OP_MVP:'L_res_trf',
           OP_SVP:'L',
           OP_RHS:'L_res_rhs',
           FORM:'LMBD_OPT',
           LIST_SPC:['ME_L','ME_Ltr','ME_Dtr','ME_Dtrdag'],
           FORM_SPC:'FOPT_L',
           MODE:'TRF',
           N_ROOTS:1})
