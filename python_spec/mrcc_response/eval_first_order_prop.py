from gecco_interface import *
from get_response_data import _response_data

_triplet = _response_data['triplet']
_inp = GeCCo_Input()

if (_triplet):
    absym_dens = -1.0
else:
    absym_dens = 1.0

# Get the restart option to skip the calculation of lower order response properties
_restart=_inp.get('calculate.properties.restart')
if(_restart == None):
    _restart=1

##########
new_target('DENS_0')
DEF_OP_FROM_OCC({LABEL:'DENS_0',JOIN:2,
                 DESCR:',;,|H,;,H|H,P;,|H,V;,|,;P,H|,P;P,|,V;P,|,;V,H|,P;V,|,V;V,'})
##########

########## defining a dummy one body hamiltonian operator 
new_target('DIPM')
DEF_HAMILTONIAN({LABEL:'DIPM',MAX_RANK:1})
##########

########## generating the formula for density
new_target('F_DENS_0')

depend('DENS_0','DIPM')
depend('F_LMBD')

REPLACE({LABEL_RES:'F_LAG_DIP',
         LABEL_IN:'F_LAG',
         OP_LIST:['H','DIPM']})

INVARIANT({LABEL_RES:'F_LAG_DIP',
           LABEL_IN:'F_LAG_DIP',
           OP_RES:'RED_LAG',
           OPERATORS:'H'})


DERIVATIVE({LABEL_RES:'F_DENS_0',
            LABEL_IN:'F_LAG_DIP',
            OP_RES:'DENS_0',
            OP_DERIV:'DIPM'})

##########

##########
new_target('DEF_ME_DENS_0')

depend('DENS_0')

_op_list={'DENS_0':'ME_DENS_0'}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:absym_dens})
##########

##########
new_target('OPT_DENS')

depend('DEF_ME_DENS_0')
depend('F_DENS_0')
depend('LIST_LMBD')
depend('DEF_ME_C0','H0','DEF_ME_T')

OPTIMIZE({LABEL_OPT:'FOPT_DENS_0',
          LABELS_IN:'F_DENS_0'})
##########

########## evaluating the density
new_target('EVAL_DENS_0')

depend ('OPT_DENS')

EVALUATE({FORM:'FOPT_DENS_0'})

##########

########## calculating the first order properties
if (_restart<2):
    new_target('EVAL_DIPM',True)
else:
    new_target('EVAL_DIPM')

depend('EVAL_DENS_0')

_eval_prop_args={}
_eval_prop_args[DENS]='ME_DENS_0'
_eval_prop_args[RANK]=1
if (_triplet):
    _eval_prop_args[TRIPLET]=True
else:
    _eval_prop_args[TRIPLET]=False

EVAL_PROP(_eval_prop_args)
