from python_interface.gecco_interface import *
from get_response_data import _response_data
from set_mrcc_response_targets import relax_ref

_triplet = _response_data['triplet']
_inp = GeCCo_Input()

# Get the name of the package GeCCo uses the integrals from 
_env = _inp.env

# Get the restart option to skip the calculation of lower order response properties
_restart=_response_data['restart']

##########
new_target('DENS_0')
DEF_OP_FROM_OCC({LABEL:'DENS_0S',JOIN:2,
                 DESCR:',;,|H,;,H|H,P;,|H,V;,|,;P,H|,P;P,|,V;P,|,;V,H|,P;V,|,V;V,'})
DEF_OP_FROM_OCC({LABEL:'DENS_0T',JOIN:2,
                 DESCR:',;,|H,;,H|H,P;,|H,V;,|,;P,H|,P;P,|,V;P,|,;V,H|,P;V,|,V;V,'})
##########

########## defining a dummy one body hamiltonian operator 
new_target('DIPM')
DEF_HAMILTONIAN({LABEL:'DIPM',MAX_RANK:1})
##########

########## generating the formula for density
new_target('F_DENS_0')

if(relax_ref):
    _lag='F_LAG'
    _depend='F_LMBD'
else:
    _lag='F_MRCC_LAG'
    _depend='F_MRCC_LAG'

depend('DENS_0','DIPM')
depend(_depend)

depend('F_LMBD')

REPLACE({LABEL_RES:'F_LAG_DIP',
         LABEL_IN:_lag,
         OP_LIST:['H','DIPM']})

INVARIANT({LABEL_RES:'F_LAG_DIP',
           LABEL_IN:'F_LAG_DIP',
           OP_RES:'RED_LAG',
           OPERATORS:'H'})


DERIVATIVE({LABEL_RES:'F_DENS_0S',
            LABEL_IN:'F_LAG_DIP',
            OP_RES:'DENS_0S',
            OP_DERIV:'DIPM'})

DERIVATIVE({LABEL_RES:'F_DENS_0T',
            LABEL_IN:'F_LAG_DIP',
            OP_RES:'DENS_0T',
            OP_DERIV:'DIPM'})

##########

##########
new_target('DEF_ME_DENS_0')

depend('DENS_0')

_op_list={'DENS_0S':'ME_DENS_0S',
          'DENS_0T':'ME_DENS_0T'}
_ab_list={'DENS_0S':+1,
          'DENS_0T':-1}

for _op in _op_list:
    DEF_ME_LIST({LIST:_op_list[_op],
                 OPERATOR:_op,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:_ab_list[_op]})
##########

##########
new_target('OPT_DENS')

depend('DEF_ME_DENS_0')
depend('F_DENS_0')
depend('LIST_LMBD')
depend('DEF_ME_C0','H0','DEF_ME_T')

OPTIMIZE({LABEL_OPT:'FOPT_DENS_0',
          LABELS_IN:['F_DENS_0S','F_DENS_0T']})
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
_eval_prop_args[RANK]=1
_eval_prop_args[DENS]='ME_DENS_0S'
_eval_prop_args[TRIPLET]=False

PRINT({STRING:'Evaluate properties from singlet density'})
EVAL_PROP(_eval_prop_args)

_eval_prop_args={}
_eval_prop_args[RANK]=1
_eval_prop_args[DENS]='ME_DENS_0T'
_eval_prop_args[TRIPLET]=True

PRINT({STRING:'Evaluate properties from triplet density'})
EVAL_PROP(_eval_prop_args)
