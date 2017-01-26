
from python_interface.gecco_interface import *
from get_response_data import _response_data, _pop_data, _cmp_data, _calc_data
from set_mrcc_response_targets import relax_ref
import math

_inp = GeCCo_Input()

# Get the name of the package GeCCo uses the integrals from 
_env = _inp.env

_orb = Orb_Info()

_s2 = _orb.get('imult')

_ms = _orb.get('ims')

_isym_0 = _orb.get('lsym')
_isym= _isym_0

if ((_ms == 0) and ((_s2-1 % 4) == 0)):
    _msc = 1
elif ((_ms == 0) and ((_s2+1 % 4) == 0)):
    _msc = -1
else:
    _msc = 0

#Getting the frequency
_freq=_response_data['freq']
#Getting the value of the restart option
_restart=_response_data['restart']
#Getting the total number of perturbation operator need to be defined 
_npop=_response_data['nPop']
#Getting total number of response calculation
_ncnt=_response_data['nCnt']
#Getting the maximum order of the properties that will be calculated
_maxord=_response_data['maxorder']
#n_par tells how many version of the same operator has to be defined 

#Defining all the perturbation operators needed for the whole calculations

_list_to_depend=[]

for ipop in xrange (0,_npop):
    _cur_ext=_pop_data['name'][ipop]+_pop_data['comp'][ipop]
    _pop_name='V'+_cur_ext

    new_target(_pop_name)

    DEF_HAMILTONIAN({LABEL:_pop_name,MAX_RANK:1})

    new_target('IMPORT_'+_pop_name)

    depend(_pop_name)

    _op_list={_pop_name:'ME_'+_pop_name}

    _irrep=_pop_data['isym'][ipop]
    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_irrep,
                     '2MS':0,
                     AB_SYM:1}) # For second order properties, the perturbations are mainly singlet. 

#Getting the type of the integrals needed while importing the integrals
    _int_type=_pop_data['int_name'][ipop]

    IMPORT({LIST:'ME_'+_pop_name,
            TYPE:_int_type,
            ENV:_env})

    PRINT_MEL({LIST:'ME_'+_pop_name})
    _list_to_depend.append('IMPORT_'+_pop_name)

    new_target('DEF_RSPNS(1)'+_cur_ext)

    DEF_SCALAR({LABEL:'RSPNS(1)'+_cur_ext})

    _op_list={'RSPNS(1)'+_cur_ext:'ME_RSPNS(1)'+_cur_ext}

    for _op in _op_list:
        DEF_ME_LIST({LIST:_op_list[_op],
                     OPERATOR:_op,
                     IRREP:_irrep,
                     '2MS':0})

new_target('IMPORT_PERT_OPS')

for ele in _list_to_depend:
    depend(ele)

# We calculate the first order response, again(!), that will be needed here
# Here we will first get the formula, and then for each of the operator we 
# will use the same formula, but changing only the list of the operators

new_target('GET_RSPNS(1)_FORM')

depend('F_MRCC_LAG')
depend('DEF_ME_T')
depend('DEF_ME_C0')
depend('DEF_ME_L')
if (relax_ref):
	depend('DEF_ME_C0_bar')

DEF_SCALAR({LABEL:'preRSPNS(1)_1'})
DEF_SCALAR({LABEL:'preRSPNS(1)_2'})
DEF_SCALAR({LABEL:'preRSPNS(1)'})

#Defining a dummy scalar operator for RSPNS(1)
DEF_SCALAR({LABEL:'RSPNS(1)'})

#Defining a dummy perturbation operator
DEF_HAMILTONIAN({LABEL:'V(1)',MAX_RANK:1})

INVARIANT({LABEL_RES:'F_preRSPNS(1)_1',
           LABEL_IN:'F_MRCC_LAG',
           OP_RES:'preRSPNS(1)_1',
           OPERATORS:'L'})

if (relax_ref):
	REPLACE({LABEL_RES:'F_preRSPNS(1)_1',
    	     LABEL_IN:'F_preRSPNS(1)_1',
        	 OP_LIST:['C0^+','C0_bar']})

DERIVATIVE({LABEL_RES:'F_preRSPNS(1)_2',
            LABEL_IN:'F_MRCC_LAG',
            OP_RES:'preRSPNS(1)_2',
            OP_DERIV:'L',
            OP_MULT:'L'})

DEF_FORMULA({LABEL:'F_preRSPNS(1)',
             FORMULA:'preRSPNS(1)=preRSPNS(1)_1+preRSPNS(1)_2'})
#            FORMULA:'preRSPNS(1)=preRSPNS(1)_2'})

EXPAND({LABEL_RES:'F_preRSPNS(1)',
        LABEL_IN:'F_preRSPNS(1)',
        INTERM:['F_preRSPNS(1)_1','F_preRSPNS(1)_2']})

REPLACE({LABEL_RES:'F_RSPNS(1)',
         LABEL_IN:'F_preRSPNS(1)',
         OP_LIST:['H','V(1)']})

INVARIANT({LABEL_RES:'F_RSPNS(1)',
           LABEL_IN:'F_RSPNS(1)',
           OP_RES:'RSPNS(1)',
           OPERATORS:'H'})

_list_to_depend=[]

for ipop in xrange (0,_npop):
    _cur_ext=_pop_data['name'][ipop]+_pop_data['comp'][ipop]
    _pop_name='V'+_cur_ext

    new_target('EVAL_RSPNS(1)'+_cur_ext)

    depend('DEF_RSPNS(1)'+_cur_ext)

    ASSIGN_ME2OP({LIST:'ME_RSPNS(1)'+_cur_ext,
                  OPERATOR:'RSPNS(1)'})

    ASSIGN_ME2OP({LIST:'ME_'+_pop_name,
                  OPERATOR:'V(1)'})

    depend('H0')
    depend('GET_RSPNS(1)_FORM')

    OPTIMIZE({LABEL_OPT:'FOPT_RSPNS(1)'+_cur_ext,
              LABELS_IN:'F_RSPNS(1)'})

    EVALUATE({FORM:'FOPT_RSPNS(1)'+_cur_ext})

    PRINT_MEL({LIST:'ME_RSPNS(1)'+_cur_ext,
#              FORMAT:'SCAL',
               COMMENT:'First order property('+str(ipop)+')'})

    _list_to_depend.append('EVAL_RSPNS(1)'+_cur_ext)

########## calculating the first order properties
if (_restart<2):
    new_target('EVAL_RSPNS(1)',True)
else:
    new_target('EVAL_RSPNS(1)')

depend('IMPORT_PERT_OPS')
for ele in _list_to_depend:
    depend(ele)

