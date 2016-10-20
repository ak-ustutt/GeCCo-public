
from python_interface.gecco_interface import *
from get_response_data import _response_data

_inp = GeCCo_Input()
_orb = Orb_Info()

# Get the name of the package GeCCo uses the integrals from 
_env = _inp.env

# For triplet perturbation, only DALTON provides the integrals, otherwise quit the calculation
if (_response_data['triplet']):
    if(not ( _env=='DALTON')):
        quit_error('Properties for triplet perturbation are possible only using DALTON')

# Get the restart option to skip the calculation of lower order response properties
_restart=_inp.get('calculate.properties.restart')
if(_restart == None):
    _restart=1

# check if we are relaxing only the cluster operators
# This would depend on options a) pure_vv:
_pure_vv=_inp.get('method.MR.pure_vv')
if(_pure_vv==None):
    _pure_vv=False
_pure_vv=bool(_pure_vv)

print 'pure_vv', _pure_vv
# b) optref:
_optref=_inp.get('calculate.solve.non_linear.optref')
if(_optref==None):
    _optref=-3

_optref=int(_optref)
print 'optref', _optref
# and c) nactel:
_nactel=_orb.get('nactel')

if(_pure_vv or (_optref==0) or (_nactel==1)):
    relax_ref=False
else:
    relax_ref=True

print 'relax_ref:', relax_ref

# first calculate the first order properties, need to solve the zeroeth order lambda parameters for that
if (_response_data['order']>=1):
    if (relax_ref):
        print 'importing solve_lambda_0 ...'
        import solve_lambda_0
        print 'importing eval_first_order_prop ...'
        import eval_first_order_prop
    else:
        print 'importing solve_lambda_0_restr ...'
        import solve_lambda_0_restr
        print 'importing eval_first_order_prop_restr ...'
        import eval_first_order_prop_restr

if (_response_data['order']>=2):
    if (relax_ref):
        print 'importing solve_t_1 ...'
        import solve_t_1
        print 'importing eval_second_order_prop ...'
        import eval_second_order_prop
    else:
        print 'importing solve_t_1_restr ...'
        import solve_t_1_restr
        print 'importing eval_second_order_prop_restr ...'
        import eval_second_order_prop_restr

if (_response_data['order']>=3):
    if (relax_ref):
        print 'importing solve_lambda_1 ...'
        import solve_lambda_1
        print 'importing eval_second_order_prop_alt ...'
        import eval_second_order_prop_alt
    else:
        print 'importing solve_lambda_0_restr ...'
        import solve_lambda_1_restr
        print 'importing eval_second_order_prop_alt_restr ...'
        import eval_second_order_prop_alt_restr

export_targets();
