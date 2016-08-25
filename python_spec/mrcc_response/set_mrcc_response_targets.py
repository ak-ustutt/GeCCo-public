
import sys,os
interface_path=os.path.join(os.getenv("GECCO_DIR"),"python_interface")
sys.path=[interface_path]+sys.path

from gecco_interface import *
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

# b) optref:
_optref=_inp.get('calculate.solve.optref')
if(_optref==None):
    _optref=-3

# and c) nactel:
_nactel=_orb.get('nactel')

if(_pure_vv or (_optref==0) or (_nactel==1)):
    relax_ref=False
else:
    relax_ref=True

# first calculate the first order properties, need to solve the zeroeth order lambda parameters for that
if (_response_data['order']>=1):
    if (relax_ref):
        import solve_lambda_0
        import eval_first_order_prop
    else:
        import solve_lambda_0_restr
        import eval_first_order_prop_restr

export_targets();
