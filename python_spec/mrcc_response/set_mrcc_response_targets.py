
import sys,os
interface_path=os.path.join(os.getenv("GECCO_DIR"),"python_interface")
sys.path=[interface_path]+sys.path

from gecco_interface import *
from get_response_data import _response_data

_inp = GeCCo_Input()
_orb = Orb_Info()


_env = _inp.env
print 'ENVIRONMENT:', _env

if (_response_data['triplet']):
    if(not ( _env=='DALTON')):
        quit_error('Properties for triplet perturbation are possible only using DALTON')


#import solve_lambda_0

export_targets();

print '--- Response Calculations are Done ---'
