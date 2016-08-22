
import sys,os
interface_path=os.path.join(os.getenv("GECCO_DIR"),"python_interface")
sys.path=[interface_path]+sys.path

from gecco_interface import *

_inp = GeCCo_Input()
_orb = Orb_Info()

# Get the informations about response calculation from the input
print '--- Starting Response Calculations ---'

_response_data={}
_order=_inp.get('method.MRCC.response.order')
if(_order == None):
    _order=1
_response_data['order']=int(_order)

_comp=_inp.get('method.MRCC.response.comp')
if(_comp == None):
    _comp='Z'
_response_data['comp']=_comp

_pert=_inp.get('method.MRCC.response.pert')
if(_pert == None):
    _pert='d'
_response_data['pert']=_pert

_freq=_inp.get('method.MRCC.response.freq')
if(_freq == None):
    _freq=0.0
_response_data['freq']=float(_freq)

_triplet=_inp.get('calculate.properties.triplet')
if(_triplet == None):
    _triplet=False
_response_data['triplet']=bool(_triplet)

print _response_data



export_targets();

print '--- Response Calculations are Done ---'
