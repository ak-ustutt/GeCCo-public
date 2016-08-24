# Get the informations about response calculation from the input
from gecco_interface import *

_inp = GeCCo_Input()

_response_data = {}
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

_use_b=_inp.get('method.MRCC.response.use_b')
if(_use_b == None):
    _use_b=True

if(_use_b):
    option=1
else:
    option=2

_response_data['option']=bool(option)

print 'Using following options for the response calculations:'
print _response_data
