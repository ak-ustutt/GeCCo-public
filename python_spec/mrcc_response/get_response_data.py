# Get the informations about response calculation from the input
from gecco_interface import *

_inp = GeCCo_Input()

#Get the maximum order of the properties we are going to calculate
_response_data = {}
_order=_inp.get('method.MRCC.response.order')
if(_order == None):
    _order=1
_response_data['order']=int(_order)

#Get the type of the perturbation, the default is to use the dipole operator
_pert=_inp.get('method.MRCC.response.pert')
if(_pert == None):
    _pert='d'
_response_data['pert']=_pert

#Get the cartesian component of the perturbation
_comp=_inp.get('method.MRCC.response.comp')
if(_comp == None):
    _comp='Z'
_response_data['comp']=_comp

#Get the frequency of the perturbation if we are going to calculate dynamic properties
_freq=_inp.get('method.MRCC.response.freq')
if(_freq == None):
    _freq=0.0
_response_data['freq']=float(_freq)

#Check if the perturation is triplet, that makes the resultant property a spin-dependent one
_triplet=_inp.get('calculate.properties.triplet')
if(_triplet == None):
    _triplet=False
_response_data['triplet']=bool(_triplet)

#Check if we are going to use the B_{\lambda c} term in to our response equations
_use_b=_inp.get('method.MRCC.response.use_b')
if(_use_b == None):
    _use_b=True

#We then change this information to call it as two different 'option'
if(_use_b):
    option=1
else:
    option=2

_response_data['option']=bool(option)

# Get the restart option to skip the calculation of lower order response properties
_restart=_inp.get('calculate.properties.restart')
if(_restart == None):
    _restart=1

_response_data['restart']=int(_restart)

print 'Using following options for the response calculations:'
print _response_data
