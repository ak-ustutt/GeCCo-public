# Get the informations about response calculation from the input
from gecco_interface import *

_inp = GeCCo_Input()

#Get the maximum order of the properties we are going to calculate
_response_data = {}
_order=_inp.get('method.MRCC.response.order')
if(_order == None):
    _order=1
_response_data['order']=int(_order)

#Get the frequency of the perturbation if we are going to calculate dynamic properties
_freq=_inp.get('method.MRCC.response.freq')
if(_freq == None):
    _freq=0.0
_response_data['freq']=float(_freq)

#Check if the perturation is triplet, that makes the resultant property a spin-dependent one
_triplet_str=_inp.get('calculate.properties.triplet')
print 'debugging 1:', _triplet_str
if(_triplet_str=='T'):
    _triplet=True
elif(_triplet_str=='F'):
    _triplet=False
else:
    _triplet=False
    
_response_data['triplet']=bool(_triplet)
print 'debugging 2:', _response_data['triplet']

#Check if we are going to use the B_{\lambda c} term in to our response equations
_use_b_str=_inp.get('method.MRCC.response.use_b')
if(_use_b_str=='T'):
    _use_b=True
elif(_use_b_str=='F'):
    _use_b=False
else:
    _use_b=True
    
#We then change this information to call it as two different 'option'
if(_use_b):
    option=1
else:
    option=2

_response_data['option']=option

# Get the restart option to skip the calculation of lower order response properties
_restart=_inp.get('calculate.properties.restart')
if(_restart == None):
    _restart=1

_response_data['restart']=int(_restart)

#Now we are reading the information about the response calculation written 
#in the file 'prop_info.gecco'

#_pop_data is a dictionary that has the information about the perturbation 
#operators influencing the response calculation

_pop_data={}
#cartisian component of the pertubation operator
_pop_data['comp']=[]
#name of the perturbation operator, for the time bieng only dipole operator ('d') is allowed
_pop_data['name']=[]
#name of the integral corresponding to the perturbation
_pop_data['int_name']=[]
#sign of the perturbation
_pop_data['sign']=[]
#symmetry of the pertubation
_pop_data['isym']=[]

#cmp_data is a dictionary that has the information about all the new components 
#of the wabefunction that has to be solved to get the response properties

_cmp_data={}
#with which perturbation operator this component of the wavefunction is linked to
_cmp_data['pop_idx']=[]
#frequency linked to the component
_cmp_data['freq']=[]
#if calculation of this component is redundant then to which previous component
_cmp_data['redun']=[]

_prop_info_name='prop_info.gecco'

lines=[line.rstrip('\n') for line in open(_prop_info_name)]

# l is the just the line number
l=0

#Now reading line by line from the file

#reading the number of operator influencing the calculation
part=lines[l].split()
_npop=int(part[1])
l=l+1
#print 'npop', _npop
for i in xrange (0,_npop):
    part=lines[l].split()
    _pop_data['comp'].append(str(part[1]))
    part=lines[l+1].split()
    _pop_data['name'].append(str(part[1]))
    if(str(part[1])!='d'):
        quit_error('Properties for ic-MRCC can be calculated only for dipole operators...')
    part=lines[l+2].split()
    _pop_data['int_name'].append(str(part[1]))
    part=lines[l+3].split()
    _pop_data['sign'].append(int(part[1]))
    part=lines[l+4].split()
    _pop_data['isym'].append(int(part[1]))
    l=l+5
#print 'pop_data', _pop_data
part=lines[l].split()
#number of the calculation asked in the input
_ncnt=int(part[0])
#maximum order of the properties that has been asked for
_maxord=int(part[1])
l=l+1
#print '_ncnt', _ncnt, '_maxord', _maxord
for i in xrange (0,_ncnt*_maxord):
    part=lines[l].split()
    _cmp_data['pop_idx'].append(int(part[1]))
    part=lines[l+1].split()
    _cmp_data['freq'].append(float(part[1]))
    part=lines[l+2].split()
    _cmp_data['redun'].append(int(part[1]))
    l=l+3
#print 'cmp_data', _cmp_data

#adding the lone parameters in to the _response_data
_response_data['nPop']=_npop
_response_data['nCnt']=_ncnt
_response_data['maxorder']=_maxord

print 'Using following options for the response calculations:'
print _response_data
