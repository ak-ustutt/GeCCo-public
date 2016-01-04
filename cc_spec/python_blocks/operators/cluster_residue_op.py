
from gecco_interface import *
from gecco_modules.NoticeUtil import *

from gecco_modules.omg_generator import *

orbitals=Orb_Info()
keywords=GeCCo_Input()

spinadapt=0
if keywords.is_keyword_set('calculate.routes.spinadapt'):
    spinadapt=int(keywords.get('calculate.routes.spinadapt'))

#-----------------------------------------------------------------#
# operators associated with T
#-----------------------------------------------------------------#

t_shape='V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH' #compatible with Matthias
t1_shape='V,H|P,V|P,H'



new_target('T-Operators')
comment('=== Cluster Operators ===')


DEF_OP_FROM_OCC({
        LABEL:'T2_ca',
        DESCR:t_shape})

ME_param={
        LIST:'T2_ca_LST',
        OPERATOR:'T2_ca',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('T2_ca_LST',info_only=True)


DEF_OP_FROM_OCC({
        LABEL:'T1_ca',
        DESCR:t1_shape})

ME_param={
        LIST:'T1_ca_LST',
        OPERATOR:'T1_ca',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('T1_ca_LST',info_only=True)





#new_target('DEF_Oges')
comment('=== Residue Operators ===')

DEF_OP_FROM_OCC({
        LABEL:'Oges',
        JOIN:2,
        DESCR:omg_generator(t_shape)})

ME_param={
        LIST:'Oges_LST',
        OPERATOR:'Oges',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('Oges_LST',info_only=True)




DEF_OP_FROM_OCC({
        LABEL:'O1',
        JOIN:2,
        DESCR:omg_generator(t1_shape)})

ME_param={
        LIST:'O1_LST',
        OPERATOR:'O1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('O1_LST',info_only=True)




new_target('DEF_LAM')
comment('=== Lambda Operators ===')
depend('T-Operators')


CLONE_OPERATOR({
        LABEL:'LAMges',
        TEMPLATE:'T2_ca',
        ADJOINT:True})

ME_param={
        LIST:'LAMges_LST',
        OPERATOR:'LAMges',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}

DEF_ME_LIST(ME_param)

debug_MEL('LAMges_LST',info_only=True)



CLONE_OPERATOR({
        LABEL:'LAM1',
        TEMPLATE:'T1_ca',
        ADJOINT:True})

ME_param={
        LIST:'LAM1_LST',
        OPERATOR:'LAM1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}

DEF_ME_LIST(ME_param)

debug_MEL('LAM1_LST',info_only=True)
