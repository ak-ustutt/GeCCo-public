
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

t2g_shape='V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH' #compatible with Matthias
t1_shape='V,H|P,V|P,H'



new_target('DEF_T2g')
comment('=== Cluster Operators ===')


DEF_OP_FROM_OCC({
        LABEL:'T2g',
        DESCR:t2g_shape})

ME_param={
        LIST:'ME_T2g',
        OPERATOR:'T2g  ',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_T2g',info_only=True)





new_target('DEF_T1 ')

DEF_OP_FROM_OCC({
        LABEL:'T1',
        DESCR:t1_shape})

ME_param={
        LIST:'ME_T1',
        OPERATOR:'T1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_T1',info_only=True)





new_target('DEF_O2g')
comment('=== Residue Operators ===')

DEF_OP_FROM_OCC({
        LABEL:'O2g',
        JOIN:2,
        DESCR:omg_generator(t2g_shape)})

ME_param={
        LIST:'ME_O2g',
        OPERATOR:'O2g',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_O2g',info_only=True)





new_target('DEF_O1')

DEF_OP_FROM_OCC({
        LABEL:'O1',
        JOIN:2,
        DESCR:omg_generator(t1_shape)})

ME_param={
        LIST:'ME_O1',
        OPERATOR:'O1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('ME_O1',info_only=True)




new_target('DEF_LAM2g')
comment('=== Lambda Operators ===')
depend('DEF_T2g')


CLONE_OPERATOR({
        LABEL:'LAM2g',
        TEMPLATE:'T2g  ',
        ADJOINT:True})

ME_param={
        LIST:'ME_LAM2g',
        OPERATOR:'LAM2g',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}

DEF_ME_LIST(ME_param)



debug_MEL('ME_LAM2g',info_only=True)





new_target('DEF_LAM1')
comment('=== Lambda Operators ===')
depend('DEF_T1')

CLONE_OPERATOR({
        LABEL:'LAM1',
        TEMPLATE:'T1',
        ADJOINT:True})

ME_param={
        LIST:'ME_LAM1',
        OPERATOR:'LAM1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}

DEF_ME_LIST(ME_param)

debug_MEL('ME_LAM1',info_only=True)
