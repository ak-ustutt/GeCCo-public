
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


c0_shape='V'*orbitals.get('nactel')+','



spinadapt=keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0



ims = int(orbitals.get('ims'))       
imult =int(orbitals.get('imult'))

if (ims == 0) and ((imult-1)%4 == 0) :
    msc = 1
elif (ims == 0) and ((imult+1)%4 == 0) :
    msc = -1
else :
    msc = 0                        # msc is the AB_sym for many operators
wf_sym=orbitals.get('lsym')

#------------------------------------------------------------------
# operator definitions for reference state
#------------------------------------------------------------------


new_target('Make_E0')
heading('=== Operators for Reference State ===')

comment('Reference energy')
DEF_SCALAR({
        LABEL:'E0'})

DEF_ME_LIST({
        LIST:'ME_E0',
        OPERATOR:'E0',
        IRREP:1,
        '2MS':0,
        AB_SYM:0
        })

new_target('Make_C0')
comment('Reference coefficients')
DEF_OP_FROM_OCC({
        LABEL:'C0',
        DESCR:c0_shape})
DEF_ME_LIST({
        LIST:'ME_C0',
        OPERATOR:'C0',
        IRREP:wf_sym,
        '2MS':ims,
        AB_SYM:msc
        })
#        S2:imult})

debug_MEL('ME_C0',info_only=True)

new_target('MAKE_D0')
comment('CASCI preconditioner')
CLONE_OPERATOR({
        LABEL:'D0',
        TEMPLATE:'C0'})
DEF_ME_LIST({
        LIST:'ME_D0',
        OPERATOR:'D0',
        IRREP:wf_sym,
        '2MS':ims,
        AB_SYM:msc
        })

depend('H0')



comment('Prepare diagonal ...')

PRECONDITIONER({
        LIST_PRC:'ME_D0',
        LIST_INP:'H0',
        MODE:'dia-H'})

debug_MEL('ME_D0')



new_target('Make_H_C0')
comment('CASCI omg')
CLONE_OPERATOR({
        LABEL:'H_C0',
        TEMPLATE:'C0'})
DEF_ME_LIST({
        LIST:'ME_H_C0',
        OPERATOR:'H_C0',
        IRREP:wf_sym,
        '2MS':ims,
        AB_SYM:msc
        })
#        S2:imult})

debug_MEL('ME_H_C0',info_only=True)


new_target('Make_GAM0')
comment('Density matrix up to third order')
DEF_OP_FROM_OCC({
        LABEL:'GAM0',
        JOIN:2,
        DESCR:',;,|,V;V,|,VV;VV,|,VVV;VVV,'})
DEF_ME_LIST({
        LIST:'GAM0_LST',
        OPERATOR:'GAM0',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

debug_MEL('GAM0_LST',info_only=True)

new_target('RefState-Operators')
depend('Make_E0')
depend('Make_C0')
depend('MAKE_D0')
depend('Make_H_C0')
depend('Make_GAM0')



#----------------------------------------------------------------------
#
new_target("DEF_A_C0")
depend("Make_H_C0")

CLONE_OPERATOR({
     LABEL:"A_C0",
     TEMPLATE:"H_C0"})
DEF_ME_LIST({LIST:"ME_A_C0",
    OPERATOR:"A_C0",
    IRREP:wf_sym, 
    '2MS':ims,
     AB_SYM:msc})

debug_MEL('ME_A_C0',info_only=True)


