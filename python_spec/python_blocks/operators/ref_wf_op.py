
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


ciroot=1
if keywords.is_keyword_set('method.MR_P.ciroot'):
    ciroot=int(keywords.get('method.MR_P.ciroot'))
print("ciroot: "+str(ciroot))

maxroot=ciroot
if keywords.is_keyword_set('method.MR_P.maxroot'):
    maxroot=int(keywords.get('method.MR_P.maxroot'))
print("maxroot: "+str(maxroot))


c0_shape='V'*orbitals.get('nactel')+','



spinadapt=keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0


triples = keywords.get('method.MR.triples')
triples = triples if triples is not None else 0


ims = int(orbitals.get('ims'))       
imult =int(orbitals.get('imult'))

s_ref=-1  # for reference state and related quantities
s_op=-1   # for operators
msc_op=0
if (spinadapt > 0):
    s_ref=imult  # currently, the program uses projection operators
                 # the solver takes care of that and reads imult from input
                 # to be fixed as we might once allow for different references
                 # with different S in the same run ...
    s_op=0        # spin-free operators
    msc_op=1      # are singlet operators

if (ims == 0) and ((imult-1)%4 == 0) :
    msc = 1
elif (ims == 0) and ((imult+1)%4 == 0) :
    msc = -1
else :
    msc = 0                        # msc is the AB_sym for operators directly related to reference symmetry
wf_sym=orbitals.get('lsym')

#------------------------------------------------------------------
# operator definitions for reference state
#------------------------------------------------------------------


new_target('Define_E0')
heading('Operators for Reference State')

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

new_target('Define_C0')
comment('Reference coefficients')
DEF_OP_FROM_OCC({
        LABEL:'C0',
        DESCR:c0_shape})
DEF_ME_LIST({
        LIST:'ME_C0',
        OPERATOR:'C0',
        IRREP:wf_sym,
        '2MS':ims,
        # S2:s_ref,  # currently not supported, done by projection in SOLVE_EVP
        AB_SYM:msc,
        MIN_REC:1,
        MAX_REC:maxroot
        })
# C00 keeps the unrelaxed reference
DEF_OP_FROM_OCC({
        LABEL:'C00',
        DESCR:c0_shape})
DEF_ME_LIST({
        LIST:'ME_C00',
        OPERATOR:'C00',
        IRREP:wf_sym,
        '2MS':ims,
        # S2:s_ref,  # currently not supported, done by projection in SOLVE_EVP
        AB_SYM:msc,
        MIN_REC:1,
        MAX_REC:maxroot
        })

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



new_target('Define_H_C0')
comment('CASCI omg')
CLONE_OPERATOR({
        LABEL:'H_C0',
        TEMPLATE:'C0'})
DEF_ME_LIST({
        LIST:'ME_H_C0',
        OPERATOR:'H_C0',
        IRREP:wf_sym,
        '2MS':ims,
        #S2:s_ref,  # currently not supported, see above
        AB_SYM:msc
        })

debug_MEL('ME_H_C0',info_only=True)


# note: we should probably have a better criterion:
new_target('Define_GAM0')
if (spinadapt > 0) :
    if (triples != '5' and triples != 'F'):
        comment('Density matrix up to fourth order')
        dstr=',;,|,V;V,|,VV;VV,|,VVV;VVV,|,VVVV;VVVV,'
    else:
        comment('Density matrix up to fifth order')
        dstr=',;,|,V;V,|,VV;VV,|,VVV;VVV,|,VVVV;VVVV,|,VVVVV;VVVVV,'
else:
    comment('Density matrix up to third order')
    dstr=',;,|,V;V,|,VV;VV,|,VVV;VVV,'
DEF_OP_FROM_OCC({
        LABEL:'GAM0',
        JOIN:2,
        DESCR:dstr})
DEF_ME_LIST({
        LIST:'GAM0_LST',
        OPERATOR:'GAM0',
        IRREP:1,
        '2MS':0,S2:s_op,
        AB_SYM:msc_op})
# for the reference density we only need one and two-body part:
dstr=',;,|,V;V,|,VV;VV,'
DEF_OP_FROM_OCC({
        LABEL:'GAM00',
        JOIN:2,
        DESCR:dstr})
DEF_ME_LIST({
        LIST:'GAM00_LST',
        OPERATOR:'GAM00',
        IRREP:1,
        '2MS':0,S2:s_op,
        AB_SYM:msc_op})

debug_MEL('GAM0_LST',info_only=True)

new_target('RefState-Operators')
depend('Define_E0')
depend('Define_C0')
depend('MAKE_D0')
depend('Define_H_C0')
depend('Define_GAM0')



#----------------------------------------------------------------------
#
new_target("DEF_A_C0")
depend("Define_H_C0")

CLONE_OPERATOR({
     LABEL:"A_C0",
     TEMPLATE:"H_C0"})
DEF_ME_LIST({LIST:"ME_A_C0",
    OPERATOR:"A_C0",
    IRREP:wf_sym, 
    '2MS':ims,
     #S2:s_ref,    # currently not supported, see above
     AB_SYM:msc})

debug_MEL('ME_A_C0',info_only=True)


