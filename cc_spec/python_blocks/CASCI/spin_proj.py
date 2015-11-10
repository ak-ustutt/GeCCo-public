from gecco_interface import *


from gecco_modules.NoticeUtil import *

orbitals=Orb_Info()
keywords=GeCCo_Input()

ims = int(orbitals.get('ims'))       
imult =int(orbitals.get('imult'))

if (ims == 0) and ((imult-1)%4 == 0) :
    msc = 1
elif (ims == 0) and ((imult+1)%4 == 0) :
    msc = -1
else :
    msc = 0                        # msc is the AB_sym for many operators


wf_sym=orbitals.get('lsym')        #symmetry of wavefunction

##################################################################
#Specials needed for spinproj >0
##################################################################
new_target('SpinProjSpec')
depend('RefState-Operators')

CLONE_OPERATOR({
        TEMPLATE:'C0',
        LABEL:'C0_sp'})

DEF_ME_LIST({
        LIST:'C0_sp_LST',
        OPERATOR:'C0_sp',
        IRREP:wf_sym,
        '2MS':ims,
        AB_SYM:msc
        })
#        S2:imult})

debug_MEL('C0_sp_LST',info_only=True)




DEF_OP_FROM_OCC({
        LABEL:'S+',
        JOIN:1,
        DESCR:'H,H|P,P|V,V'})

DEF_ME_LIST({
        LIST:'S+_LST',
        OPERATOR:'S+',
        IRREP:1,
        '2MS':2})

ADD_UNITY({
        LIST:'S+_LST',
        FAC:1.0,
        MS_SYM_SIGN:-1,
        INIT:True})

debug_MEL('S+_LST')




CLONE_OPERATOR({
        TEMPLATE:'S+',
        LABEL:'S-'})

DEF_ME_LIST({
        LIST:'S-_LST',
        OPERATOR:'S-',
        IRREP:1,
        '2MS':-2})

ADD_UNITY({
        LIST:'S-_LST',
        FAC:1.0,
        MS_SYM_SIGN:-1,
        INIT:True})

debug_MEL('S-_LST')




CLONE_OPERATOR({
        TEMPLATE:'S+',
        LABEL:'Sz'})

DEF_ME_LIST({
        LIST:'Sz_LST',
        OPERATOR:'Sz',
        IRREP:1,
        '2MS':0})

ADD_UNITY({
        LIST:'Sz_LST',
        FAC:1.0,
        MS_SYM_SIGN:0,
        INIT:True})

debug_MEL('Sz_LST')




#Formular for Spinprojection C0_sp = S^2 C0
# 1/2*(S+S- + S-S+) = {S+S-} + 1/2*(S+S-)_c + 1/2*(S-S+)_c
EXPAND_OP_PRODUCT({
        LABEL:'FORM_C0_sp',
        OP_RES:'C0_sp',
        OPERATORS: ['C0_sp','S+','S-','C0','C0_sp'],
        IDX_SV:[1,2,3,4,1],
        FAC:1.0,
        AVOID:[2,3],
        NEW:True})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_C0_sp',
        OP_RES:'C0_sp',
        OPERATORS: ['C0_sp','S+','S-','C0','C0_sp'],
        IDX_SV:[1,2,3,4,1],
        FAC:0.5,
        CONNECT:[2,3],
        NEW:False})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_C0_sp',
        OP_RES:'C0_sp',

        OPERATORS: ['C0_sp','S-','S+','C0','C0_sp'],
        IDX_SV:[1,2,3,4,1],
        FAC:0.5,
        CONNECT:[2,3],
        NEW:False})
if  ims != 0 : 
    EXPAND_OP_PRODUCT({
            LABEL:'FORM_C0_sp',
            OP_RES:'C0_sp',
            OPERATORS: ['C0_sp','Sz','Sz','C0','C0_sp'],
            IDX_SV:[1,2,3,4,1],
            FAC:1.0,
            FIX_VTX:True,
            NEW:False})

debug_FORM('FORM_C0_sp')

OPTIMIZE({
        LABEL_OPT:'FOPT_C0_sp',
        LABELS_IN:'FORM_C0_sp'})
