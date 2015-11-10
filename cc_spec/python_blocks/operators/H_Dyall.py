from gecco_interface import *
from python_modules.NoticeUtil import *


#---------------------------------------------------------------------------------
# Effective Fock operator
#---------------------------------------------------------------------------------

new_target('MAKE_F_EFF')
PRINT({STRING:'Building of effective fock matrix'})
depend('H0')
depend('MakeRefState')
depend('GAM0_CALC')





DEF_HAMILTONIAN({
        LABEL:'FOCK_EFF',
        MIN_RANK:1,
        MAX_RANK:1})

DEF_ME_LIST({
        LIST:'FOCK_EFF_LST',
        OPERATOR:'FOCK_EFF',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})

#F_eff= F + g^pv_qu*gam^u_q  implicitly limited since F_eff has only rank 1
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF',
        NEW:True,
        OP_RES:'FOCK_EFF',
        OPERATORS:['FOCK_EFF','H','FOCK_EFF'],
        IDX_SV:[1,2,1]})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF',
        NEW:False,
        OP_RES:'FOCK_EFF',
        OPERATORS:['FOCK_EFF','GAM0','H','GAM0','FOCK_EFF'],
        LABEL_DESCR:['2,,,V','2,3,,V','3,4,,V'],
        CONNECT:[1,3,3,5],
        IDX_SV:[1,2,3,2,1]
        })

debug_FORM('FORM_F_EFF')


new_target('EVAL_F_EFF')
depend('MAKE_F_EFF')

OPTIMIZE({
        LABELS_IN:'FORM_F_EFF',
        LABEL_OPT:'FOPT_F_EFF'})
EVALUATE({
        FORM:'FOPT_F_EFF'})


debug_MEL('FOCK_EFF_LST')


#because I need it without the v,v part for the preconditioner 
new_target('EVAL_F_EFF_DIAG')
depend('EVAL_F_EFF')

DEF_OP_FROM_OCC({
        LABEL:'FOCK_EFF_DIAG',
        JOIN:1,
        DESCR:'H,H|P,P'})
DEF_ME_LIST({
        LIST:'FOCK_EFF_DIAG_LST',
        OPERATOR:'FOCK_EFF_DIAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:1})

#F_eff_diag= F_eff limited by blocks of F_eff_diag
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_DIAG',
        OP_RES:'FOCK_EFF_DIAG',
        OPERATORS:['FOCK_EFF_DIAG','FOCK_EFF','FOCK_EFF_DIAG'],
        IDX_SV:[1,2,1],
        NEW:True})

OPTIMIZE({
        LABELS_IN:'FORM_F_EFF_DIAG',
        LABEL_OPT:'FOPT_F_EFF_DIAG'})
EVALUATE({
        FORM:'FOPT_F_EFF_DIAG'})

debug_MEL('FOCK_EFF_DIAG_LST')


#---------------------------------------------------------------------------------
#Building Dyalls 0th order Hamiltonian
#---------------------------------------------------------------------------------
new_target('DyallsHamiltonian')
heading("====== Building of Dyall's 0th order Hamiltonian ======")
depend('EVAL_F_EFF_DIAG')

depend('H0')





DEF_OP_FROM_OCC({
        LABEL:'HAM_D',
        DESCR:'H,H|P,P|V,V|VV,VV'})

DEF_ME_LIST({
        LIST:'HAM_D_LIST',
        OPERATOR:'HAM_D',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

debug_MEL('HAM_D_LIST',info_only=True)


 


EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D',
        OP_RES:'HAM_D',
        OPERATORS:['HAM_D','FOCK_EFF_DIAG','HAM_D'],
        IDX_SV:[1,2,1],
        NEW:True})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D',
        OP_RES:'HAM_D',
        OPERATORS:['HAM_D','H','HAM_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,V,V',
        NEW:False})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D',
        OP_RES:'HAM_D',
        OPERATORS:['HAM_D','H','HAM_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,VV,VV',
        NEW:False})

debug_FORM('FORM_HAM_D')




OPTIMIZE({
        LABELS_IN:'FORM_HAM_D',
        LABEL_OPT:'FOPT_HAM_D'})

EVALUATE({
        FORM:'FOPT_HAM_D'})

debug_MEL('HAM_D_LIST')




new_target('HAMILTONIAN_HAM_D')
depend('EVAL_F_EFF')
depend('H0')
depend('DyallsHamiltonian')


DEF_HAMILTONIAN({
        LABEL:'V',
        MIN_RANK:0, 
        MAX_RANK:2})
DEF_ME_LIST({LIST:'V_LST',
        OPERATOR:'V',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_V',
        OP_RES:'V',
        OPERATORS:['V','H','V'],
        IDX_SV:[1,2,1],
        NEW:True})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_V',
        OP_RES:'V',
        OPERATORS:['V','HAM_D','V'],
        IDX_SV:[1,2,1],
        FAC:-1.0,
        NEW:False})

OPTIMIZE({
        LABELS_IN:'FORM_V', 
        LABEL_OPT:'FOPT_V'})
EVALUATE({
        FORM:'FOPT_V'})
