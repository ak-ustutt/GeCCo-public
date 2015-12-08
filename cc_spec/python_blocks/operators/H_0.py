from gecco_interface import *
from gecco_modules.NoticeUtil import *


#---------------------------------------------------------------------------------
# Finks excitation retaining hamiltonian
#---------------------------------------------------------------------------------
new_target('MAKE_REPT_HAM')
heading('===Building of excitation retaining hamiltonian===')
#generalized from
#Fink Chemical Physics Letters 428 (2006) 461 DOI: 10.1016/j.cplett.2006.07.081
depend('H0')
depend('MakeRefState')

comment("defining REPT-hamiltonian and building formula")

DEF_OP_FROM_OCC({
        LABEL:'REPT_HAM',
        JOIN:1,
        DESCR:'H,H|V,V|P,P|HH,HH|VV,VV|PP,PP|HV,HV|HP,HP|PV,PV'})

DEF_ME_LIST({
        LIST:'REPT_HAM_LST',
        OPERATOR:'REPT_HAM',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_REPT_HAM',
        NEW:True,
        OP_RES:'REPT_HAM',
        OPERATORS:['REPT_HAM','H','REPT_HAM'],
        IDX_SV:[1,2,1]})

debug_FORM('FORM_REPL_HAM')




new_target('EVAL_REPT_HAM')
depend('MAKE_REPT_HAM')
comment("evaluate REPT-hamiltonian")

OPTIMIZE({
        LABELS_IN:'FORM_REPT_HAM',
        LABEL_OPT:'FOPT_REPT_HAM'})
EVALUATE({
        FORM:'FOPT_REPT_HAM'})


debug_MEL('REPT_HAM_LST')





#---------------------------------------------------------------------------------
# Effective Fock operator
#---------------------------------------------------------------------------------
new_target('MAKE_FOCK_EFF')
heading('===Building of effective fock matrix===')
depend('H0')
depend('MakeRefState')
depend('GAM0_CALC')

comment("defining effective fock operator and building formula")

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
depend('MAKE_FOCK_EFF')
comment("evaluate effective Fock operator")

OPTIMIZE({
        LABELS_IN:'FORM_F_EFF',
        LABEL_OPT:'FOPT_F_EFF'})
EVALUATE({
        FORM:'FOPT_F_EFF'})


debug_MEL('FOCK_EFF_LST')


#because I need it without the v,v part for the preconditioner 
new_target("Make_F_EFF_INACT")
depend('EVAL_F_EFF')

comment("defining effective Fock operator without V,V part and defining Fomula for that")

DEF_OP_FROM_OCC({
        LABEL:'FOCK_EFF_INACT',
        JOIN:1,
        DESCR:'H,H|P,P'})
DEF_ME_LIST({
        LIST:'FOCK_EFF_INACT_LST',
        OPERATOR:'FOCK_EFF_INACT',
        IRREP:1,
        '2MS':0,
        AB_SYM:1})

#F_eff_diag= F_eff limited by blocks of F_eff_diag
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_INACT',
        OP_RES:'FOCK_EFF_INACT',
        OPERATORS:['FOCK_EFF_INACT','FOCK_EFF','FOCK_EFF_INACT'],
        IDX_SV:[1,2,1],
        NEW:True})




new_target('EVAL_F_EFF_INACT')
depend("Make_F_EFF_INACT")

comment("Evaluate effective Fock operator (inactive part only)")

OPTIMIZE({
        LABELS_IN:'FORM_F_EFF_INACT',
        LABEL_OPT:'FOPT_F_EFF_INACT'})
EVALUATE({
        FORM:'FOPT_F_EFF_INACT'})

debug_MEL('FOCK_EFF_INACT_LST')


#---------------------------------------------------------------------------------
#Building Dyalls 0th order Hamiltonian
#---------------------------------------------------------------------------------
new_target('Make_HAM_D')
heading("====== Building of Dyall's 0th order Hamiltonian ======")
depend('EVAL_F_EFF_INACT')
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
        OPERATORS:['HAM_D','FOCK_EFF_INACT','HAM_D'],
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




new_target("EVAL_HAM_D")
#this target should contain the evaluation of FORM_HAM_D but currently there are targets that expect that in 
#Make_HAM_D
depend('Make_HAM_D')






#-------------------------------------------------------------------------------
#the fluctuation hamiltonian for H_0=Ham_D 
#-------------------------------------------------------------------------------
#currently never used
new_target('Make_V')
depend('EVAL_F_EFF')
depend('H0')
depend('Make_HAM_D')

comment("Create perturbing hamiltonian V=H-HAM_D")

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

new_target("EVAL_V")
depend("Make_V")
comment("Evaluate the perturbing operator")
OPTIMIZE({
        LABELS_IN:'FORM_V', 
        LABEL_OPT:'FOPT_V'})
EVALUATE({
        FORM:'FOPT_V'})
