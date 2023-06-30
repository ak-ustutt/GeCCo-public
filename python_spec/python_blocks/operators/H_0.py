from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

i_am='H_0'


#---------------------------------------------------------------------------------
# Finks excitation retaining hamiltonian
#---------------------------------------------------------------------------------
new_target('MAKE_REPT_HAM')
heading('Building of excitation retaining hamiltonian')
#generalized from
#Fink Chemical Physics Letters 428 (2006) 461 DOI: 10.1016/j.cplett.2006.07.081
depend('H0')
depend('MakeRefState')

comment("defining REPT-hamiltonian and building formula")

DEF_OP_FROM_OCC({
        LABEL:'REPT_HAM',
        JOIN:1,
        DESCR:',|H,H|V,V|P,P|HH,HH|VV,VV|PP,PP|HV,HV|HP,HP|PV,PV'})

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
# Simple Fock operator 
#---------------------------------------------------------------------------------
new_target('MAKE_FOCK_REF')

heading('Building of simple fock matrix')
depend('H0')
depend('MakeRefState')
depend('GAM00_CALC')

DEF_OP_FROM_OCC({
        LABEL:'FOCK_REF',
        JOIN:1,
        DESCR:',|H,H|P,P|H,P|P,H'})

DEF_ME_LIST({
        LIST:'ME_FOCK_REF',
        OPERATOR:'FOCK_REF',
        IRREP:1,
        '2MS':0,
	    AB_SYM:1})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_FOCK_REF',
        NEW:True,
        OP_RES:'FOCK_REF',
        OPERATORS:['FOCK_REF','GAM00','H','GAM00','FOCK_REF'],
        AVOID:[1,5,2,4],
        IDX_SV:[1,2,3,2,1]
        })
OPTIMIZE({
        LABEL_OPT:"FOPT_FOCK_REF",
        LABELS_IN:"FORM_FOCK_REF"})

EVALUATE({
        FORM:"FOPT_FOCK_REF"})

#---------------------------------------------------------------------------------
# Effective Fock operator, diagonal part only
#---------------------------------------------------------------------------------
new_target('MAKE_FOCK_EFF_D')
heading('Building of effective fock matrix (diagonal blocks)')
depend('H0')
depend('MakeRefState')
depend('GAM00_CALC')

comment("defining effective fock operator and building formula")

DEF_OP_FROM_OCC({
        LABEL:'FOCK_EFF_D',
        DESCR:'H,H|V,V|P,P'
        })

DEF_SCALAR({
        LABEL:'FOCK_EFF_D_EXP'})

DEF_ME_LIST({
        LIST:'FOCK_EFF_D_LST',
        OPERATOR:'FOCK_EFF_D',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})

DEF_ME_LIST({
        LIST:'FOCK_EFF_D_EXP_LST',
        OPERATOR:'FOCK_EFF_D_EXP',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})


#F_eff= F + g^pv_qu*gam^u_q  implicitly limited since F_eff has only rank 1
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_D',
        NEW:True,
        OP_RES:'FOCK_EFF_D',
        OPERATORS:['FOCK_EFF_D','H','FOCK_EFF_D'],
        IDX_SV:[1,2,1]})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_D',
        NEW:False,
        OP_RES:'FOCK_EFF_D',
        OPERATORS:['FOCK_EFF_D','GAM00','H','GAM00','FOCK_EFF_D'],
        LABEL_DESCR:['2,,,V'],
        CONNECT:[2,3],
        AVOID:[1,2,1,4,2,5,4,5],
        IDX_SV:[1,2,3,2,1]
        })

# expectation value <F_EFF_D> - the active energy is sufficient
# note that GAM0 is used here, not GAM00 ... the operator is 
# defined for GAM00 but the exp. value for the actual density

EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_D_EXP',
        NEW:True,
        FAC:1.0,
        OP_RES:'FOCK_EFF_D_EXP',
        OPERATORS:['GAM0','FOCK_EFF_D','GAM0'],
        LABEL_DESCR:['1,,,V'],
        CONNECT:[1,2],
        IDX_SV:[1,2,1]
        })

debug_FORM('FORM_F_EFF_D')
debug_FORM('FORM_F_EFF_D_EXP')


new_target('EVAL_F_EFF_D')
depend('MAKE_FOCK_EFF_D')
comment("evaluate effective Fock operator")

OPTIMIZE({
        LABELS_IN:['FORM_F_EFF_D','FORM_F_EFF_D_EXP'],
        LABEL_OPT:'FOPT_F_EFF_D'})
EVALUATE({
        FORM:'FOPT_F_EFF_D'})


debug_MEL('FOCK_EFF_D_LST')
debug_MEL('FOCK_EFF_D_EXP_LST')

#---------------------------------------------------------------------------------
# Effective Fock operator
#---------------------------------------------------------------------------------
new_target('MAKE_FOCK_EFF')
heading('Building of effective fock matrix')
depend('H0')
depend('MakeRefState')
depend('GAM00_CALC')

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

DEF_SCALAR({
        LABEL:'FOCK_EFF_EXP'})

DEF_ME_LIST({
        LIST:'FOCK_EFF_EXP_LST',
        OPERATOR:'FOCK_EFF_EXP',
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
        OPERATORS:['FOCK_EFF','GAM00','H','GAM00','FOCK_EFF'],
        LABEL_DESCR:['2,,,V'],
        CONNECT:[2,3],
        AVOID:[1,2,1,4,2,5,4,5],
        IDX_SV:[1,2,3,2,1]
        })

# expectation value <F_EFF> - the active energy is sufficient
# regarding GAM00 and GAM0 see comments above

EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_EXP',
        NEW:True,
        FAC:1.0,
        OP_RES:'FOCK_EFF_EXP',
        OPERATORS:['GAM0','FOCK_EFF','GAM0'],
        LABEL_DESCR:['1,,,V'],
        CONNECT:[1,2],
        IDX_SV:[1,2,1]
        })

debug_FORM('FORM_F_EFF')
debug_FORM('FORM_F_EFF_EXP')


new_target('EVAL_F_EFF')
depend('MAKE_FOCK_EFF')
comment("evaluate effective Fock operator")

OPTIMIZE({
        LABELS_IN:['FORM_F_EFF','FORM_F_EFF_EXP'],
        LABEL_OPT:'FOPT_F_EFF'})
EVALUATE({
        FORM:'FOPT_F_EFF'})


debug_MEL('FOCK_EFF_LST')
debug_MEL('FOCK_EFF_EXP_LST')


#because for DYALL I need it without the v,v part for the preconditioner 
new_target("Make_F_EFF_INACT")
depend('EVAL_F_EFF')

comment("defining effective Fock operator without V,V part and defining Fomula for that")

DEF_OP_FROM_OCC({
        LABEL:'FOCK_EFF_INACT',
        JOIN:1,
        DESCR:',|H,H|P,P'})
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
heading("Building of Dyall's 0th order Hamiltonian")
depend('EVAL_F_EFF_INACT')
depend('H0')


DEF_OP_FROM_OCC({
        LABEL:'HAM_D',
        DESCR:',|H,H|P,P|V,V|VV,VV'})

DEF_ME_LIST({
        LIST:'HAM_D_LIST',
        OPERATOR:'HAM_D',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

debug_MEL('HAM_D_LIST')

DEF_SCALAR({
        LABEL:'HAM_D_EXP'})

DEF_ME_LIST({
        LIST:'HAM_D_EXP_LIST',
        OPERATOR:'HAM_D_EXP',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})



 
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D',
        OP_RES:'HAM_D',
        OPERATORS:['HAM_D','FOCK_EFF_INACT','HAM_D'],
        IDX_SV:[1,2,1],
        NEW:True})

# The scalar part of HAM_D
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D',
        OP_RES:'HAM_D',
        OPERATORS:['HAM_D','H','HAM_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,,',
        NEW:False})


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

# expectation value <HAM_D> - the active energy is sufficient
# regarding GAM00 and GAM0 see comments above

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_D_EXP',
        NEW:True,
        FAC:1.0,
        OP_RES:'HAM_D_EXP',
        OPERATORS:['GAM0','HAM_D','GAM0'],
        CONNECT:[1,2],
        AVOID:[1,3],
        IDX_SV:[1,2,1]
        })

debug_FORM('FORM_HAM_D_EXP')


OPTIMIZE({
        LABELS_IN:['FORM_HAM_D','FORM_HAM_D_EXP'],
        LABEL_OPT:'FOPT_HAM_D'})

EVALUATE({
        FORM:'FOPT_HAM_D'})

debug_MEL('HAM_D_LIST')
debug_MEL('HAM_D_EXP_LIST')




new_target("EVAL_HAM_D")
#this target should contain the evaluation of FORM_HAM_D but currently there are targets that expect that in 
#Make_HAM_D
depend('Make_HAM_D')



# for the extended DYALL I need some more parts
new_target("Make_F_EFF_4DX")
depend('EVAL_F_EFF')

comment("defining effective Fock operator without V,V part and defining Fomula for that")

DEF_OP_FROM_OCC({
        LABEL:'FOCK_EFF_4DX',
        JOIN:1,
        DESCR:'H,H|P,P|H,V|V,H|H,P|P,H|V,P|P,V'})
DEF_ME_LIST({
        LIST:'FOCK_EFF_4DX_LST',
        OPERATOR:'FOCK_EFF_4DX',
        IRREP:1,
        '2MS':0,
        AB_SYM:1})

#F_eff_diag= F_eff limited by blocks of F_eff_diag
EXPAND_OP_PRODUCT({
        LABEL:'FORM_F_EFF_4DX',
        OP_RES:'FOCK_EFF_4DX',
        OPERATORS:['FOCK_EFF_4DX','FOCK_EFF','FOCK_EFF_4DX'],
        IDX_SV:[1,2,1],
        NEW:True})




new_target('EVAL_F_EFF_4DX')
depend("Make_F_EFF_4DX")

comment("Evaluate effective Fock operator (inactive part only)")

OPTIMIZE({
        LABELS_IN:'FORM_F_EFF_4DX',
        LABEL_OPT:'FOPT_F_EFF_4DX'})
EVALUATE({
        FORM:'FOPT_F_EFF_4DX'})

debug_MEL('FOCK_EFF_4DX_LIST')



#---------------------------------------------------------------------------------
# Dyalls 0th order Hamiltonian (extended version with off-diagonal Fock)
#---------------------------------------------------------------------------------
new_target('EVAL_HAM_DX')
heading("Building of Dyall's 0th order Hamiltonian (extended)")
depend('EVAL_F_EFF_4DX')
depend('H0')


DEF_OP_FROM_OCC({
        LABEL:'HAM_DX',
        DESCR:'H,H|P,P|H,V|V,H|H,P|P,H|V,P|P,V|V,V|VV,VV'})

DEF_ME_LIST({
        LIST:'HAM_DX_LIST',
        OPERATOR:'HAM_DX',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

debug_MEL('HAM_DX_LIST')


 
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_DX',
        OP_RES:'HAM_DX',
        OPERATORS:['HAM_DX','FOCK_EFF_4DX','HAM_DX'],
        IDX_SV:[1,2,1],
        NEW:True})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_DX',
        OP_RES:'HAM_DX',
        OPERATORS:['HAM_DX','H','HAM_DX'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,V,V',
        NEW:False})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_DX',
        OP_RES:'HAM_DX',
        OPERATORS:['HAM_DX','H','HAM_DX'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,VV,VV',
        NEW:False})

debug_FORM('FORM_HAM_DX')


OPTIMIZE({
        LABELS_IN:'FORM_HAM_DX',
        LABEL_OPT:'FOPT_HAM_DX'})

EVALUATE({
        FORM:'FOPT_HAM_DX'})

debug_MEL('HAM_DX_LIST')





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


# -------------------
# Some Hamiltonians in between Dyall's and Fink's
#
#
# Select below which couplings will be added to the Hamiltonian
ham_classes=[
    'HP,HP',
    'HH,HH',
    'PP,PP'
    ]

extra_classes = ''
if ham_classes:
    extra_classes = '|' + '|'.join(ham_classes)

#---------------------------------------------------------------------------------
# Extended Dyalls 0th order Hamiltonian
#
# Does not have couplings between active and inactive parts
# so it still uses f_eff for the one-electron inactive parts
#
#---------------------------------------------------------------------------------

new_target('Make_HAM_EXT_D')
heading("Building of Extended Dyall's 0th order Hamiltonian")
depend('EVAL_F_EFF_INACT')
depend('H0')

DEF_OP_FROM_OCC({
        LABEL:'HAM_EXT_D',
        DESCR:',|H,H|P,P|V,V|VV,VV' + extra_classes})

DEF_ME_LIST({
        LIST:'HAM_EXT_D_LIST',
        OPERATOR:'HAM_EXT_D',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

debug_MEL('HAM_EXT_D_LIST')

DEF_SCALAR({
        LABEL:'HAM_EXT_D_EXP'})

DEF_ME_LIST({
        LIST:'HAM_EXT_D_EXP_LIST',
        OPERATOR:'HAM_EXT_D_EXP',
        IRREP:1,
        '2MS':0,
	AB_SYM:1})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_EXT_D',
        OP_RES:'HAM_EXT_D',
        OPERATORS:['HAM_EXT_D','FOCK_EFF_INACT','HAM_EXT_D'],
        IDX_SV:[1,2,1],
        NEW:True})

# The scalar part of HAM_EXT_D
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_EXT_D',
        OP_RES:'HAM_EXT_D',
        OPERATORS:['HAM_EXT_D','H','HAM_EXT_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,,',
        NEW:False})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_EXT_D',
        OP_RES:'HAM_EXT_D',
        OPERATORS:['HAM_EXT_D','H','HAM_EXT_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,V,V',
        NEW:False})
EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_EXT_D',
        OP_RES:'HAM_EXT_D',
        OPERATORS:['HAM_EXT_D','H','HAM_EXT_D'],
        IDX_SV:[1,2,1],
        LABEL_DESCR:'2,,VV,VV',
        NEW:False})
for hamC in ham_classes:
    EXPAND_OP_PRODUCT({
            LABEL:'FORM_HAM_EXT_D',
            OP_RES:'HAM_EXT_D',
            OPERATORS:['HAM_EXT_D','H','HAM_EXT_D'],
            IDX_SV:[1,2,1],
            LABEL_DESCR:'2,,' + hamC,
            NEW:False})

debug_FORM('FORM_HAM_EXT_D')

# expectation value <HAM_D> - the active energy is sufficient
# regarding GAM00 and GAM0 see comments above

EXPAND_OP_PRODUCT({
        LABEL:'FORM_HAM_EXT_D_EXP',
        NEW:True,
        FAC:1.0,
        OP_RES:'HAM_EXT_D_EXP',
        OPERATORS:['GAM0','HAM_EXT_D','GAM0'],
        CONNECT:[1,2],
        AVOID:[1,3],
        IDX_SV:[1,2,1]
        })

debug_FORM('FORM_HAM_EXT_D_EXP')

OPTIMIZE({
        LABELS_IN:['FORM_HAM_EXT_D','FORM_HAM_EXT_D_EXP'],
        LABEL_OPT:'FOPT_HAM_EXT_D'})

EVALUATE({
        FORM:'FOPT_HAM_EXT_D'})

debug_MEL('HAM_EXT_D_LIST')
debug_MEL('HAM_EXT_D_EXP_LIST')

new_target("EVAL_HAM_EXT_D")
#this target should contain the evaluation of FORM_HAM_EXT_D but currently there are targets that expect that in
#Make_HAM_EXT_D
depend('Make_HAM_EXT_D')


#---------------------------------------------------------------------------------
# Simplified Finks excitation retaining hamiltonian
#
# It has the couplings between active and inactive parts
# Thus it uses the Hamiltonian element for the one-electrons
# inactive parts
#
#---------------------------------------------------------------------------------
new_target('MAKE_SIMP_REPT_HAM')
heading('Building of simplified excitation retaining hamiltonian')
depend('H0')
depend('MakeRefState')

comment("defining simpl REPT-hamiltonian and building formula")

DEF_OP_FROM_OCC({
        LABEL:'SIMP_REPT_HAM',
        JOIN:1,
        DESCR:',|H,H|V,V|P,P|VV,VV|HV,HV|PV,PV' + extra_classes})

DEF_ME_LIST({
        LIST:'SIMP_REPT_HAM_LST',
        OPERATOR:'SIMP_REPT_HAM',
        IRREP:1,
        '2MS':0,
        AB_SYM:1})

EXPAND_OP_PRODUCT({
        LABEL:'FORM_SIMP_REPT_HAM',
        NEW:True,
        OP_RES:'SIMP_REPT_HAM',
        OPERATORS:['SIMP_REPT_HAM','H','SIMP_REPT_HAM'],
        IDX_SV:[1,2,1]})

debug_FORM('FORM_SIMP_REPL_HAM')

new_target('EVAL_SIMP_REPT_HAM')
depend('MAKE_SIMP_REPT_HAM')
comment("evaluate SIMP_REPT-hamiltonian")

OPTIMIZE({
        LABELS_IN:'FORM_SIMP_REPT_HAM',
        LABEL_OPT:'FOPT_SIMP_REPT_HAM'})
EVALUATE({
        FORM:'FOPT_SIMP_REPT_HAM'})

debug_MEL('SIMP_REPT_HAM_LST')
