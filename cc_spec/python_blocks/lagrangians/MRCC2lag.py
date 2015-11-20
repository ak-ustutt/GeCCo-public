from gecco_interface import *
from gecco_modules.NoticeUtil import *
import gecco_modules.string_to_form as stf

i_am="MRCC2lag.py"

lag_type=4
if keywords.is_keyword_set('method.MRCC2.lagrangian'):
    lag_type=int(keywords.get('method.MRCC2.lagrangian'))
print("lagrangian:",lag_type,type(lag_type))

known_hamiltonians={"DYALL","REPT","F_EFF"}
hamiltonian="DYALL"
if keywords.is_keyword_set('method.MRCC2.hamiltonian'):
    hamiltonian=str(keywords.get('method.MRCC2.hamiltonian')).strip()
print("hamiltonian: ", hamiltonian, type(hamiltonian))

if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

#--------------------------------------------------------------------------------------#
#Define the lagrangian
#--------------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG2')

depend('T-Operators')

depend('MakeRefState')

depend('H0')

if hamiltonian=="DYALL":
    depend('EVAL_HAM_D')
elif hamiltonian=="REPT":
    depend('EVAL_REPT_HAM')
elif hamiltonian=="F_EFF":
    depend('EVAL_F_EFF')


DEF_SCALAR({
        LABEL:'PT_LAG'})
DEF_ME_LIST({
        LIST:'PT_LAG_LST',
        OPERATOR:'PT_LAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

#linear lagrangian
if lag_type >= 1 :
    #energy expression 
    LAG=stf.Formula("FORM_PT_LAG2:PT_LAG="\
                "<C0^+*("\
                    "H"\
                    #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
                    "+(H*T1_ca)"\
                    "+(H*T2_ca)"\
                ")*C0>"
    )
    #T1 projection
    LAG.append(
        "<C0^+*(T1_ca^+)*("\
            "H"\
            "+[H,T1_ca]"\
            "+[H,T2_ca]"\
        ")*C0>"
    )
    #T2 projection
    LAG.append(
        "<C0^+*(T2_ca^+)*("\
            "H"\
            "+[H,T1_ca]"\
        ")*C0>"             
    )
    if hamiltonian=="DYALL":
        LAG.append("<C0^+*(T2_ca^+)*([HAM_D,T2_ca])*C0>")
    elif hamiltonian=="REPT":
        LAG.append("<C0^+*(T2_ca^+)*([REPT_HAM,T2_ca])*C0>")
    elif hamiltonian=="F_EFF":
        LAG.append("<C0^+*(T2_ca^+)*([F_EFF,T2_ca])*C0>")
#quadratic lagrangian: linear lagrangian+something
#something:
if lag_type >= 2 :
    #energy expression    
    LAG.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            "1/2((H*T1_ca)*T1_ca)"\

            #these terms where not included in CC2-T1q 
            "+1/2((H*T2_ca)*T1_ca)+1/2((H*T1_ca)*T2_ca)"\
        ")*C0>"             
    )
    #T1 projection
    LAG.append(

        "<C0^+*(T1_ca^+)*("\
            "1/2[[H,T1_ca],T1_ca]"\
            "+1/2[[H,T2_ca],T1_ca]+1/2[[H,T1_ca],T2_ca]"\
        ")*C0>"             
    )
    #T2 projection
    LAG.append(
        "<C0^+*(T2_ca^+)*("\
            "1/2[[H,T1_ca],T1_ca]"\
        ")*C0>"             
    )



#kubic lagrangian: quadratic lagrangian+something
#something:

if lag_type >= 3 :
    #energy expression    
    LAG.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            #these terms where not included in CC2-T1c
            "1/6(((H*T1_ca)*T1_ca)*T1_ca)"\
            "+1/6(((H*T2_ca)*T1_ca)*T1_ca)+1/6(((H*T1_ca)*T2_ca)*T1_ca)+1/6(((H*T1_ca)*T1_ca)*T2_ca)"\
        ")*C0>"             
    )
    #T1 projection
    LAG.append(
        "<C0^+*(T1_ca^+)*("\
            "1/6[[[H,T1_ca],T1_ca],T1_ca]"\

            #these terms where not included in CC2-T1c
            "+1/6[[[H,T2_ca],T1_ca],T1_ca]+1/6[[[H,T1_ca],T2_ca],T1_ca]+1/6[[[H,T1_ca],T1_ca],T2_ca]"\
        ")*C0>"             
    )
    #T2 projection
    LAG.append(
        "<C0^+*(T2_ca^+)*("\
            "1/6[[[H,T1_ca],T1_ca],T1_ca]"\
        ")*C0>"             
    )

#quartic lagrangian: cubic lagrangian+something
#something:
if lag_type >= 4 :
    #energy expression    
    LAG.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            #these terms where not included in CC2
            "1/24((((H*T1_ca)*T1_ca)*T1_ca)*T1_ca)"\
            "+1/24((((H*T2_ca)*T1_ca)*T1_ca)*T1_ca)+1/24((((H*T1_ca)*T2_ca)*T1_ca)*T1_ca)+1/24((((H*T1_ca)*T1_ca)*T2_ca)*T1_ca)+1/24((((H*T1_ca)*T1_ca)*T1_ca)*T2_ca)"\
        ")*C0>"             
    )

    #T1 projection
    LAG.append(
        "<C0^+*(T1_ca^+)*("\
            "1/24[[[[H,T1_ca],T1_ca],T1_ca],T1_ca]"\

            #these terms where not included in CC2
            "+1/24[[[[H,T2_ca],T1_ca],T1_ca],T1_ca]+1/24[[[[H,T1_ca],T2_ca],T1_ca],T1_ca]+1/24[[[[H,T1_ca],T1_ca],T2_ca],T1_ca]+1/24[[[[H,T1_ca],T1_ca],T1_ca],T2_ca]"\
        ")*C0>"             
    )
    #T2 projection
    LAG.append(
        "<C0^+*(T2_ca^+)*("\
            "1/24[[[[H,T1_ca],T1_ca],T1_ca],T1_ca]"\
        ")*C0>"             
    )
if not 0<lag_type<5 :
    raise Exception("MRCCPT_split unknown lagrangian type\nlag_type="+str(lag_type))



for item in LAG.show():
    print item
LAG.set_rule()

mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_LAG2'})


debug_FORM('FORM_PT_LAG2')



#Make the Derivative with respect to T2_ca^+ and factor that out. 
#PT_LAG=E+A => PT_LAG=E_1,2,3+T1_ca^+*O1+T2_ca^+*Oges
DERIVATIVE({LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_Ampl_T1',
        OP_RES:'O1',
        OP_DERIV:'T1_ca^+'})

DERIVATIVE({LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_Ampl_T2',
        OP_RES:'Oges',
        OP_DERIV:'T2_ca^+'})


debug_FORM('FORM_PT_Ampl_T2')

debug_FORM('FORM_PT_Ampl_T1')

# For some reason the reorder formulas are necessary. 
REORDER_FORMULA({
        LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_LAG2'})

REORDER_FORMULA({
        LABEL_IN:'FORM_PT_Ampl_T1',
        LABEL_RES:'FORM_PT_Ampl_T1'})
FACTOR_OUT({
        LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_LAG2',
        INTERM:'FORM_PT_Ampl_T1'})

REORDER_FORMULA({
        LABEL_IN:'FORM_PT_Ampl_T2',
        LABEL_RES:'FORM_PT_Ampl_T2'})
FACTOR_OUT({
        LABEL_IN:'FORM_PT_LAG2',
        LABEL_RES:'FORM_PT_LAG2',
        INTERM:'FORM_PT_Ampl_T2'})


debug_FORM('FORM_PT_LAG2')

debug_FORM('FORM_PT_Ampl_T2')

debug_FORM('FORM_PT_Ampl_T1')

OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG2',
        LABELS_IN:['FORM_PT_Ampl_T2','FORM_PT_Ampl_T1','FORM_PT_LAG2']})

