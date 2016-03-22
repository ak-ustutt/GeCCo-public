from gecco_interface import *
from gecco_modules.NoticeUtil import *
import gecco_modules.string_to_form as stf

i_am="MRCC2lag.py"

lag_type=4
if keywords.is_keyword_set('method.MRCC2.lagrangian'):
    lag_type=int(keywords.get('method.MRCC2.lagrangian'))
print("lagrangian:",lag_type,type(lag_type))

known_hamiltonians=["DYALL","REPT","F_EFF"]
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

depend('DEF_T2g')
depend('DEF_T1')

depend('DEF_LAM2g')
depend('DEF_LAM1')

depend('DEF_O2g')
depend('DEF_O1')


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

DEF_SCALAR({
        LABEL:'PT_LAG_A1'})

DEF_SCALAR({
        LABEL:'PT_LAG_A2'})



#linear lagrangian
if lag_type >= 1 :
    #energy expression 
    LAG_E=stf.Formula("FORM_PT_LAG:PT_LAG="\
                "<C0^+*("\
                    "H"\
                    #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
                    "+(H*T1)"\
                    "+(H*T2g)"\
                ")*C0>"
    )

    #T1 amplitudes
    LAG_A1=stf.Formula("FORM_PT_LAG_A1:PT_LAG_A1="\
        "<C0^+*(LAM1)*("\
            "H"\
            "+[H,T1]"\
            "+[H,T2g]"\
        ")*C0>"
    )
    #T2 amplitudes
    LAG_A2=stf.Formula("FORM_PT_LAG_A2:PT_LAG_A2="\
        "<C0^+*(LAM2g)*("\
            "H"\
            "+[H,T1]"\
        ")*C0>"             
    )
    if hamiltonian=="DYALL":
        LAG_A2.append("<C0^+*(LAM2g)*([HAM_D,T2g])*C0>")
    elif hamiltonian=="REPT":
        LAG_A2.append("<C0^+*(LAM2g)*([REPT_HAM,T2g])*C0>")
    elif hamiltonian=="F_EFF":
        LAG_A2.append("<C0^+*(LAM2g)*([FOCK_EFF,T2g])*C0>")
#quadratic lagrangian: linear lagrangian+something
#something:
if lag_type >= 2 :
    #energy expression    
    LAG_E.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            "1/2((H*T1)*T1)"\

            #these terms where not included in SRCC2-T1q 
            "+1/2((H*T2g)*T1)+1/2((H*T1)*T2g)"\
        ")*C0>"             
    )
    #T1 amplitudes
    LAG_A1.append(

        "<C0^+*(LAM1)*("\
            "1/2[[H,T1],T1]"\
            "+1/2[[H,T2g],T1]+1/2[[H,T1],T2g]"\
        ")*C0>"             
    )
    #T2 amplitudes
    LAG_A2.append(
        "<C0^+*(LAM2g)*("\
            "1/2[[H,T1],T1]"\
        ")*C0>"             
    )



#kubic lagrangian: quadratic lagrangian+something
#something:

if lag_type >= 3 :
    #energy expression    
    LAG_E.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            #these terms where not included in SRCC2-T1c
            "1/6(((H*T1)*T1)*T1)"\
            "+1/6(((H*T2g)*T1)*T1)+1/6(((H*T1)*T2g)*T1)+1/6(((H*T1)*T1)*T2g)"\
        ")*C0>"             
    )
    #T1 amplitudes
    LAG_A1.append(
        "<C0^+*(LAM1)*("\
            "1/6[[[H,T1],T1],T1]"\

            #these terms where not included in SRCC2-T1c
            "+1/6[[[H,T2g],T1],T1]+1/6[[[H,T1],T2g],T1]+1/6[[[H,T1],T1],T2g]"\
        ")*C0>"             
    )
    #T2 amplitudes
    LAG_A2.append(
        "<C0^+*(LAM2g)*("\
            "1/6[[[H,T1],T1],T1]"\
        ")*C0>"             
    )

#quartic lagrangian: cubic lagrangian+something
#something:
if lag_type >= 4 :
    #energy expression    
    LAG_E.append(
        "<C0^+*("\
            #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
            #these terms where not included in SRCC2
            "1/24((((H*T1)*T1)*T1)*T1)"\
            "+ 1/24((((H*T2g)*T1)*T1)*T1) + 1/24((((H*T1)*T2g)*T1)*T1) + 1/24((((H*T1)*T1)*T2g)*T1) + 1/24((((H*T1)*T1)*T1)*T2g)"\
        ")*C0>"             
    )

    #T1 projection
    LAG_A1.append(
        "<C0^+*(LAM1)*("\
            "1/24[[[[H,T1],T1],T1],T1]"\

            #these terms where not included in CC2
            "+ 1/24[[[[H,T2g],T1],T1],T1] + 1/24[[[[H,T1],T2g],T1],T1] + 1/24[[[[H,T1],T1],T2g],T1] + 1/24[[[[H,T1],T1],T1],T2g]"\
        ")*C0>"             
    )
    #T2 projection
    LAG_A2.append(
        "<C0^+*(LAM2g)*("\
            "1/24[[[[H,T1],T1],T1],T1]"\
        ")*C0>"             
    )
if not 0<lag_type<5 :
    raise Exception("MRCCPT_split unknown lagrangian type\nlag_type="+str(lag_type))


#LAG_E.append("<C0^+*(T1^+)*O1*C0>")
#LAG_E.append("<C0^+*(T2g^+)*O2g*C0>")



for item in LAG_E.show():
    print item
print "LAG_E finished"
LAG_E.set_rule()

for item in LAG_A1.show():
    print item
print "LAG_A1 finished"
LAG_A1.set_rule()

for item in LAG_A2.show():
    print item
print "LAG_A2 finished"
LAG_A2.set_rule()







FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG',
        LABEL_IN:'FORM_PT_LAG',
        INTERM:'FORM_GAM0'})

mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

debug_FORM('FORM_PT_LAG')



FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG_A1',
        LABEL_IN:'FORM_PT_LAG_A1',
        INTERM:'FORM_GAM0'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A1',
        LABEL_RES:'FORM_PT_LAG_A1'})

debug_FORM('FORM_PT_LAG_A1')


FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG_A2',
        LABEL_IN:'FORM_PT_LAG_A2',
        INTERM:'FORM_GAM0'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A2',
        LABEL_RES:'FORM_PT_LAG_A2'})

debug_FORM('FORM_PT_LAG_A2')





#Make the Derivative with respect to LAM  
DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A1',
        LABEL_RES:'FORM_PT_LAG_Amp1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A2',
        LABEL_RES:'FORM_PT_LAG_Amp2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

debug_FORM('FORM_PT_LAG_Amp1')

debug_FORM('FORM_PT_LAG_Amp2')



OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG2',
        LABELS_IN:['FORM_PT_LAG_Amp2','FORM_PT_LAG_Amp1','FORM_PT_LAG']})

