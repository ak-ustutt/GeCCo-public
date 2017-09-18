
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

i_am="MRCC2lag.py"

lagrangian = keywords.get('method.MRCC2.lagrangian')
lag_type = int(lagrangian) if lagrangian is not None else 4

known_hamiltonians=["DYALL","REPT","F_EFF","FULL"]
hamiltonian = keywords.get('method.MRCC2.hamiltonian')
hamiltonian=str(hamiltonian).strip() if hamiltonian is not None else "DYALL"

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

depend('GAM0_CALC')
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

LAG_E  = stf.GenForm(label="FORM_PT_LAG_E",      OP_res="PT_LAG") # energy part of the lagrangian
LAG_A1 = stf.GenForm(label="FORM_PT_LAG_A1_RAW", OP_res="PT_LAG_A1") # single excitation equations 
LAG_A2 = stf.GenForm(label="FORM_PT_LAG_A2_RAW", OP_res="PT_LAG_A2") # double excitation equations

if lag_type >= 1 :
    LAG_E +=  "<C0^+*("\
                    "H"\
                    "+(H*T1)"\
                    "+(H*T2g)"\
                ")*C0>"
    LAG_A1 += "<C0^+*(LAM1)*("\
                    "H"\
                    "+[H,T1]"\
                    "+[H,T2g]"\
              ")*C0>"
              
    LAG_A2 += "<C0^+*(LAM2g)*("\
                    "H"\
                    "+[H,T1]"\
              ")*C0>"      
              
    if hamiltonian=="DYALL":
        LAG_A2.append("<C0^+*(LAM2g)*([HAM_D,T2g])*C0>")
    elif hamiltonian=="REPT":
        LAG_A2.append("<C0^+*(LAM2g)*([REPT_HAM,T2g])*C0>")
    elif hamiltonian=="F_EFF":
        LAG_A2.append("<C0^+*(LAM2g)*([FOCK_EFF,T2g])*C0>")
    elif hamiltonian=="FULL":
        LAG_A2.append("<C0^+*(LAM2g)*([H,T2g])*C0>")

if lag_type >= 2 :
#no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
    LAG_E +="<C0^+*("\
                    "1/2((H*T1)*T1)"\
                    "+1/2((H*T2g)*T1)+1/2((H*T1)*T2g)"\
            ")*C0>"

    LAG_A1 += "<C0^+*(LAM1)*("\
                    "1/2[[H,T1],T1]"\
                    "+1/2[[H,T2g],T1]+1/2[[H,T1],T2g]"\
              ")*C0>"

    LAG_A2 += "<C0^+*(LAM2g)*("\
                    "1/2[[H,T1],T1]"\
              ")*C0>"             

if lag_type >= 3 :
#no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
    LAG_E += "<C0^+*("\
                    "1/6(((H*T1)*T1)*T1)"\
                    "+1/6(((H*T2g)*T1)*T1)+1/6(((H*T1)*T2g)*T1)+1/6(((H*T1)*T1)*T2g)"\
             ")*C0>"
          
    LAG_A1 += "<C0^+*(LAM1)*("\
                    "1/6[[[H,T1],T1],T1]"\
                    "+1/6[[[H,T2g],T1],T1]+1/6[[[H,T1],T2g],T1]+1/6[[[H,T1],T1],T2g]"\
              ")*C0>"

    LAG_A2 += "<C0^+*(LAM2g)*("\
                    "1/6[[[H,T1],T1],T1]"\
              ")*C0>"             

if lag_type >= 4 :
#no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
#these terms where not included in SRCC2
    LAG_E += "<C0^+*("\
                    "1/24((((H*T1)*T1)*T1)*T1)"\
                    "+ 1/24((((H*T2g)*T1)*T1)*T1) + 1/24((((H*T1)*T2g)*T1)*T1) + 1/24((((H*T1)*T1)*T2g)*T1) + 1/24((((H*T1)*T1)*T1)*T2g)"\
             ")*C0>"             

    LAG_A1 += "<C0^+*(LAM1)*("\
                    "1/24[[[[H,T1],T1],T1],T1]"\
                    "+ 1/24[[[[H,T2g],T1],T1],T1] + 1/24[[[[H,T1],T2g],T1],T1] + 1/24[[[[H,T1],T1],T2g],T1] + 1/24[[[[H,T1],T1],T1],T2g]"\
              ")*C0>"

    LAG_A2 += "<C0^+*(LAM2g)*("\
                    "1/24[[[[H,T1],T1],T1],T1]"\
              ")*C0>"
if not 0<lag_type<5 :
    raise Exception("MRCC2 unknown lagrangian type\nlag_type="+str(lag_type))


LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()


REORDER_FORMULA({LABEL_RES:"FORM_PT_LAG_E",
                 LABEL_IN:"FORM_PT_LAG_E"})


debug_FORM('FORM_PT_LAG_E')



FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG_INT',
        LABEL_IN:'FORM_PT_LAG_E',
        INTERM:'FORM_GAM0'})



mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_INT',
        LABEL_RES:'FORM_PT_LAG'})

debug_FORM('FORM_PT_LAG')



FACTOR_OUT({
        LABEL_IN:'FORM_PT_LAG_A1_RAW',
        LABEL_RES:'FORM_PT_LAG_A1_INT',
        INTERM:'FORM_GAM0'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A1_INT',
        LABEL_RES:'FORM_PT_LAG_A1'})

debug_FORM('FORM_PT_LAG_A1')


FACTOR_OUT({
        LABEL_IN:'FORM_PT_LAG_A2_RAW',
        LABEL_RES:'FORM_PT_LAG_A2_INT',
        INTERM:'FORM_GAM0'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A2_INT',
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

REORDER_FORMULA({LABEL_RES:"FORM_PT_LAG_Amp2",
                     LABEL_IN:"FORM_PT_LAG_Amp2"})
REORDER_FORMULA({LABEL_RES:"FORM_PT_LAG_Amp1",
                     LABEL_IN:"FORM_PT_LAG_Amp1"})

OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG2',
        LABELS_IN:['FORM_GAM0','FORM_PT_LAG_Amp2','FORM_PT_LAG_Amp1','FORM_PT_LAG']})

#----------------------------------------------------------------------
# formula for reference relaxiation  
new_target("FOPT_OMG_C0")
depend("DEF_FORM_PT_LAG2")

depend("DEF_A_C0")



DERIVATIVE({
    LABEL_IN:'FORM_PT_LAG_E',
    LABEL_RES:'FORM_A_C0',
    OP_RES:'A_C0',
    OP_DERIV:'C0^+'})

debug_FORM('FORM_A_C0')

OPTIMIZE({
        LABEL_OPT:'FOPT_OMG_C0',
        LABELS_IN:'FORM_A_C0'})


