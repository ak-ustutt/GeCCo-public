
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import ref_relaxation

i_am="MRCC2lag.py"

lagrangian = keywords.get('method.MRCC2.lagrangian')
lag_type = int(lagrangian) if lagrangian is not None else 4

lagrangian = keywords.get('method.MRCC2.maxcom_en')
maxcom_en = int(lagrangian) if lagrangian is not None else lag_type
lagrangian = keywords.get('method.MRCC2.maxcom_res1')
maxcom_res1 = int(lagrangian) if lagrangian is not None else lag_type
lagrangian = keywords.get('method.MRCC2.maxcom_res2')
maxcom_res2 = int(lagrangian) if lagrangian is not None else lag_type

MRCCSD_mode = keywords.is_keyword_set('method.MRCC2.MRCCSD_mode')

lin_psD = keywords.get('method.MRCC2.lin_psD')
lin_psD = int(lin_psD) if lin_psD is not None else 0

known_hamiltonians=["DYALL","REPT","F_EFF","FULL"]
hamiltonian = keywords.get('method.MRCC2.hamiltonian')
hamiltonian=str(hamiltonian).strip() if hamiltonian is not None else "DYALL"



if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

#--------------------------------------------------------------------------------------#
#Define the lagrangian
#--------------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG2')

depend('DEF_T')
#depend('DEF_T2g')
#depend('DEF_T1')

depend('DEF_LAM')
#depend('DEF_LAM2g')
#depend('DEF_LAM1')

depend('DEF_O')
#depend('DEF_O2g')
#depend('DEF_O1')

depend('DEF_ToX')

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

DEF_SCALAR({
        LABEL:'PT_LAG_A'})





def create_lag_E(label, OP_res, maxcom, MRCCSD_mode = False):
    LAG_E  = stf.GenForm(label=label, OP_res=OP_res) # energy part of the lagrangian
    #no commutator necessary since T*H has open hole or particle lines (if T is not purely active)
    if maxcom >= 1 :
        LAG_E +=  "<C0^+*("\
                      "H"\
                      "+(H*To0)"\
                      "+(H*To1)"\
                  ")*C0>"
    if maxcom >= 2 :
        LAG_E += "<C0^+*("\
                      "1/2((H*To0)*To0)"\
                      "+1/2((H*To1)*To0)+1/2((H*To0)*To1)"\
                 ")*C0>"
        if MRCCSD_mode:
            LAG_E += "<C0^+*("\
                         "1/2((H*To1)*To1)"\
                     ")*C0>"
    if maxcom >= 3 :
        if MRCCSD_mode:
            LAG_E += "<C0^+*("\
                         "1/6( ( ( H*(To0+To1) )*(To0+To1) )*(To0+To1) )"\
                     ")*C0>"
        else:
            LAG_E += "<C0^+*("\
                         "1/6(((H*To0)*To0)*To0)"\
                         "+1/6(((H*To1)*To0)*To0)+1/6(((H*To0)*To1)*To0)+1/6(((H*To0)*To0)*To1)"\
                     ")*C0>"
    if maxcom >= 4 :
        if MRCCSD_mode:
            LAG_E += "<C0^+*("\
                         "1/24( ( ( ( H*(To0+To1) )*(To0+To1) )*(To0+To1) )*(To0+To1) )"\
                     ")*C0>"
        else:
            LAG_E += "<C0^+*("\
                          "1/24((((H*To0)*To0)*To0)*To0)"\
                          "+ 1/24((((H*To1)*To0)*To0)*To0) + 1/24((((H*To0)*To1)*To0)*To0) + 1/24((((H*To0)*To0)*To1)*To0) + 1/24((((H*To0)*To0)*To0)*To1)"\
                     ")*C0>"             
    if not 0<maxcom<5 :
        raise Exception("MRCC2 unknown maxcommutator for energy\nmaxcom="+str(maxcom))
    return LAG_E

def create_lag_A1(label, OP_res, maxcom, MRCCSD_mode = False, lin_psD = 0):
    LAG_A1  = stf.GenForm(label=label, OP_res=OP_res) # residual part of the lagrangian
    if lin_psD > 0 and maxcom > 2 :
         raise Exception("linearized pseudo doubles not yet implemented for maxcom > 2!")
    if lin_psD > 0 and MRCCSD_mode :
         raise Exception("contradictory flags: lin_psD and MRCCSD_mode")
    if maxcom >= 1 :
        LAG_A1 +=  "<C0^+*(LAMo0)*("\
                      "H"\
                      "+[H,To0]"\
                      "+[H,To1]"\
                  ")*C0>"
    if maxcom >= 2 :
        if lin_psD == 0:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                      "1/2[[H,To0],To0]"\
                      "+1/2[[H,To1],To0]+1/2[[H,To0],To1]"\
                  ")*C0>"
        elif lin_psD == 1:
            LAG_A1 += "<C0^+*(LAM1)*("\
                      "1/2[[H,To0],To0]"\
                      "+1/2[[H,To1],To0]+1/2[[H,To0],To1]"\
                  ")*C0>"
            LAG_A1 += "<C0^+*(LAM2ps)*("\
                      " 1/2[[H,T1],T1]"\
                      "+1/2[[H,T2ps],T1]+1/2[[H,T1],T2ps]"\
                      "+1/2[[H,To1],T1]+1/2[[H,T1],To1]"\
                  ")*C0>"
        elif lin_psD == 2:
            LAG_A1 += "<C0^+*(LAM1)*("\
                      " 1/2[[H,T1],T1]"\
                      "+1/2[[H,T1],T2ps]+1/2[[H,T2ps],T1]"\
                      "+1/2[[H,T1],To1]+1/2[[H,To1],T1]"\
                  ")*C0>"
            LAG_A1 += "<C0^+*(LAM2ps)*("\
                      " 1/2[[H,T1],T1]"\
                      "+1/2[[H,T2ps],T1]+1/2[[H,T1],T2ps]"\
                      "+1/2[[H,To1],T1]+1/2[[H,T1],To1]"\
                  ")*C0>"
        if MRCCSD_mode:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                          "1/2[[H,To1],To1]"\
                      ")*C0>"
  
    if maxcom >= 3 :
        if MRCCSD_mode:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                         "1/6[[[H,To0+To1],To0+To1],To0+To1]"\
                      ")*C0>"

        else:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                          "1/6[[[H,To0],To0],To0]"\
                          "+1/6[[[H,To1],To0],To0]+1/6[[[H,To0],To1],To0]+1/6[[[H,To0],To0],To1]"\
                      ")*C0>"
    if maxcom >= 4 :
        if MRCCSD_mode:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                          "1/24[[[[H,To0+To1],To0+To1],To0+To1],To0+To1]"\
                      ")*C0>"             
        else:
            LAG_A1 += "<C0^+*(LAMo0)*("\
                          "1/24[[[[H,To0],To0],To0],To0]"\
                          "+ 1/24[[[[H,To1],To0],To0],To0] + 1/24[[[[H,To0],To1],To0],To0] + 1/24[[[[H,To0],To0],To1],To0] + 1/24[[[[H,To0],To0],To0],To1]"\
                     ")*C0>"             
    if not 0<maxcom<5 :
        raise Exception("MRCC2 unknown lagrangian type\nmaxcom="+str(maxcom))
    return LAG_A1

def create_lag_A2(label, OP_res, maxcom, hamiltonian, MRCCSD_mode = False, lin_psD = 0):
    LAG_A2  = stf.GenForm(label=label, OP_res=OP_res) # residual part of the lagrangian
    if lin_psD > 0 and maxcom > 2 :
         raise Exception("linearized pseudo doubles not yet implemented for maxcom > 2!")
    if lin_psD > 0 and MRCCSD_mode :
         raise Exception("contradictory flags: lin_psD and MRCCSD_mode")
    if maxcom >= 1 :
        LAG_A2 +=  "<C0^+*(LAMo1)*("\
                      "H"\
                      "+[H,To0]"\
                  ")*C0>"
    if hamiltonian=="DYALL":
        LAG_A2.append("<C0^+*(LAMo1)*([HAM_D,To1])*C0>")
    elif hamiltonian=="REPT":
        LAG_A2.append("<C0^+*(LAMo1)*([REPT_HAM,To1])*C0>")
    elif hamiltonian=="F_EFF":
        LAG_A2.append("<C0^+*(LAMo1)*([FOCK_EFF,To1])*C0>")
    elif hamiltonian=="FULL":
        LAG_A2.append("<C0^+*(LAMo1)*([H,To1])*C0>")
        
    if maxcom >= 2 :
        if MRCCSD_mode:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/2[[H,To0+To1],To0+To1]"\
                      ")*C0>"
        elif lin_psD == 0:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/2[[H,To0],To0]"\
                      ")*C0>"
        elif lin_psD >= 1:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/2[[H,T1],T1]"\
                      ")*C0>"

    if maxcom >= 3 :
        if MRCCSD_mode:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/6[[[H,To0+To1],To0+To1],To0+To1]"\
                      ")*C0>"
        else:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/6[[[H,To0],To0],To0]"\
                      ")*C0>"
    if maxcom >= 4 :
        if MRCCSD_mode:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/24[[[[H,To0+To1],To0+To1],To0+To1],To0+To1]"\
                      ")*C0>"
        else:
            LAG_A2 += "<C0^+*(LAMo1)*("\
                          "1/24[[[[H,To0],To0],To0],To0]"\
                      ")*C0>"             
    if not 0<maxcom<5 :
        raise Exception("MRCC2 unknown lagrangian type\nmaxcom="+str(maxcom))
    return LAG_A2



    
create_lag_E("FORM_PT_LAG_E", "PT_LAG", maxcom_en, MRCCSD_mode).set_rule()
create_lag_A1("FORM_PT_LAG_A1_RAW", "PT_LAG_A1", maxcom_res1, MRCCSD_mode, lin_psD).set_rule()
create_lag_A2("FORM_PT_LAG_A2_RAW", "PT_LAG_A2", maxcom_res2, hamiltonian, MRCCSD_mode, lin_psD).set_rule()

# joint A1 and A2 into a single expansion
DEF_FORMULA({LABEL:'FORM_PT_LAG_A_RAW',FORMULA:'PT_LAG_A=PT_LAG_A1+PT_LAG_A2'})
EXPAND({LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW',
        INTERM:['FORM_PT_LAG_A1_RAW','FORM_PT_LAG_A2_RAW']})

# -> replace by actual T and LAM operators ...
EXPAND({LABEL_IN:'FORM_PT_LAG_E',
        LABEL_RES:'FORM_PT_LAG_E',
        INTERM:'FORM_To0'})
debug_FORM('FORM_PT_LAG_E')

EXPAND({LABEL_IN:'FORM_PT_LAG_E',
        LABEL_RES:'FORM_PT_LAG_E',
        INTERM:'FORM_To1'})
debug_FORM('FORM_PT_LAG_E')

REPLACE({LABEL_RES:'FORM_PT_LAG_A_RAW',
         LABEL_IN:'FORM_PT_LAG_A_RAW',
         OP_LIST:['T2ps','T2g']})
REPLACE({LABEL_RES:'FORM_PT_LAG_A_RAW',
         LABEL_IN:'FORM_PT_LAG_A_RAW',
         OP_LIST:['LAM2ps','LAM2g']})
EXPAND({LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW',
        INTERM:'FORM_To0'})
EXPAND({LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW',
        INTERM:'FORM_To1'})
EXPAND({LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW',
        INTERM:'FORM_LAMo0'})
EXPAND({LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW',
        INTERM:'FORM_LAMo1'})
SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_RAW'})
debug_FORM('FORM_PT_LAG_A_RAW',only_this=True)


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
        LABEL_IN:'FORM_PT_LAG_A_RAW',
        LABEL_RES:'FORM_PT_LAG_A_INT',
        INTERM:'FORM_GAM0'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A_INT',
        LABEL_RES:'FORM_PT_LAG_A'})

debug_FORM('FORM_PT_LAG_A')



#Make the Derivative with respect to LAM  
DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_LAG_Amp1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A',
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


#-----
ref_relaxation.make_form_for_optref_minus3('FORM_PT_LAG_E', 'DEF_FORM_PT_LAG2')

