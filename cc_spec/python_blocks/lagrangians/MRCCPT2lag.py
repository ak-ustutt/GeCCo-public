from gecco_interface import *
from gecco_modules.NoticeUtil import *
import gecco_modules.string_to_form as stf

i_am="MRCCPT2lag.py"

known_hamiltonians={"DYALL","REPT","F_EFF"}
hamiltonian="DYALL"
if keywords.is_keyword_set('method.MRCCPT2.hamiltonian'):
    hamiltonian=str(keywords.get('method.MRCCPT2.hamiltonian')).strip()
print("hamiltonian: ", hamiltonian, type(hamiltonian))

if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

third_ord_energy=False
if keywords.is_keyword_set('method.MRCCPT2.3rd_E'):
    if (keywords.get('method.MRCCPT2.3rd_E') == "T"):
        third_ord_energy=True
    elif(keywords.get('method.MRCCPT2.3rd_E') == "F"):
        third_ord_energy=False
    else :
        raise Exception(i_am+": unrecognised value for opion 3rd_E:"+str(hamiltonian))
print("3rd_E ", third_ord_energy, type(third_order_energy))

#------------------------------------------------------------------------------#
#Define the MRCCPT(2) lagrangian
#------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG')


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

#dummy for residue of amplitude equation
DEF_SCALAR({
        LABEL:'PT_LAG_A'})


#Energy equation no 
LAG_E=stf.Formula("FORM_PT_LAG:PT_LAG="\
                  "<C0^+*(H+H*T2g)C0>")

LAG_A=stf.Formula("FORM_PT_LAG_A:PT_LAG="\
                  "<C0^+*(LAM2g)*(H)C0>")

if hamiltonian=="DYALL":
    LAG_A.append("<C0^+*(LAM2g)*([HAM_D,T2g])*C0>")
elif hamiltonian=="REPT":
    LAG_A.append("<C0^+*(LAM2g)*([REPT_HAM,T2g])*C0>")
elif hamiltonian=="F_EFF":
    LAG_A.append("<C0^+*(LAM2g)*([FOCK_EFF,T2g])*C0>")


if (third_ord_energy):
    if hamiltonian=="DYALL":
        LAG_A.append("<C0^+*(T2g^+)*((H-HAM_D)*T2g)*C0>")
    elif hamiltonian=="REPT":
        LAG_A.append("<C0^+*(T2g^+)*((H-REPT_HAM)*T2g)*C0>")
    elif hamiltonian=="F_EFF":
        LAG_A.append("<C0^+*(T2g^+)*((H-FOCK_EFF)*T2g)*C0>")


# optional penality term
#LAG_E.append("<C0^+*(T2g^+)*O2g*C0>")




for item in LAG_E.show():
    print item
LAG_E.set_rule()

print("LAG_E finished")


for item in LAG_A.show():
    print item
LAG_A.set_rule()

print("LAG_A finished")


FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG',
        LABEL_IN:'FORM_PT_LAG',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG_A',
        LABEL_IN:'FORM_PT_LAG_A',
        INTERM:'FORM_GAM0'})



mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_LAG_A'})



debug_FORM('FORM_PT_LAG')


#Make the Derivative with respect to LAM2g.
DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_Amp',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

debug_FORM('FORM_PT_Amp')



OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG',
        LABELS_IN:['FORM_PT_Amp','FORM_PT_LAG']})

