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

#------------------------------------------------------------------------------#
#Define the MRCCPT(2) lagrangian
#------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG')

depend('T-Operators')
depend('DEF_LAM')
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
        LABEL:'PT_LAG_A'})


#Energy equation no 
LAG_E=stf.Formula("FORM_PT_LAG:PT_LAG="\
                  "<C0^+*(H+H*T2_ca)C0>")

LAG_A=stf.Formula("FORM_PT_LAG_A:PT_LAG="\
                  "<C0^+*(LAMges)*(H)C0>")



if hamiltonian=="DYALL":
    LAG_A.append("<C0^+*(LAMges)*([HAM_D,T2_ca])*C0>")
elif hamiltonian=="REPT":
    LAG_A.append("<C0^+*(LAMges)*([REPT_HAM,T2_ca])*C0>")
elif hamiltonian=="F_EFF":
    LAG_A.append("<C0^+*(LAMges)*([FOCK_EFF,T2_ca])*C0>")

LAG_E.append("<C0^+*(T2_ca^+)*Oges*C0>")

for item in LAG_E.show():
    print item
LAG_E.set_rule()

comment("LAG_E finished")

for item in LAG_A.show():
    print item
LAG_A.set_rule()
comment("LAG_A finished")


mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_LAG_A'})



debug_FORM('FORM_PT_LAG')



#Make the Derivative with respect to LAMges.
DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_Amp',
        OP_RES:'Oges',
        OP_DERIV:'LAMges'})



debug_FORM('FORM_PT_Amp')

# For some reason the reorder formulas are necessary. 


OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG',
        LABELS_IN:['FORM_PT_Amp','FORM_PT_LAG']})

