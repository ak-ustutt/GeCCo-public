from gecco_interface import *
from gecco_modules.NoticeUtil import *
import gecco_modules.string_to_form as stf

i_am="MRCCPT2lag.py"

known_hamiltonians={"DYALL","REPT","F_EFF"}
hamiltonian="DYALL"
if keywords.is_keyword_set('method.MRCC2.hamiltonian'):
    hamiltonian=str(keywords.get('method.MRCC2.hamiltonian')).strip()
print("hamiltonian: ", hamiltonian, type(hamiltonian))

if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

#------------------------------------------------------------------------------#
#Define the MRCCPT(2) lagrangian
#------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG')

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

#Energy equation no 
LAG=stf.Formula("FORM_PT_LAG:PT_LAG="\
                "<C0^+*(H+H*T2_ca)C0>")

LAG.append("<C0^+*(T2_ca^+)*(H)C0>")

if hamiltonian=="DYALL":
    LAG.append("<C0^+*(T2_ca^+)*([HAM_D,T2_ca])*C0>")
elif hamiltonian=="REPT":
    LAG.append("<C0^+*(T2_ca^+)*([REPT_HAM,T2_ca])*C0>")
elif hamiltonian=="F_EFF":
    LAG.append("<C0^+*(T2_ca^+)*([F_EFF,T2_ca])*C0>")


for item in LAG.show():
    print item
LAG.set_rule()

mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})



debug_FORM('FORM_PT_LAG')



#Make the Derivative with respect to T2_ca^+ and factor that out. 
#PT_LAG=E+A => PT_LAG=E_1,2,3+T1_ca^+*O1+T2_ca^+*Oges
DERIVATIVE({LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_Ampl',
        OP_RES:'Oges',
        OP_DERIV:'T2_ca^+'})


debug_FORM('FORM_PT_Ampl_T1')

# For some reason the reorder formulas are necessary. 
REORDER_FORMULA({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

REORDER_FORMULA({
        LABEL_IN:'FORM_PT_Ampl',
        LABEL_RES:'FORM_PT_Ampl'})
FACTOR_OUT({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG',
        INTERM:'FORM_PT_Ampl'})


debug_FORM('FORM_PT_LAG')

debug_FORM('FORM_PT_Ampl')


OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG',
        LABELS_IN:['FORM_PT_Ampl','FORM_PT_LAG']})

