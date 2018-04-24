"""Implementation of icMRCC lagragian, used to generate and control equations and formula



History:

Based on icMRCClag.py

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import ref_relaxation

i_am="ITFlag.py"

new_target('DEF_FORM_MRCC_LAG')
heading('Defining the icMRCC Lagrangian')

depend('DEF_T2g')
depend('DEF_T1')

depend('DEF_LAM2g')
depend('DEF_LAM1')

depend('DEF_O2g')
depend('DEF_O1')

DEF_SCALAR({
        LABEL:'MRCC_LAG'})

DEF_ME_LIST({
        LIST:'MRCC_LAG_LST',
        OPERATOR:'MRCC_LAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

DEF_SCALAR({
        LABEL:'MRCC_LAG_A1'})

DEF_SCALAR({
        LABEL:'MRCC_LAG_A2'})


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")


LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))



LAG_E.append(_refexp("[H,T1]"))
LAG_E.append(_refexp("[H,T2g]"))


LAG_A1.append(_L1_refexp("[H,T1]"))
LAG_A1.append(_L1_refexp("[H,T2g]"))


LAG_A2.append(_L2_refexp("[H,T1]"))
LAG_A2.append(_L2_refexp("[H,T2g]"))


LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()


#comment("debug form MRCC_LAG_E")
#debug_FORM('FORM_MRCC_LAG_E', True,mode='LONG')

#Make the Derivative with respect to LAM  
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_Amp1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})
                     
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})


OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})

TRANSLATE_ITF({
        LABEL:'FOPT_MRCC_LAG',
        OUTPUT:'Hello.txt'})


#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
