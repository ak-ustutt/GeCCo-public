#
# An implementation of the icMRCC Lagrangian, with term by term analysis.
# 
#
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

i_am="icMRCCSDlag.py"



new_target('DEF_FORM_MRCC_LAG')

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


LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _refexp("LAM1*H"))
LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _refexp("LAM2g*H"))



LAG_E.append(_refexp("(H*T1)+(H*T2g)"))

LAG_A1.append(_refexp("(LAM1*H*T1)+(LAM1*H*T2g)"))
LAG_A1.append(_refexp("(LAM1*T1*H)+(LAM1*T2g*H)"))

LAG_A1.append(_refexp("(LAM2g*H*T1)+(LAM2g*H*T2g)"))
LAG_A1.append(_refexp("(LAM2g*T1*H)+(LAM2g*T2g*H)"))




LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()


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
        LABEL_OPT:'FOPT_MRCC_LAG2',
        LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})


