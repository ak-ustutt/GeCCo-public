"""Simple implementation of icMRCC lagragian, used to generate and control equations and formula

History:

Based on icMRCClag.py

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

i_am="ITFlag.py"

new_target('DEF_FORM_ITF_LAG')
heading('Defining the icMRCC Lagrangian used in ITF translator')

depend('DEF_T2g')
depend('DEF_LAM2g')
depend('DEF_O2g')

DEF_SCALAR({
        LABEL:'MRCC_LAG'})

DEF_ME_LIST({
        LIST:'MRCC_LAG_LST',
        OPERATOR:'MRCC_LAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

DEF_SCALAR({
        LABEL:'MRCC_LAG_A2'})


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")

LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))

# Use en_type/res_type to generate icMRCCSD equations
# Default is linear icMRCCSD
if keywords.is_keyword_set('method.ITF.en_type'):
    if (keywords.get('method.ITF.en_type') == '1'):
        LAG_E.append(_refexp("[H,T2g]"))
    elif(keywords.get('method.ITF.en_type') == '2'):
        LAG_E.append(_refexp("[H,T2g]"))
        LAG_E.append(_refexp("[[H,T2g],T2g]"))
    elif(keywords.get('method.ITF.en_type') == '3'):
        LAG_E.append(_refexp("[H,T2g]"))
        LAG_E.append(_refexp("[[H,T2g],T2g]"))
        LAG_E.append(_refexp("[[[H,T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.en_type') == '4'):
        LAG_E.append(_refexp("[H,T2g]"))
        LAG_E.append(_refexp("[[H,T2g],T2g]"))
        LAG_E.append(_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_E.append(_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
    else:
        raise Exception(i_am+": unrecognised value for en_type, must be {1,2,3,4}")
else:
    LAG_E.append(_refexp("[H,T2g]"))

if keywords.is_keyword_set('method.ITF.res_type'):
    if (keywords.get('method.ITF.res_type') == '1'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
    elif(keywords.get('method.ITF.res_type') == '2'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '3'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '4'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '5'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[H,T2g],T2g],T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '6'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[H,T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '7'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[H,T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g],T2g]"))
    elif(keywords.get('method.ITF.res_type') == '8'):
        LAG_A2.append(_L2_refexp("[H,T2g]"))
        LAG_A2.append(_L2_refexp("[[H,T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[H,T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[H,T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[H,T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g],T2g]"))
        LAG_A2.append(_L2_refexp("[[[[[[[[H,T2g],T2g],T2g],T2g],T2g],T2g],T2g],T2g]"))
    else:
        raise Exception(i_am+": unrecognised value for res_type, must be {1,2,3,4,5,6,7,8}")
else:
    LAG_A2.append(_L2_refexp("[H,T2g]"))

print("en_type: ", keywords.get('method.ITF.en_type'))
print("res_type: ", keywords.get('method.ITF.res_type'))


LAG_E.set_rule()
LAG_A2.set_rule()


DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_E']})


# Translate optmised formulae into ITF algo code
TRANSLATE_ITF({
        LABEL:'FOPT_MRCC_LAG',
        OUTPUT:'gecco.itfaa',
        TITLE:'gecco.formulae'})
