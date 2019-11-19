"""An implementation of the icMRCCSD Lagrangian, done term by term and separated T1 and T2.



History:

Yuri august 2017: Creation based on MRCC2lag.py. Implementation up to maxcom=2

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import ref_relaxation

i_am="icMRCCSDlag.py"

new_target('DEF_FORM_MRCC_LAG')
heading('Defining the icMRCC Lagrangian')

depend('DEF_T')
#depend('DEF_T2g')
#depend('DEF_T1')

depend('DEF_LAM')
#depend('DEF_LAM2g')
#depend('DEF_LAM1')

depend('DEF_O')
#depend('DEF_O2g')
#depend('DEF_O1')

depend('GAM0_CALC')

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


CLONE_OPERATOR({LABEL:'INTpp',TEMPLATE:'T2g'})
DEF_ME_LIST({LIST:'ME_INTpp',OPERATOR:'INTpp',IRREP:1,'2MS':0,AB_SYM:+1})


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")

def _L3_refexp(x):
    return _refexp("-LAM2g(" + x + ")")

def _L4_refexp(x):
    return _refexp("-LAM1(" + x + ")")

LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
LAG_A1.append(_L4_refexp("H"))

LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))
LAG_A2.append(_L3_refexp("H"))

#LAG_E.append(_refexp("[H,T1]"))
#LAG_E.append(_refexp("[H,T2g]"))

#LAG_A1.append(_L1_refexp("[H,T1]"))
#LAG_A1.append(_L1_refexp("[H,T2g]"))
#LAG_A1.append(_L1_refexp("(1/2)*[[H,T1],T2g]"))

#LAG_A2.append(_L2_refexp("[H,T1]"))
#LAG_A2.append(_L2_refexp("[H,T2g]"))
#LAG_A2.append(_L2_refexp("(1/2)*[[H,T2g],T2g]"))


#LAG_E.append(_refexp("(H*T1)+(H*T2g)"))
#LAG_E.append(_refexp("-(T1*H)-(T2g*H)"))
#
#LAG_E.append(_refexp("(1/2)*(H*T1 *T1 )"))
#LAG_E.append(_refexp("(1/2)*(H*T1 *T2g)"))
#LAG_E.append(_refexp("(1/2)*(H*T2g*T1 )"))
#LAG_E.append(_refexp("(1/2)*(H*T2g*T2g)"))
#
#LAG_E.append(_refexp("-(T1 *H*T1 )"))
#LAG_E.append(_refexp("-(T1 *H*T2g)"))
#LAG_E.append(_refexp("-(T2g*H*T1 )"))
#LAG_E.append(_refexp("-(T2g*H*T2g)"))
#
#LAG_E.append(_refexp("(1/2)*(T1 *T1 *H)"))
#LAG_E.append(_refexp("(1/2)*(T1 *T2g*H)"))
#LAG_E.append(_refexp("(1/2)*(T2g*T1 *H)"))
#LAG_E.append(_refexp("(1/2)*(T2g*T2g*H)"))
#
#
#
#LAG_A1.append(_L1_refexp("(H*T1)+(H*T2g)"))
#LAG_A1.append(_L1_refexp("-(T1*H)-(T2g*H)"))
#
#LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T1 )"))
#LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T2g)"))
#LAG_A1.append(_L1_refexp("(1/2)*(H*T2g*T1 )"))
#LAG_A1.append(_L1_refexp("(1/2)*(H*T2g*T2g)"))
#
#LAG_A1.append(_L1_refexp("-(T1 *H*T1 )"))
#LAG_A1.append(_L1_refexp("-(T1 *H*T2g)"))
#LAG_A1.append(_L1_refexp("-(T2g*H*T1 )"))
#LAG_A1.append(_L1_refexp("-(T2g*H*T2g)"))
#
#LAG_A1.append(_L1_refexp("(1/2)*(T1 *T1 *H)"))
#LAG_A1.append(_L1_refexp("(1/2)*(T1 *T2g*H)"))
#LAG_A1.append(_L1_refexp("(1/2)*(T2g*T1 *H)"))
#LAG_A1.append(_L1_refexp("(1/2)*(T2g*T2g*H)"))
#
#
#
#LAG_A2.append(_L2_refexp("(H*T1)+(H*T2g)"))
#LAG_A2.append(_L2_refexp("-(T1*H)-(T2g*H)"))
#
#LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T1 )"))
#LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T2g)"))
#LAG_A2.append(_L2_refexp("(1/2)*(H*T2g*T1 )"))
#LAG_A2.append(_L2_refexp("(1/2)*(H*T2g*T2g)"))
#
#LAG_A2.append(_L2_refexp("-(T1 *H*T1 )"))
#LAG_A2.append(_L2_refexp("-(T1 *H*T2g)"))
#LAG_A2.append(_L2_refexp("-(T2g*H*T1 )"))
#LAG_A2.append(_L2_refexp("-(T2g*H*T2g)"))
#
#LAG_A2.append(_L2_refexp("(1/2)*(T1 *T1 *H)"))
#LAG_A2.append(_L2_refexp("(1/2)*(T1 *T2g*H)"))
#LAG_A2.append(_L2_refexp("(1/2)*(T2g*T1 *H)"))
#LAG_A2.append(_L2_refexp("(1/2)*(T2g*T2g*H)"))



LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()

# E ============================
EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,P,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,P,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,V,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,PP,VH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
                   OPERATORS:['C0^+','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   CONNECT:[2,3],
                   LABEL_DESCR:["3,,PV,HH"]})


# T:aacc
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_E',NEW:False,OP_RES:'MRCC_LAG',
#                   OPERATORS:['C0^+','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4],
#                   CONNECT:[2,3],
#                   LABEL_DESCR:["2,,HH,VV,", "3,,VV,HH"]})



# A1 ============================
EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,H,P"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,P,", "4,,P,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,P,", "4,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,P,", "4,,P,V"]})

#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   CONNECT:[2,3, 3,4],
#                   LABEL_DESCR:["2,,H,P,", "4,,PP,VH"]})



EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,V,P"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,V,P,", "4,,P,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,V,P,", "4,,P,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,4, 3,4],
                   LABEL_DESCR:["2,,V,P,", "4,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,4, 3,4],
                   LABEL_DESCR:["2,,V,P,", "4,,PP,VH"]})



EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,H,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,V,", "4,,P,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,V,", "4,,P,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,V,", "4,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,H,V,", "4,,PP,VH"]})



# T:ac
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,,V", "2,4,H,", "3,4,,P", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,V", "3,4,H,P", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,3,,V", "2,3,H,", "3,4,H,P", "2,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,3,,V", "2,3,,V", "2,4,H,", "3,4,,P", "3,5,,V"]})


# T:ec
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,", "2,4,,P", "3,4,H,", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,,P", "2,4,H,", "3,4,,P", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,", "2,4,,P", "3,4,HH,P", "1,5,,V"]})



# A2 ============================
EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,HH,PP"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,HH,PP,", "4,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 3,4],
                   LABEL_DESCR:["2,,HH,PP,", "4,,P,H"]})

#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   CONNECT:[2,3, 3,4],
#                   LABEL_DESCR:["2,,HH,PP,", "4,,PP,VH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   LABEL_DESCR:["1,3,,V", "2,3,HH,", "2,4,,PP", "3,4,H,", "4,5,,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   LABEL_DESCR:["1,3,,V", "2,3,H,P", "2,4,H,P", "3,4,,P", "4,5,,V"]})



# T:eeac
EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,VH,PP"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 2,4, 3,4],
                   LABEL_DESCR:["2,,VH,PP,", "4,,PP,VH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 2,4, 3,4],
                   LABEL_DESCR:["2,,VH,PP,", "4,,P,H"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 2,4, 3,4],
                   LABEL_DESCR:["2,,VH,PP,", "4,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 2,4, 3,4],
                   LABEL_DESCR:["2,,VH,PP,", "4,,P,V"]})

EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   CONNECT:[2,3, 2,4, 3,4],
                   LABEL_DESCR:["2,,VH,PP,", "4,,V,H"]})



#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,2,,V", "2,3,H,", "2,4,,PP", "3,4,H,", "4,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,2,,V", "2,3,,P", "2,4,H,P", "3,4,,P", "4,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,2,,V", "2,3,H,P", "2,4,,P", "3,4,H,P", "4,5,,V"]})




# T:eacc
EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','C0'],
                   IDX_SV   :[1, 2, 3, 4],
                   LABEL_DESCR:["2,,HH,PV"]})




# T:aacc
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','C0'],
#                   IDX_SV   :[1, 2, 3, 4],
#                   LABEL_DESCR:["2,,HH,VV"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   CONNECT:[2,3, 3,4],
#                   LABEL_DESCR:["2,,HH,VV,", "4,,PP,HH"]})
#
##EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
##                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
##                   IDX_SV   :[1, 2, 3, 4, 5],
##                   CONNECT:[2,3, 3,4],
##                   LABEL_DESCR:["2,,HH,VV,", "4,,P,H"]})
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,VV", "2,4,H,","3,4,,P","1,5,,V"]})
#
##EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
##                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
##                   IDX_SV   :[1, 2, 3, 4, 5],
##                   CONNECT:[2,3, 3,4],
##                   LABEL_DESCR:["2,,HH,VV,", "4,,VV,HH"]})



# Doesn't exist?
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,,HH,PP,", "4,,P,V"]})






# T:eecc
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,", "2,4,H,PP", "3,4,H,", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,,P", "2,4,HH,P", "3,4,,P", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,HH", "2,4,,PP", "3,4,HH,", "1,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,3,,V", "2,3,H,", "2,4,H,PP", "3,4,H,", "3,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,3,,V", "2,3,,P", "2,4,HH,P", "3,4,,P", "3,5,,V"]})
#
#EXPAND_OP_PRODUCT({LABEL:'FORM_MRCC_LAG_A2',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,H,P", "2,4,H,P", "3,4,H,P", "1,5,,V"]})


EXPAND_OP_PRODUCT({LABEL:'CCD_FORM_PP',NEW:True,OP_RES:'INTpp',
                   OPERATORS:['INTpp','H','T2g','INTpp'],
                   IDX_SV   :[1,2,3,1],
                   LABEL_DESCR:["2,,PP,PP","3,,PP,HH"]})

EXPAND_OP_PRODUCT({LABEL:'CCD_FORM_PP',NEW:False,OP_RES:'INTpp',
                   OPERATORS:['INTpp','H','T2g','INTpp'],
                   IDX_SV   :[1,2,3,1],
                   LABEL_DESCR:["2,,PP,PP","3,,PP,VH"]})

PRINT_FORMULA({LABEL:'CCD_FORM_PP',MODE:'SHORT'})

# K4E
# T:eecc
#EXPAND_OP_PRODUCT({LABEL:'K4E',NEW:True,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["2,3,,PP", "2,4,HH,", "3,4,,PP", "1,5,,V"]})

EXPAND_OP_PRODUCT({LABEL:'K4E',NEW:True,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   LABEL_DESCR:["2,3,H,PP", "2,4,H,", "3,4,,P", "1,5,,V"]})

# T:eeac
#EXPAND_OP_PRODUCT({LABEL:'K4E',NEW:False,OP_RES:'MRCC_LAG_A2',
#                   OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,2,,V", "2,3,,PP", "2,4,H,", "3,4,,PP", "4,5,,V"]})

# R:eeac K3E
EXPAND_OP_PRODUCT({LABEL:'K4E',NEW:False,OP_RES:'MRCC_LAG_A2',
                   OPERATORS:['C0^+','LAM2g','H','T1','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   LABEL_DESCR:["1,2,,V", "2,3,,PP", "2,4,H,", "3,4,,P", "3,5,,V"]})

# K3E
# T:ec
EXPAND_OP_PRODUCT({LABEL:'K4E1',NEW:True,OP_RES:'MRCC_LAG_A1',
                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
                   IDX_SV   :[1, 2, 3, 4, 5],
                   LABEL_DESCR:["2,3,,P", "2,4,H,", "3,4,H,PP", "1,5,,V"]})

#EXPAND_OP_PRODUCT({LABEL:'K4E1',NEW:False,OP_RES:'MRCC_LAG_A1',
#                   OPERATORS:['C0^+','LAM1','H','T2g','C0'],
#                   IDX_SV   :[1, 2, 3, 4, 5],
#                   LABEL_DESCR:["1,3,,V", "2,3,,P", "2,4,H,", "3,4,,PP", "4,5,,V"]})




PRINT_FORMULA({LABEL:'K4E',MODE:'SHORT'})
FACTOR_OUT({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',INTERM:'K4E'})

PRINT_FORMULA({LABEL:'K4E1',MODE:'SHORT'})
FACTOR_OUT({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',INTERM:'K4E1'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})


#FACTOR_OUT({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',INTERM:'CCD_FORM_PP'})
#PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})


#comment("debug form MRCC_LAG_E")
#debug_FORM('FORM_MRCC_LAG_E', True,mode='LONG')

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_E',
        LABEL_RES:'FORM_MRCC_LAG_E',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'FORM_GAM0'})


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
        #LABELS_IN:['CCD_FORM_PP','FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})
        LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})


TRANSLATE_ITF({
        LABEL:'FOPT_MRCC_LAG',
        OUTPUT:'icmrcc_mrccsd_11_doublet.itfaa',
        TITLE:'icmrcc_mrccsd_11_doublet.formulae',
        MULTI:True,
        PROCESS:True,
        KEXT:False,
        INIT_RES:True,
        ITIN:True})

#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
