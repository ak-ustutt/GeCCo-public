#
# Set the multistate coupling term for the icMRCC lagrangian:
#
#  - sum_i sum_(j != i) 
# <0| C0_i^+ L_i exp(-T_i) exp(-T_j) C0_j |0> <0| C0_j^+ exp(-T_i) H exp(-T_i) C0_i |0>
#
# i and j refer to state indices. The uncoupled part of the Lagrangian is just the summation
# of the tradidional Lagrangian for each state.
# The truncation of the BCH expansion is controled by the following input keywords:
# method.MRCC.maxcom_res_cpl and method.MRCC.maxcom_en_cpl
# The first for the L- and the second for the H-part. If not present, by
# method.MRCC.maxcom_res and method.MRCC.maxcom_en
#
# Yuri, Nov 2014
#
from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

import math

# Set some parameters fom input
#
multistate = inp.get('method.MR.multistate') == 'T'
ciroot = inp.get('method.MR.ciroot')
if (ciroot == None):
    ciroot = 1
else:
    ciroot = int(ciroot)

n_states = ciroot

if (not( multistate)):
    quit_error("MRCC_Lagrangian_coupl_term.py being called for multistate = F")
if (n_states < 2):
    quit_error("MRCC_Lagrangian_coupl_term.py being called for n_states < 2")

# maximum level in commutator expansion
maxcom_default = 2
maxcom_en_default = 2

maxcom = inp.get('method.MRCC.maxcom_res_cpl')
if (maxcom == None):
    maxcom = inp.get('method.MRCC.maxcom_res')
    if (maxcom == None):
        maxcom = maxcom_defaut
maxcom = int(maxcom)

maxcom_en = inp.get('method.MRCC.maxcom_en_cpl')
if (maxcom_en == None):
    maxcom_en = inp.get('method.MRCC.maxcom_en')
    if (maxcom_en == None):
        maxcom_en = maxcom_en_default
maxcom_en = int(maxcom_en)


# Set some variables for operator labels
#
lag_label = 'F_MRCC_LAG'
res = 'NORM'
top = 'T'
hop = 'H'
lop = 'L'
cop = 'C0'

# A function to set the avoid list, where n1 and n2
# are the number of operators: every operator
# from 1 to n1 avoids the operators from n1+1 to n1+n2
#
def set_avoid_list( n1, n2):
    avoid = []
    for i1 in range(1, n1+1):
        for i2 in range(n1+1, n1+n2+1):
            avoid.extend([i1, i2])
    return avoid

# The factor of the BCH expansion
def set_factor(n, k):
    fac = 1.0/(math.factorial(k)*math.factorial(n-k))
    if (k % 2 == 1):
        fac = -fac
    return fac

new_target('F_MRCC_LAG_coupl')
depend('F_MRCC_LAG')

for i_state in range(1, n_states+1):
    state_label = "" if (i_state == 1) else "_" + str(i_state)
    cdg_i = cop + state_label + '^+'
    lop_i = lop + state_label
    top_i = top + state_label
    cop_i = cop + state_label

    for j_state in range(1, n_states+1):
        if (j_state == i_state):
            continue

        state_label = "" if (j_state == 1) else "_" + str(j_state)
        cdg_j = cop + state_label + '^+'
        top_j = top + state_label
        cop_j = cop + state_label

        for nL in range(0, maxcom+1):
            nop_L = 3 + nL
            for nH in range(0, maxcom_en+1):
                nop_H = 3 + nH

                for kL in range(0, nL+1):
                    fac_L = set_factor(nL, kL)
                    for kH in range(0, nH+1):
                        fac_H = set_factor(nH, kH)

                        op_list = [cdg_i, lop_i] + [top_i]*kL  + [top_j]*(nL-kL) + [cop_j] +\
                                  [cdg_j] + [top_i]*kH + [hop] + [top_i]*(nH-kH) + [cop_i]
                        
                        EXPAND_OP_PRODUCT({LABEL: lag_label,
                                           OP_RES: res,
                                           OPERATORS: op_list,
                                           IDX_SV: range(1, nop_L + nop_H + 1),
                                           AVOID: set_avoid_list(nop_L, nop_H),
                                           NEW: False,
                                           FIX_VTX: True,
                                           FAC: -fac_L*fac_H})

SUM_TERMS({LABEL_IN:'F_MRCC_LAG',LABEL_RES:'F_MRCC_LAG'})

# dbg
#PRINT_FORMULA({LABEL:'F_MRCC_LAG'})
#ABORT({})
# dbg

export_targets();
