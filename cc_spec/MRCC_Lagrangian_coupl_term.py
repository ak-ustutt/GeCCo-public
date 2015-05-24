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
#
# Also: set targets for projref: reference of nth state will be relaxed with
# the condition that it remains orthogonal to the first (n-1)th reference
# wave function
#
#
# Yuri, Nov 2014
#
from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

from BCH_fac import set_BCH_factor

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

if inp.get('method.MR.pure_vv') == True:
    quit_error("MRCC_Lagrangian_coupl_term.py: pure_vv=T not implemented for coupling term.")

msc = 1
if (orb.get('ims') != 0):
    msc = 0

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

# The next two targets will be constructed simultaneously:

#  (a) the lagrangian coupling term:
new_target('F_MRCC_LAG_coupl')
depend('F_MRCC_LAG')

#  (b) the Heff_ij, to be used as intermediate:
new_target('F_MS_Heff_int')
depend('E(MR)')

print 'F_MS_Heff_int'

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

        heff_ij_op='MS_Heff_' + str( i_state)  + "_" + str( j_state)

        # Elements of the multistate effective Hamiltonian...
        modify_target( 'F_MS_Heff_int')

        CLONE_OPERATOR({LABEL:heff_ij_op,
                        TEMPLATE:'E(MR)'})

        DEF_ME_LIST({LIST:'ME_' + heff_ij_op,
             OPERATOR:heff_ij_op,
             IRREP:1,
             '2MS':0,
             AB_SYM:msc})

        # and the formulas:
        for nL in range(0, maxcom+1):
            nop_L = 3 + nL
            for nH in range(0, maxcom_en+1):
                nop_H = 3 + nH

                for kL in range(0, nL+1):
                    fac_L = set_BCH_factor(nL, kL)
                    # Since pure_vv=F, there is no T on the left side of H. Otherwise use range( 1, nH+1)
                    for kH in range(0, 1):
                        fac_H = set_BCH_factor(nH, kH)

                        op_list_Heff_int = [cdg_j] + [top_i]*kH + [hop] + [top_i]*(nH-kH) + [cop_i]

                        op_list = [cdg_i, lop_i] + [top_i]*kL  + [top_j]*(nL-kL) + [cop_j] +\
                                  op_list_Heff_int
                        
                        if (nL == 0):
                            modify_target( 'F_MS_Heff_int')
                            EXPAND_OP_PRODUCT({LABEL: 'F_' + heff_ij_op,
                                               OP_RES: heff_ij_op,
                                               OPERATORS: op_list_Heff_int,
                                               IDX_SV: range(1, nop_H + 1),
                                               NEW: nH==0,
                                               FIX_VTX: True,
                                               FAC: fac_H})

                        modify_target( 'F_MRCC_LAG_coupl')
                        EXPAND_OP_PRODUCT({LABEL: lag_label,
                                           OP_RES: res,
                                           OPERATORS: op_list,
                                           IDX_SV: range(1, nop_L + nop_H + 1),
                                           AVOID: set_avoid_list(nop_L, nop_H),
                                           NEW: False,
                                           FIX_VTX: True,
                                           FAC: -fac_L*fac_H})

        modify_target( 'F_MS_Heff_int')
        SUM_TERMS({LABEL_IN:'F_' + heff_ij_op, LABEL_RES:'F_' + heff_ij_op})

        # dbg
        #PRINT_FORMULA({LABEL:'F_' + heff_ij_op})
        #ABORT({})
        # dbg

modify_target( 'F_MRCC_LAG_coupl')
SUM_TERMS({LABEL_IN:'F_MRCC_LAG', LABEL_RES:'F_MRCC_LAG'})

# dbg
#PRINT_FORMULA({LABEL:'F_MRCC_LAG'})
#ABORT({})
# dbg

# Project reference during its optimization: -|C0_1><C0|C0_1> - |C0_2><C0|C0_2> ...
# Special formula for solve_evp in optref=-3
new_target( 'F_MS_C0_prj')
depend( 'C0')

for i_state in range( 2, n_states + 1):
    for i_state2 in range( 1, i_state):
        EXPAND_OP_PRODUCT( {LABEL: 'F_C0_' + str( i_state) + '_prj',
                            OP_RES: 'C0',
                            OPERATORS:['C0', 'C0_'+str( i_state2), 'C0^+', 'C0_'+str( i_state2), 'C0'],
                            IDX_SV:   [   1,                    2,      3,                    4,    1],
                            AVOID: [1,4 , 3,5],
                            FAC: -1.0,
                            NEW: (i_state2 == 1)})

new_target( 'FOPT_MS_C0_prj')
depend( 'F_MS_C0_prj', 'DEF_ME_C0')

for i_state in range( 2, n_states + 1):
    OPTIMIZE( {LABEL_OPT: 'FOPT_C0_' + str( i_state) + '_prj',
               LABELS_IN: 'F_C0_' + str( i_state) + '_prj'})

export_targets();
