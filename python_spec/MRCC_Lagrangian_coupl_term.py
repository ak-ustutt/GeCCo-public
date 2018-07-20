#
# Set the coupling term for the MRCC lagrangian:
#
# MS-icMRCC or (ic)SU-MRCC:
#  - sum_i sum_(j != i) 
# <0| C0_i^+ L_i exp(-T_i) exp(T_j) C0_j |0> <0| C0_j^+ exp(-T_i) H exp(-T_i) C0_i |0>
#
# (ic)Mk-MRCC:
#  sum_i sum_(j != i)
# <0| C0_i^+ L_i exp(-T_i) exp(T_j) C0_i |0> <0| C0_i^+ exp(-T_j) H exp(-T_j) C0_j |0> C_j/C_i
#
# (ic)sr-MRCC:
# (ic)SU-MRCC + (ic)Mk-MRCC
#
# for (ic)BW-MRCC see BW_MRCC_Lagrangian.py
#
# i and j refer to state indices.
# The uncoupled part of the Lagrangian is just the summation
# of the tradidional Lagrangian for each state, see set_mrcc_lagrangian.f
#
# The truncation of the BCH expansion is controled by the following input keywords:
# method.MRCC.maxcom_res_cpl and method.MRCC.maxcom_en_cpl
# The first for the L- and the second for the H-part. If not present, by
# method.MRCC.maxcom_res and method.MRCC.maxcom_en
#
#
# Yuri, Nov 2014
#       Ago 2015 -> include Mk and sr MRCC methods
#

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *

inp = GeCCo_Input()
orb = Orb_Info()

from python_interface.gecco_modules.BCH_fac import set_BCH_factor

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

MRCC_type = inp.get('method.MRCC.type')
if (MRCC_type == None):
    MRCC_type = "SU"
if (MRCC_type == "BW"):
    quit_error( 'BW-MRCC has its own Lagrangian file: BW_MRCC_Lagrangian.py')
if (MRCC_type != "SU" and MRCC_type != "Mk" and MRCC_type != "sr"):
    quit_error( 'Unknown MRCC type.')

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
cbar = 'C0_bar'

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

#-----------------------
# The next two targets will be constructed simultaneously,
# since Heff is in intermediate from the Lagrangian:

#  (a) the lagrangian coupling term:
new_target('F_MRCC_LAG_coupl')

#  (b) the Heff_ij, to be used as intermediate:
new_target('F_MS_Heff_int')

# Pay attention and always use
# the modify_target function!
#-----------------------

modify_target( 'F_MRCC_LAG_coupl')
depend( lag_label)
depend( 'C_MS', 'C0_non_orth')
# insist on creating effective H first!
depend('F_MS_Heff_int')

modify_target('F_MS_Heff_int')
depend('E(MR)', 'C0_non_orth')

for i_state in range(1, n_states+1):
    state_label = "" if (i_state == 1) else "_" + str(i_state)
    cdg_i = cop + state_label + '^+'
    lop_i = lop + state_label
    top_i = top + state_label
    cop_i = cop + state_label
    cbardg_i = cbar + state_label + '^+'

    for j_state in range(1, n_states+1):
        if (j_state == i_state):
            continue

        state_label = "" if (j_state == 1) else "_" + str(j_state)
        cdg_j = cop + state_label + '^+'
        top_j = top + state_label
        cop_j = cop + state_label
        cbardg_j = cbar + state_label + '^+'

        MS_coef_ji = 'C_MS_' + str(j_state) + "_" + str(i_state)

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
        # first do Heff, then the coupling term
        modify_target( 'F_MS_Heff_int')
        nop_L = 3 
        for nH in range(0, maxcom_en+1):
            nop_H = 3 + nH

            for kL in range(0, 1):
               # Since pure_vv=F, there is no T on the left side of H. Otherwise use range( 0, nH+1)
               for kH in range(0, 1):
                   fac_H = set_BCH_factor(nH, kH)

                   if (MRCC_type == "SU" or MRCC_type == "sr"):

                       op_list_Heff_int = [cbardg_j] + [top_i]*kH + [hop] + [top_i]*(nH-kH) + [cop_i]
                   if (MRCC_type == 'Mk' or MRCC_type == "sr"):

                       op_list_Heff_int = [cbardg_i] + [top_j]*kH + [hop] + [top_j]*(nH-kH) + [cop_j]

                   EXPAND_OP_PRODUCT({LABEL: 'F_' + heff_ij_op,
                                      OP_RES: heff_ij_op,
                                      OPERATORS: op_list_Heff_int,
                                      IDX_SV: range(1, nop_H + 1),
                                      NEW: nH==0,
                                      FIX_VTX: True,
                                      FAC: fac_H})

        SUM_TERMS({LABEL_IN:'F_' + heff_ij_op, LABEL_RES:'F_' + heff_ij_op})

        # dbg
        #PRINT_FORMULA({LABEL:'F_' + heff_ij_op, MODE:'SHORT'})
        #ABORT({})
        # dbg

        # and now the coupling term
        modify_target( 'F_MRCC_LAG_coupl')
        for nL in range(0, maxcom+1):
            nop_L = 3 + nL
            nop_H = 1

            for kL in range(0, nL+1):
               fac_L = set_BCH_factor(nL, kL)
               # Since pure_vv=F, there is no T on the left side of H. Otherwise use range( 0, nH+1)

               if (MRCC_type == "SU" or MRCC_type == "sr"):

                   op_list = [cdg_i, lop_i] + [top_i]*kL + [top_j]*(nL-kL) + [cop_j] +\
                             [heff_ij_op] 

                   EXPAND_OP_PRODUCT({LABEL: lag_label,
                                      OP_RES: res,
                                      OPERATORS: op_list,
                                      IDX_SV: range(1, nop_L + nop_H + 1),
                                      AVOID: set_avoid_list(nop_L, nop_H),
                                      NEW: False,
                                      FIX_VTX: True,
                                      FAC: -fac_L})

               if (MRCC_type == 'Mk' or MRCC_type == "sr"):

                   op_list = [cdg_i, lop_i] + [top_i]*kL + [top_j]*(nL-kL) + [cop_i] +\
                                      [heff_ij_op] + [MS_coef_ji]

                   EXPAND_OP_PRODUCT({LABEL: lag_label,
                                      OP_RES: res,
                                      OPERATORS: op_list,
                                      IDX_SV: range(1, nop_L + nop_H + 1 + 1),
                                      AVOID: set_avoid_list(nop_L, nop_H),
                                      NEW: False,
                                      FIX_VTX: True,
                                      FAC: fac_L})

SUM_TERMS({LABEL_IN: lag_label, LABEL_RES: lag_label})

# dbg
#PRINT_FORMULA({LABEL: lag_label, MODE:'SHORT'})
#ABORT({})
# dbg

export_targets();
