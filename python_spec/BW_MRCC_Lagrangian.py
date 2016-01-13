#
# Set the BW-MRCC Lagrangian:
#
# sum_i <0| C0_i^+ L_i H exp(T_i) C0_i |0> - E <0| C0_i^+ L_i exp(T_i) C0_i |0>
#
# i refer to state indice.
#
# The truncation of the BCH expansion is controled by the following input keyword:
# method.MRCC.maxcom_res_cpl
# If not present, by
# method.MRCC.maxcom_res
#
# Yuri, Ago 2015
#

import sys,os
sys.path=sys.path+[os.getenv("GECCO_DIR")+"/python_interface"]
from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

from gecco_modules.BCH_fac import set_BCH_factor

# Set some parameters fom input
#
multistate = inp.get('method.MR.multistate') == 'T'
ciroot = inp.get('method.MR.ciroot')
if (ciroot == None):
    ciroot = 1
else:
    ciroot = int(ciroot)

n_states = ciroot

# maximum level in commutator expansion
maxcom_default = 2
maxcom = inp.get('method.MRCC.maxcom_res')
if (maxcom == None):
    maxcom = maxcom_defaut
maxcom = int( maxcom)

maxcom_cpl = inp.get('method.MRCC.maxcom_res_cpl')
if (maxcom_cpl == None):
    maxcom_cpl = maxcom
maxcom_cpl = int( maxcom_cpl)


# Set some variables for operator labels
lag_label = 'F_BW_MRCC_LAG'
res = 'NORM'
top = 'T'
hop = 'H'
lop = 'L'
cop = 'C0'
eop = 'E_MS'

new_target( 'F_BW_MRCC_LAG')
depend( res, top, hop, lop, cop, eop)

for i_state in range(1, n_states+1):
    state_label = "" if (i_state == 1) else "_" + str(i_state)
    cdg_i = cop + state_label + '^+'
    lop_i = lop + state_label
    top_i = top + state_label
    cop_i = cop + state_label

    # Part with Hamiltonian
    for nL in range(0, maxcom+1):
        fac = set_BCH_factor(nL, 0)
        nop = 4 + nL

        op_list = [cdg_i, lop_i, hop] + [top_i]*nL + [cop_i]

        EXPAND_OP_PRODUCT({LABEL: lag_label,
                           OP_RES: res,
                           OPERATORS: op_list,
                           IDX_SV: range(1, nop + 1),
                           NEW: i_state==1 and nL==0,
                           FIX_VTX: True,
                           FAC: fac})

    # Part with energy
    for nL in range(0, maxcom_cpl+1):
        fac = set_BCH_factor(nL, 0)
        nop = 4 + nL

        op_list = [eop, cdg_i, lop_i] + [top_i]*nL + [cop_i]

        EXPAND_OP_PRODUCT({LABEL: lag_label,
                           OP_RES: res,
                           OPERATORS: op_list,
                           IDX_SV: range(1, nop + 1),
                           NEW: False,
                           FIX_VTX: True,
                           FAC: -fac})


SUM_TERMS({LABEL_IN:'F_BW_MRCC_LAG', LABEL_RES:'F_BW_MRCC_LAG'})

# dbg
#PRINT_FORMULA({LABEL:'F_BW_MRCC_LAG'})
#ABORT({})
# dbg

export_targets();
