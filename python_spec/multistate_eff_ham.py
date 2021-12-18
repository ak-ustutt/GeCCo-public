#
# Set targets for multistate effective Hamiltonian
# 
# yuri, nov 2014
#       sep 2015 -> extend to non-orthonormal basis in the P-space
#

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *

inp = GeCCo_Input()
orb = Orb_Info()

from python_interface.gecco_modules.BCH_fac import set_BCH_factor

multistate=inp.get('method.MR.multistate') == 'T'
coupled_states=inp.get('method.MRCC.coupled_states')
coupled_states=coupled_states=='T' or coupled_states==None
ciroot=int(inp.get('method.MR.ciroot'))
skip=inp.is_keyword_set('calculate.skip_E')

msc = 1
if (orb.get('ims') != 0):
    msc = 0

n_states = 1
if (multistate):
    n_states = ciroot
    if (n_states < 2):
        quit_error("Please, do not use multistate=T for one state!")

maxcom_en = inp.get('method.MRCC.maxcom_en')
if (maxcom_en == None):
    maxcom_en = 2
maxcom_en = int(maxcom_en)

# Symmetrize Heff?
Heff_symm = inp.get('method.MRCC.Heff_symm') == 'T'

# Assume orthonormal basis?
assume_orth = inp.get('method.MRCC.assume_orth') == 'T'

#==================================================
# packed overlap matrix in P-space:
# S = < 0| C0^+ C0_1 |0>
#
# Used for MRCC theories with non-orthonormal
# basis in the P-space
#
# The elemetns of S are evaluated with
# the same formula but advancing the state
new_target('pack_P_S')
depend('E(MR)')
CLONE_OPERATOR({LABEL:'pack_P_S',TEMPLATE:'E(MR)'})

DEF_ME_LIST({LIST:'ME_pack_P_S',
             OPERATOR:'pack_P_S',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})

# and its formula
new_target('F_pack_P_S')
depend('C0', 'pack_P_S')

EXPAND_OP_PRODUCT({LABEL: 'F_pack_P_S',
                   OP_RES:  'pack_P_S',
                   OPERATORS: ['C0^+', 'C0_1'],
                   IDX_SV:    [1     , 2],
                   NEW: True,
                   FIX_VTX: True,
                   FAC: 1.0})

# Optimize it
new_target('FOPT_pack_P_S')
depend('F_pack_P_S')
OPTIMIZE({LABEL_OPT:'FOPT_pack_P_S',
          LABELS_IN:   'F_pack_P_S'})

# and evaluate it!
new_target('EVAL_pack_P_S')
depend('FOPT_pack_P_S', 'SOLVE_REF')

SET_STATE({LIST:['ME_C0','ME_pack_P_S'],
           OPERATORS:'C0_1',
           USE1:True,
           ISTATE:1})

for i_state in range( 1, n_states*n_states+1):

    EVALUATE({FORM:'FOPT_pack_P_S'})

# dbg
#    PRINT_MEL({LIST:'ME_pack_P_S',
#               COMMENT:'ME_pack_P_S'})
# dbgend

    ADV_STATE({LISTS:['ME_pack_P_S'],
               N_ROOTS:n_states*n_states})
    ADV_STATE({LISTS:'ME_C0',
               N_ROOTS:n_states})
    if (i_state%n_states == 0):
        ADV_STATE({OPERATORS:'C0_1',
                   N_ROOTS:n_states,
                   USE1:True})

#==================================================
# Inverted overlap matrix
#
# operators, ME-lists and evaluation of inverse of
# the overlap matrix in the P-space
#

# Create operators and ME-lists
new_target( 'P_Sinv')
depend('pack_P_S')

DEF_ME_LIST({LIST:'ME_pack_P_Sinv',
             OPERATOR:'pack_P_S',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})

ASSIGN_ME2OP({LIST:'ME_pack_P_S',
              OPERATOR:'pack_P_S'})

for i_state in range(1, n_states+1):
    state_label = "_" + str(i_state)
    CLONE_OPERATOR({LABEL:'P_Sinv' + state_label,
                    TEMPLATE:'E(MR)'})

    DEF_ME_LIST({LIST:'ME_P_Sinv' + state_label,
                 OPERATOR:'P_Sinv' + state_label,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:msc,
                 MIN_REC:1,
                 MAX_REC:n_states,
                 REC:1})

# Invert and copy to desired ME-lists
new_target( 'invert_P_S')
depend( 'EVAL_pack_P_S', 'P_Sinv')

INV_PACKED_OP({LIST_IN:'ME_pack_P_S',
               LIST_OUT:'ME_pack_P_Sinv',
               N_ROOTS:n_states})

j_state = 1
state_label = '_1'
for i_state in range(1, n_states*n_states+1):

    SCALE_COPY({LIST_INP:'ME_pack_P_Sinv',
                LIST_RES:'ME_P_Sinv' + state_label,
                FAC:1.0})

# dbg
#    PRINT_MEL({LIST:'ME_pack_P_Sinv',
#               COMMENT:'ME_pack_P_Sinv'})
#    PRINT_MEL({LIST:'ME_P_Sinv' + state_label,
#               COMMENT:'ME_P_Sinv' + state_label})
# dbg

    ADV_STATE({LISTS:'ME_pack_P_Sinv',
               N_ROOTS:n_states*n_states})
    ADV_STATE({LISTS:'ME_P_Sinv' + state_label,
               N_ROOTS:n_states})

    if (i_state%n_states == 0):
        j_state += 1
        state_label = '_' + str( j_state)


#==================================================
# The modified reference in non-orthonormal basis
# for the model space:
#
# <C0_bar| = sum_k P_Sinv_k <C0_k|
#
# and the corresponding spread MEl, C0_bar_i,
# that have the property:
#
# <C0_bar_i|C0_j> = delta_ij
#
new_target('C0_non_orth')
depend('C0')
if (not(assume_orth)):
    depend('P_Sinv')

for i_state in range(0, n_states+1):
    state_label = "" if (i_state == 0) else "_" + str(i_state)

    CLONE_OPERATOR({LABEL:'C0_bar' + state_label,
                    TEMPLATE:'C0'})

    DEF_ME_LIST({LIST:'ME_C0_bar' + state_label,
                 OPERATOR:'C0_bar' + state_label,
                 IRREP:orb.get('lsym'),
                 '2MS':orb.get('ims'),
                 AB_SYM:msc,
                 MIN_REC:1,
                 MAX_REC:n_states if (i_state == 0) else 1,
                 REC:1})

if (not(assume_orth)):
    for k_state in range(1, n_states+1):
        state_label = "_" + str(k_state)
        EXPAND_OP_PRODUCT({LABEL: 'F_C0_bar',
                           OP_RES: 'C0_bar',
                           OPERATORS: ['C0_bar','P_Sinv' + state_label, 'C0' + state_label,'C0_bar'],
                           IDX_SV:    [1       , 2                    , 3                 ,1],
                           NEW: k_state == 1,
                           FIX_VTX: True,
                           FAC: 1.0})

    # dbg
    #PRINT_FORMULA({LABEL:'F_C0_bar'})
    #ABORT({})
    # dbg end

    OPTIMIZE({LABEL_OPT:'FOPT_C0_bar',
              LABELS_IN:   'F_C0_bar'})

# Evaluate it!
new_target('EVAL_C0_non_orth')
depend('C0_non_orth')
if (not(assume_orth)):
    depend('invert_P_S')

for i_state in range( 1, n_states+1):

    if (assume_orth):
        SCALE_COPY({LIST_INP:'ME_C0',
                    LIST_RES:'ME_C0_bar',
                    FAC:1.0})
    else:
        EVALUATE({FORM:'FOPT_C0_bar'})

    ADV_STATE({LISTS:['ME_C0_bar'],
               N_ROOTS:n_states})
    if (assume_orth):
        ADV_STATE({LISTS:['ME_C0'],
                   N_ROOTS:n_states})
    else:
        for j_state in range( 1, n_states+1):
            state_label = "_" + str(j_state)
            ADV_STATE({LISTS:['ME_P_Sinv' + state_label],
                       N_ROOTS:n_states})

lists_to_spread = []
for i_state in range( 1, n_states+1):
    lists_to_spread.append( 'ME_C0_bar_' + str(i_state))

SPREAD_MEL({LIST_IN:'ME_C0_bar',
            LIST_OUT:lists_to_spread})

# dbg
#for i_state in range( 1, n_states+1):
#    PRINT_MEL({LIST:'ME_C0_bar',
#               COMMENT:'C0_bar'})
#    PRINT_MEL({LIST:'ME_C0_bar_'+str(i_state),
#               COMMENT:'C0_bar'})
#    ADV_STATE({LISTS:['ME_C0_bar'],
#               N_ROOTS:n_states})
#ABORT({})
# dbgend


#==================================================
# Multistate packed effective Hamiltonian:
# similar to the E(MR) but with C0_1, and C0_bar^+
# Heff = < 0| C0_bar^+ e^-T H e^T C0_1 |0>
#
# T must run with C0_1
# C0_bar = C0 for orthonormal basis
#
# The elements of Heff are evaluated with
# the same formula but advancing the records

# The operator and its ME-list
new_target('pack_Heff_MS')
depend('E(MR)')

CLONE_OPERATOR({LABEL:'pack_Heff_MS',TEMPLATE:'E(MR)'})

DEF_ME_LIST({LIST:'ME_pack_Heff_MS',
             OPERATOR:'pack_Heff_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})

# Its formula
new_target('F_pack_Heff_MS')
depend('F_MRCC_E', 'pack_Heff_MS', 'EVAL_C0_non_orth')

# This INVARIANT is only to change the OP_RES.
# Its a good idea implement this in the REPLACE rule
INVARIANT({LABEL_RES:'F_pack_Heff_MS',
           LABEL_IN:'F_MRCC_E_1',
           OP_RES:'pack_Heff_MS',
           OPERATORS:['D']})

REPLACE({LABEL_RES:'F_pack_Heff_MS',
         LABEL_IN:'F_pack_Heff_MS',
         OP_LIST:['C0','C0_1']})

REPLACE({LABEL_RES:'F_pack_Heff_MS',
         LABEL_IN:'F_pack_Heff_MS',
         OP_LIST:['C0^+','C0_bar^+']})

# dbg
#PRINT_FORMULA({LABEL:'F_pack_Heff_MS'})
#ABORT({})
# dbg end

# Optimize it
new_target('FOPT_pack_Heff_MS')
depend('F_pack_Heff_MS')
OPTIMIZE({LABEL_OPT:'FOPT_pack_Heff_MS',
          LABELS_IN:'F_pack_Heff_MS'})

# and evaluate it!
new_target('EVAL_pack_Heff_MS')
depend('FOPT_pack_Heff_MS', 'SOLVE_MRCC')

for i_state in range( 1, n_states*n_states+1):
    EVALUATE({FORM:'FOPT_pack_Heff_MS'})

    ADV_STATE({LISTS:'ME_pack_Heff_MS',
               N_ROOTS:n_states*n_states})
    ADV_STATE({LISTS:'ME_C0_bar',
               N_ROOTS:n_states})
    if (i_state%n_states == 0):
        ADV_STATE({OPERATORS:'C0_1',
                   N_ROOTS:n_states,
                   USE1:True})
        ADV_STATE({OPERATORS:'T',
                   N_ROOTS:n_states})

if (Heff_symm):
    DEF_ME_LIST({LIST:'ME_symm_buf',
                 OPERATOR:'pack_Heff_MS',
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:msc,
                 MIN_REC:1,
                 MAX_REC:n_states*n_states,
                 REC:1})

    DEF_ME_LIST({LIST:'ME_symm_buf_2',
                 OPERATOR:'pack_Heff_MS',
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:msc,
                 MIN_REC:1,
                 MAX_REC:n_states*n_states,
                 REC:1})

    ASSIGN_ME2OP({LIST:'ME_pack_Heff_MS',
                  OPERATOR:'pack_Heff_MS'})

    for i_state in range( 1, n_states*n_states+1):
        i_tr = n_states*((i_state-1)%n_states) + (i_state-1)/n_states + 1
        SET_STATE({ISTATE:i_state,
                   LISTS:'ME_pack_Heff_MS'})
        SET_STATE({ISTATE:i_tr,
                   LISTS:'ME_symm_buf'})
        SET_STATE({ISTATE:i_state,
                   LISTS:'ME_symm_buf_2'})
        SCALE_COPY({LIST_INP:'ME_pack_Heff_MS',
                    LIST_RES:'ME_symm_buf',
                    FAC:0.5})
        SCALE_COPY({LIST_INP:'ME_pack_Heff_MS',
                    LIST_RES:'ME_symm_buf_2',
                    FAC:0.5})

    for i_state in range( 1, n_states*n_states+1):

        SET_STATE({ISTATE:i_state,
                   LISTS:['ME_pack_Heff_MS','ME_symm_buf','ME_symm_buf_2']})

        ADD({LIST_SUM:'ME_pack_Heff_MS',
             LISTS:['ME_symm_buf','ME_symm_buf_2'],
             FAC:[1.0,1.0]})

        PRINT_MEL({LIST:'ME_pack_Heff_MS',
                   COMMENT:'ME_pack_Heff_MS'})

    DELETE_ME_LIST({LIST:'ME_symm_buf'})
    DELETE_ME_LIST({LIST:'ME_symm_buf_2'})

SET_STATE({ISTATE:1,
           LISTS:'ME_pack_Heff_MS'})

#==================================================
# Multi-state wave function coefficients
new_target('C_MS')
depend('C0')
CLONE_OPERATOR({LABEL:'C_MS',TEMPLATE:'C0'})

DEF_ME_LIST({LIST:'ME_C_MS',
             OPERATOR:'C_MS',
             IRREP:orb.get('lsym'),
             '2MS':orb.get('ims'),
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states,
             REC:1})

for i_state in range( 1, n_states+1):
    for j_state in range( 1, n_states+1):

        ij = '_' + str(i_state) + '_' + str(j_state)

        DEF_HAMILTONIAN({LABEL:'C_MS'+ij,
                         MIN_RANK:0,
                         MAX_RANK:0})

        DEF_ME_LIST({LIST: 'ME_C_MS' + ij,
                     OPERATOR:'C_MS' + ij,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:0,
                     MIN_REC:1,
                     MAX_REC:n_states,
                     REC:1})

        for k_state in range( 1, n_states+1):

            SET_MEL({LIST: 'ME_C_MS' + ij,
                     IDX_LIST: [1],
                     VAL_LIST: [1.0]})

            ADV_STATE({LISTS:'ME_C_MS' + ij,
                       N_ROOTS:n_states})

# Multi-state energies
new_target('E_MS')
DEF_HAMILTONIAN({LABEL:'E_MS',
                 MIN_RANK:0,
                 MAX_RANK:0})

DEF_ME_LIST({LIST:'ME_E_MS',
             OPERATOR:'E_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:0,
             MIN_REC:1,
             MAX_REC:n_states,
             REC:1})

#==================================================
# Diagonalise the effective Hamiltonian
new_target( 'SOLVE_Heff_MS', multistate and not(skip))
depend( 'EVAL_pack_Heff_MS', 'C_MS', 'E_MS')

EVP_PACKED_OP({LIST_IN:'ME_pack_Heff_MS',
               LIST_E:'ME_E_MS',
               LIST_EVEC:'ME_C_MS',
               N_ROOTS:n_states})

# dbg
#for i_state in range(1, n_states+1):
#    PRINT_MEL({LIST:'ME_E_MS',
#               COMMENT:'Multistate energy for state ' + str(i_state)})
#    PRINT_MEL({LIST:'ME_C_MS',
#               COMMENT:'Multistate CI coefficients for state ' + str(i_state)})
#
#    ADV_STATE({LISTS:['ME_E_MS','ME_C_MS'],
#               N_ROOTS:n_states})
# dbg end

export_targets();
