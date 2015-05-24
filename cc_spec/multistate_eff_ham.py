#
# Set targets for multistate effective Hamiltonian
# 
# yuri, nov 2014

from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

from BCH_fac import set_BCH_factor

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


# Multistate packed effective Hamiltonian:
# similar to the E(MR) but with C0_1,
# Heff = < 0| C0^+ e^-T H e^T C0_1 |0>
#
# T must run with C0_1
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
depend('F_MRCC_E', 'pack_Heff_MS')
# This INVARIANT is only to change the OP_RES.
# Its a good idea implement this in the REPLACE rule
INVARIANT({LABEL_RES:'F_pack_Heff_MS',
           LABEL_IN:'F_MRCC_E_1',
           OP_RES:'pack_Heff_MS',
           OPERATORS:['D']})

REPLACE({LABEL_RES:'F_pack_Heff_MS',
         LABEL_IN:'F_pack_Heff_MS',
         OP_LIST:['C0','C0_1']})
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
    ADV_STATE({LISTS:'ME_C0',
               N_ROOTS:n_states})
    if (i_state%n_states == 0):
        ADV_STATE({OPERATORS:'C0_1',
                   N_ROOTS:n_states,
                   USE1:True})
        ADV_STATE({OPERATORS:'T',
                   N_ROOTS:n_states})




# Multistate packed coupling-states Hamiltonian.
# These are the matrix elements of the Hamiltonian
# in the e^T C0 |0> basis:
#
# Hcpl = < 0| C0^+ e^T_2^+ H e^T C0_1 |0>
#
# T   must run with C0_1
# T_2 must run with C0
#
# The elemetns of Hcpl are evaluated with
# the same formula but advancing the records

# The operator and its ME-list
new_target('pack_Hcpl_MS')
depend('E(MR)')
CLONE_OPERATOR({LABEL:'pack_Hcpl_MS',TEMPLATE:'E(MR)'})

DEF_ME_LIST({LIST:'ME_pack_Hcpl_MS',
             OPERATOR:'pack_Hcpl_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})

# Its formula
new_target('F_pack_Hcpl_MS')
depend('C0', 'T', 'H', 'pack_Hcpl_MS')

for n in range(0, maxcom_en+1):
    for k in range(0, n+1):
        EXPAND_OP_PRODUCT({LABEL: 'F_pack_Hcpl_MS',
                           OP_RES: 'pack_Hcpl_MS',
                           OPERATORS: ['C0^+'] + ['T_2^+']*k + ['H'] + ['T']*(n-k) + ['C0_1'],
                           IDX_SV: range(1, n+3 + 1),
                           NEW: n==0,
                           FIX_VTX: True,
                           FAC: set_BCH_factor(n, k, False)})

SUM_TERMS({LABEL_IN: 'F_pack_Hcpl_MS', LABEL_RES: 'F_pack_Hcpl_MS'})

# dbg
#PRINT_FORMULA({LABEL:'F_pack_Hcpl_MS'})
#ABORT({})
# dbg end

# Optimize it
new_target('FOPT_pack_Hcpl_MS')
depend('F_pack_Hcpl_MS')
OPTIMIZE({LABEL_OPT:'FOPT_pack_Hcpl_MS',
          LABELS_IN:'F_pack_Hcpl_MS'})

# and evaluate it!
new_target('EVAL_pack_Hcpl_MS')
depend('FOPT_pack_Hcpl_MS', 'SOLVE_MRCC')

SET_STATE({OPERATORS:'T_2',
           ISTATE:1})

for i_state in range( 1, n_states*n_states+1):
    EVALUATE({FORM:'FOPT_pack_Hcpl_MS'})

    ADV_STATE({LISTS:['ME_pack_Hcpl_MS'],
               N_ROOTS:n_states*n_states})
    ADV_STATE({LISTS:'ME_C0',
               N_ROOTS:n_states})
    ADV_STATE({OPERATORS:'T_2',
               N_ROOTS:n_states})
    if (i_state%n_states == 0):
        ADV_STATE({OPERATORS:'C0_1',
                   N_ROOTS:n_states,
                   USE1:True})
        ADV_STATE({OPERATORS:'T',
                   N_ROOTS:n_states})

# Back to the right state
SET_STATE({OPERATORS:'T_2',
           ISTATE:2})





# Multistate packed overlap matrix:
# Smat = < 0| C0^+ e^T_2^+ e^T C0_1 |0>
#
# T   must run with C0_1
# T_2 must run with C0
#
# The elemetns of Smat are evaluated with
# the same formula but advancing the records

# The operator and its ME-list
new_target('pack_Smat_MS')
depend('E(MR)')
CLONE_OPERATOR({LABEL:'pack_Smat_MS',TEMPLATE:'E(MR)'})

DEF_ME_LIST({LIST:'ME_pack_Smat_MS',
             OPERATOR:'pack_Smat_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})

# and its formula
new_target('F_pack_Smat_MS')
depend('C0', 'T', 'H', 'pack_Smat_MS')

for n in range(0, maxcom_en+1):
    for k in range(0, n+1):
        EXPAND_OP_PRODUCT({LABEL: 'F_pack_Smat_MS',
                           OP_RES: 'pack_Smat_MS',
                           OPERATORS: ['C0^+'] + ['T_2^+']*k + ['T']*(n-k) + ['C0_1'],
                           IDX_SV: range(1, n+2 + 1),
                           NEW: n==0,
                           FIX_VTX: True,
                           FAC: set_BCH_factor(n, k, False)})

SUM_TERMS({LABEL_IN: 'F_pack_Smat_MS', LABEL_RES: 'F_pack_Smat_MS'})

# dbg
#PRINT_FORMULA({LABEL:'F_pack_Smat_MS'})
#ABORT({})
# dbg end

# Optimize it
new_target('FOPT_pack_Smat_MS')
depend('F_pack_Smat_MS')
OPTIMIZE({LABEL_OPT:'FOPT_pack_Smat_MS',
          LABELS_IN:'F_pack_Smat_MS'})

# and evaluate it!
new_target('EVAL_pack_Smat_MS')
depend('FOPT_pack_Smat_MS', 'SOLVE_MRCC')

SET_STATE({OPERATORS:'T_2',
           ISTATE:1})

for i_state in range( 1, n_states*n_states+1):
    EVALUATE({FORM:'FOPT_pack_Smat_MS'})

    ADV_STATE({LISTS:['ME_pack_Smat_MS'],
               N_ROOTS:n_states*n_states})
    ADV_STATE({LISTS:'ME_C0',
               N_ROOTS:n_states})
    ADV_STATE({OPERATORS:'T_2',
               N_ROOTS:n_states})
    if (i_state%n_states == 0):
        ADV_STATE({OPERATORS:'C0_1',
                   N_ROOTS:n_states,
                   USE1:True})
        ADV_STATE({OPERATORS:'T',
                   N_ROOTS:n_states})

# Back to the right state
SET_STATE({OPERATORS:'T_2',
           ISTATE:2})



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


# Diagonalise the effective Hamiltonian
new_target( 'SOLVE_Heff_MS', multistate and not(skip))
depend( 'EVAL_pack_Heff_MS', 'C_MS', 'E_MS')

EVP_PACKED_OP({LIST_IN:'ME_pack_Heff_MS',
               LIST_E:'ME_E_MS',
               LIST_EVEC:'ME_C_MS',
               N_ROOTS:n_states})



# Diagonalise the coupling state Hamiltonian
new_target( 'SOLVE_Hcpl_MS', multistate and not(skip))
depend( 'EVAL_pack_Hcpl_MS', 'EVAL_pack_Smat_MS', 'C_MS', 'E_MS')

EVP_PACKED_OP({LIST_IN:'ME_pack_Hcpl_MS',
               LIST_S:'ME_pack_Smat_MS',
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
