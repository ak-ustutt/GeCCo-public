#
# Set targets for multistate effective Hamiltonian
# 
# yuri, nov 2014

from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

multistate=inp.get('method.MR.multistate') == 'T'
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


# Packed Effective Hamiltonian for multistate...
new_target('pack_Heff_MS')
depend('E(MR)')
CLONE_OPERATOR({LABEL:'pack_Heff_MS',TEMPLATE:'E(MR)'})

# and its ME list
new_target('DEF_ME_pack_Heff_MS')
depend('pack_Heff_MS')
DEF_ME_LIST({LIST:'ME_pack_Heff_MS',
             OPERATOR:'pack_Heff_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states*n_states,
             REC:1})


# Multi-state wave function coefficients...
new_target('C_MS')
depend('C0')
CLONE_OPERATOR({LABEL:'C_MS',TEMPLATE:'C0'})

# and its ME list
new_target('DEF_ME_C_MS')
depend('C_MS')
DEF_ME_LIST({LIST:'ME_C_MS',
             OPERATOR:'C_MS',
             IRREP:orb.get('lsym'),
             '2MS':orb.get('ims'),
             AB_SYM:msc,
             MIN_REC:1,
             MAX_REC:n_states,
             REC:1})


# Multi-state energies ...
new_target('E_MS')
DEF_HAMILTONIAN({LABEL:'E_MS',
                 MIN_RANK:0,
                 MAX_RANK:0})

# and its ME list
new_target('DEF_ME_E_MS')
depend('E_MS')
DEF_ME_LIST({LIST:'ME_E_MS',
             OPERATOR:'E_MS',
             IRREP:1,
             '2MS':0,
             AB_SYM:0,
             MIN_REC:1,
             MAX_REC:n_states,
             REC:1})


# Multistate packed effective Hamiltonian:
# similar to the E(MR) but with C0_1,
# < 0| C0^+ e^-T H e^T C0_1 |0>

# Set its formula...
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

# optimize it...
new_target('FOPT_pack_Heff_MS')
depend('F_pack_Heff_MS','DEF_ME_pack_Heff_MS')
OPTIMIZE({LABEL_OPT:'FOPT_pack_Heff_MS',
          LABELS_IN:'F_pack_Heff_MS'})

# and evaluate:
new_target('EVAL_pack_Heff_MS')
depend('FOPT_pack_Heff_MS','SOLVE_MRCC')
for i_state in range( 1, n_states*n_states+1):
    EVALUATE({FORM:'FOPT_pack_Heff_MS'})
    # dbg
    #PRINT_MEL({LIST:'ME_pack_Heff_MS',
    #           COMMENT:'packed H_eff (' + str(i_state) + '):'})
    # dbg end

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


# Diagonalise the effective Hamiltonian
new_target('SOLVE_Heff_MS',multistate and not(skip))
depend('EVAL_pack_Heff_MS','DEF_ME_C_MS','DEF_ME_E_MS')

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
