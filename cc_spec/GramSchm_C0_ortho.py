#
# Set targets for projref: reference of nth state will be relaxed with
# the condition that it remains orthogonal to the first (n-1)th reference
# wave function
#
# Yuri, Mar 2015
#
from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

# Set some parameters fom input
#
multistate = inp.get('method.MR.multistate') == 'T'
ciroot = inp.get('method.MR.ciroot')
if (ciroot == None):
    ciroot = 1
else:
    ciroot = int(ciroot)

n_states = ciroot


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


    # dbg
    #PRINT_FORMULA({LABEL:'F_C0_' + str( i_state) + '_prj'})
    # dbg


new_target( 'FOPT_MS_C0_prj')
depend( 'F_MS_C0_prj', 'DEF_ME_C0')

for i_state in range( 2, n_states + 1):
    OPTIMIZE( {LABEL_OPT: 'FOPT_C0_' + str( i_state) + '_prj',
               LABELS_IN: 'F_C0_' + str( i_state) + '_prj'})

export_targets();
