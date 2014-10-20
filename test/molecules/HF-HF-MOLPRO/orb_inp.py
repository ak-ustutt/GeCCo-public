from gecco_interface import *

# Get the access for input file and orbital informations by:
inp = GeCCo_Input()
orb = Orb_Info()


#
# Look how to get the informations in the input
#
new_target('INPUT_INFO',True)

PRINT({STRING:'>>> Information from input:'})

shell_names = inp.get( 'orb_space.shell.type')
for i in range( 0, len(shell_names)):
    shell = inp.get( 'orb_space.shell.def')[i]
    PRINT({STRING: shell_names[i] + ' shell: ' + str( shell)})





#
# And how to get informations about the orbitals
#
new_target('ORB_INFO',True)

n_active_el = orb.get( 'nactel')
PRINT({STRING:'Number of active electrons: ' + str( n_active_el)})

export_targets();
