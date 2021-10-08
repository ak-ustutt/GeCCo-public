
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk


##
#@date 24.11.15
#@author arne bargholz


####################################################################
####################################################################
#global switches
####################################################################
####################################################################

ntest=000

##################################################################
#Things to do before calculation starts
##################################################################
new_target("general considerations",True)

ims = int(orbitals.get('ims'))   
notice("ims",ims)    
imult =int(orbitals.get('imult'))
notice("imult",imult)

if (ims == 0) and ((imult-1)%4 == 0) :
    msc = 1
elif (ims == 0) and ((imult+1)%4 == 0) :
    msc = -1
else :
    msc = 0                        # msc is the AB_sym for many operators

wf_sym=orbitals.get('lsym')        #symmetry of wavefunction


import python_blocks.operators.cluster_residue_op 
import python_blocks.operators.ref_wf_op
import python_blocks.operators.H_0
#Important note: prepared for MR_P preparation of wavefunction
import python_blocks.icMR_orthogonalization.seq_orthogonalization
import python_blocks.preconditioner.diag_precon
import python_blocks.lagrangians.MRCCPT2lag
import python_blocks.solve.MRCCPT2solve

#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#

modify_target('do_all')
depend('SOLVE_MRCCPT2')
