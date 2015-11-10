print("python start"+"-"*20)

from gecco_interface import *
from python_blocks.Arnes_helfer import *

from gecco_modules.omg_generator import *
import gecco_modules.default_keywords as dk


##
#@date 5.6.15
#@author arne bargholz





####################################################################
####################################################################
#global switches
####################################################################
####################################################################



##################################################################
#Things to do before calculation starts
##################################################################
new_target("general considerations",True)

spinadapt=0
if keywords.is_keyword_set('calculate.routes.spinadapt'):
    spinadapt=int(keywords.get('calculate.routes.spinadapt'))
    notice("spinadapt",spinadapt)
  # ms value of wavefunction



imult =int(orbitals.get('imult'))     



ims = int(orbitals.get('ims'))       
imult =int(orbitals.get('imult'))

if (ims == 0) and ((imult-1)%4 == 0) :
    msc = 1
elif (ims == 0) and ((imult+1)%4 == 0) :
    msc = -1
else :
    msc = 0                        # msc is the AB_sym for many operators

wf_sym=orbitals.get('lsym')        #symmetry of wavefunction


if orbitals.get('nactel')<2:
    ABORT({STRING:'I am expecting at least 2 active electrons'})

###################################################################
###################################################################
# calculate the reference state to get C0
###################################################################
###################################################################
import operators.cluster_residue_op 
import operators.ref_wf_op
import operators.H_0
import CASCI.spin_proj
import CASCI.MR
import icMR_orthogonalization.seq_orthogonalization
import preconditioner.split_precon
import lagrangians.MRCC2lag
import solve.MRCC2solve







#######################################################################################
#Define Operators for PT
######################################################################################


###################################################################
# ... do it ...
###################################################################

new_target('do-all',True)
#depend('MakeRefState','MakeOrthBasis')
depend('SOLVE_MRCCPT2')


#debug_MEL('X_TRM_LIST')
#debug_MEL('ME_Dtr')

#debug_MEL('ME_A')

#debug_MEL('A_TRF_LST')


#debug_FORM('F_Atr')
#debug_FORM('FORM_A_TRF_FINAL')

#debug_MEL('ME_Dinv')

#debug_MEL('ME_D')
#debug_MEL('GSLST')

