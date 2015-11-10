
from gecco_interface import *
from gecco_modules.NoticeUtil import*

from gecco_modules.omg_generator import *
import gecco_modules.default_keywords as dk


##
#@date 10.11.15
#@author arne bargholz





####################################################################
####################################################################
#global switches
####################################################################
####################################################################

ntest=000

#creating an objects to access orbital information and input information




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
#Important note: prepared for MR_P preparation of wavefunction
import icMR_orthogonalization.seq_orthogonalization
import preconditioner.diag_precon
import lagrangians.MRCC2lag
import solve.MRCC2solve






###################################################################
###################################################################
#setting up the lagrangian (from external file)
###################################################################
###################################################################

#######################################################################################
#Define Operators for PT
######################################################################################


###################################################################
# ... do it ...
###################################################################

modify_target('do all')
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

