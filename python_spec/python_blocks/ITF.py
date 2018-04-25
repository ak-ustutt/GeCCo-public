"""Simple implementation of the icMRCC method based on code in icMRCC.py by Yuri Aoto

Usage:

To generate and test formulae for ITF translator

method
 MR_P
 ITF en_type=2,res_type=2
calculate
 solve non_linear optref=0,

"""

from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk



#----------------------------------------------------------------#
# Things to do before calculation starts
#---------------------------------------------------------------##
new_target("general considerations",True)

import operators.cluster_residue_op 
import operators.ref_wf_op
import operators.H_0
#Important note: prepared for MR_P preparation of wavefunction
import icMR_orthogonalization.seq_orthogonalization
import preconditioner.diag_precon
import lagrangians.ITFlag
import solve.ITFsolve


#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#
modify_target('do_all')
depend('SOLVE_ITF')
