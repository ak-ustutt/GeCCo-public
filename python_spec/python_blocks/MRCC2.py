
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk


##
#@date 10.11.15
#@author arne bargholz





#------------------------------------------------------------------#
#------------------------------------------------------------------#
#global switches
#------------------------------------------------------------------#
#------------------------------------------------------------------#

ntest=000




#----------------------------------------------------------------#
#Things to do before calculation starts
#---------------------------------------------------------------##
new_target("general considerations",True)





import operators.cluster_residue_op 
import operators.ref_wf_op
import operators.H_0
#Important note: prepared for MR_P preparation of wavefunction
import icMR_orthogonalization.seq_orthogonalization
import preconditioner.diag_precon
import lagrangians.MRCC2lag
import solve.MRCC2solve





#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#
if (not keywords.is_keyword_set("method.MRCC2.excite")):
    modify_target('do_all')
    depend('SOLVE_MRCCPT2')
if (keywords.is_keyword_set("method.relaxref")):
    modify_target('do_all')
    depend('SOLVE_MRCC2ref')
    
