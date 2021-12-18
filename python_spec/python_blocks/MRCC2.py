
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





#----------------------------------------------------------------#
#Things to do before calculation starts
#---------------------------------------------------------------##
new_target("general considerations",True)





import python_blocks.operators.cluster_residue_op 
import python_blocks.operators.ref_wf_op
import python_blocks.operators.H_0
#Important note: prepared for MR_P preparation of wavefunction
import python_blocks.icMR_orthogonalization.seq_orthogonalization
import python_blocks.preconditioner.diag_precon
import python_blocks.lagrangians.MRCC2lag
import python_blocks.solve.MRCC2solve





#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#
if (not keywords.is_keyword_set("method.MRCC2.excite")):
    modify_target('do_all')
    depend('MAKE_MRCC2')
