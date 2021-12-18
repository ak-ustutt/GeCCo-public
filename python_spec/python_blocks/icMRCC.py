"""An implementation of the icMRCCSD theory

Usage:

For a calculation with fixed reference, the GeGGo input needs (possibly among other things):

method
 MR_P
 MRCC_new
calculate
 solve non_linear optref=0,

For a calculation relaxing the reference:

method
 MR project=4,prc_traf=F
calculate
 solve non_linear optref=-3
 skip_E


TODO:


For details, see:
  solve/icMRCCSDsolve.py
  lagrangians/icMRCCSDlag.py

History:

Yuri Aoto, august 2017: Creation based on MRCC2.py

"""

from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk



#----------------------------------------------------------------#
# Things to do before calculation starts
#---------------------------------------------------------------##
new_target("general considerations",True)

import python_blocks.operators.cluster_residue_op
import python_blocks.operators.ref_wf_op
import python_blocks.operators.H_0
#Important note: prepared for MR_P preparation of wavefunction
import python_blocks.icMR_orthogonalization.seq_orthogonalization
import python_blocks.preconditioner.diag_precon
import python_blocks.lagrangians.icMRCCSDlag
import python_blocks.solve.icMRCCSDsolve


#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#
modify_target('do_all')
depend('SOLVE_MRCC')

