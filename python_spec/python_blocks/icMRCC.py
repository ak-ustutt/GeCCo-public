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

# perturbative correction requested?
word = keywords.get('method.MR.pertCorr')
if word is None:
    pertCorr =  False
else:
    if word == "F":
        pertCorr = False
    elif word == "T":
        pertCorr = True
    else:
        quit_error('pertCorr must be T or F, found: '+word)

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
import python_blocks.lagrangians.icMRCCSDlag  # includes some preparations of (T) if requested
import python_blocks.solve.icMRCCSDsolve
if pertCorr:
    import python_blocks.solve.icMRCCpTsolve


#-----------------------------------------------------------------#
# ... do it ...
#-----------------------------------------------------------------#
modify_target('do_all')
depend('SOLVE_MRCC')
if pertCorr:
    depend('SOLVE_MRCC_PT')

