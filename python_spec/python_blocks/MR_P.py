
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk


##
#@date 10.11.15
#@author arne bargholz

import operators.cluster_residue_op 
import operators.ref_wf_op
import operators.H_0

import CASCI.spin_proj
import CASCI.MR


modify_target("do_all")
depend("MakeRefState")
