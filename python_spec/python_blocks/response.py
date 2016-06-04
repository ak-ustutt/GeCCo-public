from gecco_interface import *
from gecco_modules.NoticeUtil import*

from gecco_modules.omg_generator import *
import gecco_modules.default_keywords as dk

import MRCC2

import operators.linear_response_op
import lagrangians.MRCCPT2_LRlag
import preconditioner.precon_c0
import solve.MRCC2_LR_solve as solve


modify_target('do all')
for target in solve._solve_eqn_arr:
    depend(target)
