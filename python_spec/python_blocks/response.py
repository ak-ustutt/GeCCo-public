from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

from python_interface.gecco_modules.omg_generator import *
import python_interface.gecco_modules.default_keywords as dk

import python_blocks.MRCC2

import python_blocks.operators.linear_response_op
import python_blocks.lagrangians.MRCC2_LRlag
import python_blocks.preconditioner.precon_c0
import python_blocks.solve.MRCC2_LR_solve as solve


modify_target('do_all')
for target in solve._solve_eqn_arr:
    depend(target)
