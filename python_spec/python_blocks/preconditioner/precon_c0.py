from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

new_target("DEF_PRECON_C0")
depend("MakeRefState")

CLONE_OPERATOR({LABEL:"PRECONC0",
                TEMPLATE:"C0"
            })
