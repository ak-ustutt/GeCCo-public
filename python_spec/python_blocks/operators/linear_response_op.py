from python_interface.gecco_interface import *


new_target('DEF_RESPONSE_OPs')
depend('DEF_T2g', 'DEF_T1', 'DEF_O2g', 'DEF_O1', 'DEF_LAM2g', 'DEF_LAM1')
depend('RefState-Operators')

_op_list={'R1_prime_q':'T1',
          'R1_q':'T1',
          'R2g_prime_q':'T2g',
          'R2g_q':'T2g',
          'R_mu':'C0',
          'AR1_rspns_q':'O1',
          'AR2g_rspns_q':'O2g',
          'AR_rspns_mu':'H_C0',
          'SR1_rspns_q':'O1',
          'SR2g_rspns_q':'O2g',
          'SR_rspns_mu':'C0'}




for _op in _op_list:
    CLONE_OPERATOR({LABEL:_op,
                    TEMPLATE:_op_list[_op]})
DEF_SCALAR({LABEL:'DUMMY'})
    
DEF_SCALAR({LABEL:'Scal_S1'})
DEF_SCALAR({LABEL:'Scal_S2'})
DEF_SCALAR({LABEL:'RED_LAG'})


new_target("MAKE_E_MRCC2")

DEF_SCALAR({LABEL:'E_MRCC2'})
DEF_ME_LIST
