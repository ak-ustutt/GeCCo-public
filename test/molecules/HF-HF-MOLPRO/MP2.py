#
# The MP2 tutorial as a tutorial for the python interface to GeCCo
#
# First one has to import the interface module:

import sys,os
sys.path=sys.path+[os.getenv("GECCO_DIR")+"/python_interface"]
from gecco_interface import *

# One can define targets, like this:
new_target('MP2_HAM')

# And rules like this:
rule('DEF_HAMILTONIAN',
     {'LABEL':'FOCK_2',
      'MIN_RANK':1,
      'MAX_RANK':1})

# or this
DEF_HAMILTONIAN({LABEL:'FOCK',
                 MIN_RANK:1,
                 MAX_RANK:1})
#
# The first always works. Is the "real" way. The second is a wrapper
# to make the code more appealing, but it depends if the rules and
# arguments names are writen in the appropriate files
#


#
# This is a rule of the target just defined. The first argument is 
# the name of the rule and the second is a dictionary: a fundamental
# data structure in python (a hash), that contains pairs <key>:<value>
# It is just perfect for our needs
#

new_target('MP2_OPS')
#
# Now we defined a new target and the following rules
# are rules of this target.
#
DEF_SCALAR({LABEL:'LMP2'})
DEF_SCALAR({LABEL:'EMP2'})
DEF_EXCITATION({LABEL:'T2',MIN_RANK:2,MAX_RANK:2})
CLONE_OPERATOR({LABEL:'O2',TEMPLATE:'T2'})


# If for some reason we need to add a new rule for an old
# target, we can do it like this:
modify_target('MP2_HAM')
# and the rules now are referet to this target.
#
DEF_HAMILTONIAN({LABEL:'PHI',
                 MIN_RANK:2,
                 MAX_RANK:2})


new_target('MP2_FORM')
# This is how dependencies are added:
depend('MP2_HAM','MP2_OPS')
# The following rule has lists. As all python lists,
# they are defined with square brackets:
EXPAND_OP_PRODUCT({LABEL:'MP2_FORM',
                   NEW:True,
                   OP_RES:'LMP2',
                   OPERATORS:['T2^+','FOCK','T2'],
                   IDX_SV:   [1     ,2     ,3]})

# One can, for instance, define the list previouly
op_list = ['T2^+','PHI']
EXPAND_OP_PRODUCT({LABEL:['MP2_FORM'],
                   NEW:False,
                   OP_RES:'LMP2',
                   OPERATORS:op_list,
                   IDX_SV:   [1     ,2    ]})
EXPAND_OP_PRODUCT({LABEL:'MP2_FORM',
                   NEW:False,
                   OP_RES:'LMP2',
                   OPERATORS:['PHI','T2'],
                   IDX_SV:   [1     ,2    ]})
depend('H')
REPLACE({LABEL_RES:'MP2_FORM',
         LABEL_IN:'MP2_FORM',
         OP_LIST:['FOCK','H','PHI','H']})
PRINT_FORMULA({LABEL:['MP2_FORM']})

new_target('MP2_RES')
depend('MP2_FORM')
DERIVATIVE({LABEL_RES:'MP2_RES',
            LABEL_IN:'MP2_FORM',
            OP_RES:'O2',
            OP_DERIV:'T2^+'})
PRINT_FORMULA({LABEL:'MP2_RES'})

new_target('MP2_EN')
depend('MP2_FORM','MP2_RES')
FACTOR_OUT({LABEL_RES:'MP2_EN',
            LABEL_IN:'MP2_FORM',
            INTERM:'MP2_RES'})
PRINT_FORMULA({LABEL:'MP2_EN'})


new_target('MP2_LISTS')
depend('MP2_OPS')

# Just exploring some python features.
# The same effect of
#DEF_ME_LIST({LIST:'ME_LMP2',
#             OPERATOR:'LMP2',
#             IRREP:1,
#             '2MS':2,
#             AB_SYM:+1})
#
#DEF_ME_LIST({LIST:'ME_T2',
#             OPERATOR:'T2',
#             IRREP:1,
#             '2MS':2,
#             AB_SYM:+1})
#
#DEF_ME_LIST({LIST:'ME_O2',
#             OPERATOR:'O2',
#             IRREP:1,
#             '2MS':2,
#             AB_SYM:+1})
#
# can be achieved with

me_list_basis={IRREP:1,
               '2MS':0,
               AB_SYM:+1}
for op in ['LMP2','T2','O2']:          # loop over the elements of the list
    DEF_ME_LIST(dict({LIST:'ME_'+op, # merge this dictionary ...
                      OPERATOR:op},
                     **me_list_basis)) # with this, that os common for all ME lists

new_target('MP2_OPT')
depend('MP2_EN','MP2_RES','MP2_LISTS','H0')
OPTIMIZE({LABEL_OPT:'MP2_OPT',
          LABELS_IN:['MP2_RES','MP2_EN']})

new_target('DIAG')
depend('H0','MP2_OPS')
CLONE_OPERATOR({LABEL:'D2',
                TEMPLATE:'T2'})
DEF_ME_LIST(dict({LIST:'DIAG',    # merge this dictionary ...
                  OPERATOR:'D2'},
                 **me_list_basis))  # with this, that os common for all ME lists
PRECONDITIONER({LIST_PRC:'DIAG',
                LIST_INP:'H0'})


new_target('MY_TARGET',True)
# or:
# Required()
# or:
# Required(True)

depend('DIAG')
depend('MP2_OPT')

PRINT({STRING:"Now solving the equations ..."})

# Instead of a direct assignment of arguments-values, like
#
#rule('SOLVE_NLEQ',
#     {'LIST_OPT':['ME_T2'],
#      'LIST_RESID':['ME_O2'],
#      'LIST_PRC':['DIAG'],
#      'MODE':'DIA',
#      'LIST_E':['ME_LMP2'],
#      'FORM':['MP2_OPT']})
#
# one can set the dict before:
#
nleq_args={}
nleq_args[LIST_OPT] = 'ME_T2'
nleq_args[LIST_RESID] = 'ME_O2'
nleq_args[LIST_PRC] = 'DIAG'
nleq_args[MODE] = 'DIA'
nleq_args[LIST_E] = 'ME_LMP2'
nleq_args[FORM] = 'MP2_OPT'

SOLVE_NLEQ( nleq_args)

#
# Of course this is more useful in a non trivial rule,
# if you have some contiditions to one or another case.
#

PRINT_MEL({LIST:'ME_LMP2'})
PRINT({STRING:"All over ..."})

#
# This last statement is very important. It writes on disk the 
# informations about the targets, to be accessed by GeCCo
#
export_targets();
