from python_interface.gecco_interface import *

new_target('TEST_ADD_UNITY',True)



DEF_OP_FROM_OCC({
        LABEL:"DUMMY_1",
        DESCR:'P,P|PP,PP|V,V|VV,VV|H,H|HH,HH'
})

SET_HERMITIAN({
        LABEL:"DUMMY_1",
        CA_SYMMETRY:+1})

DEF_ME_LIST({
        LIST:'ME_DUMMY_1',
        OPERATOR:'DUMMY_1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1,
        DIAG_TYPE:1,
        MAX_REC:3,
        MIN_REC:1,
        REC:2
})

DEF_ME_LIST({
        LIST:'ME_DUMMY_2',
        OPERATOR:'DUMMY_1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1,
        DIAG_TYPE:1,
        MAX_REC:3,
        MIN_REC:1,
        REC:2
})


ADD_UNITY({
        LIST:'ME_DUMMY_1',
        FAC:0.5,
        INIT:True,
        MS_SYM_SIGN:1})

PRINT({STRING:"Mode square"})


SCALE_COPY({LIST_RES:'ME_DUMMY_2',
           LIST_INP:'ME_DUMMY_1',
           FAC:3,
           MODE:'square'
})


PRINT_MEL({LIST:'ME_DUMMY_2'})

PRINT({STRING:"Mode prc-thresh"})

SCALE_COPY({LIST_RES:'ME_DUMMY_2',
           LIST_INP:'ME_DUMMY_1',
           FAC:0.8,
           MODE:'prc-thresh'
})


PRINT_MEL({LIST:'ME_DUMMY_2'})

PRINT({STRING:"Mode scale"})

SCALE_COPY({LIST_RES:'ME_DUMMY_2',
           LIST_INP:'ME_DUMMY_1',
           FAC:2.0,
           MODE:'scale'
})


PRINT_MEL({LIST:'ME_DUMMY_2'})

PRINT({STRING:"Mode precond"})

# PReparing "preconditioner"
SCALE_COPY({LIST_RES:'ME_DUMMY_2',
           LIST_INP:'ME_DUMMY_2',
           FAC:0.0,
           MODE:'scale'
})

SCALE_COPY({LIST_RES:'ME_DUMMY_2',
           LIST_INP:'ME_DUMMY_2',
           FAC:0.5,
           MODE:'prc-thresh'
})

SCALE_COPY({LIST_RES:'ME_DUMMY_1',
           LIST_INP:'ME_DUMMY_2',
           FAC:1.0,
           MODE:'precond'
})


PRINT_MEL({LIST:'ME_DUMMY_1'})


#
#DEF_OP_FROM_OCC({
#        LABEL:"DUMMY_2",
#        DESCR:'P,H|H,V|H,P'
#})


#DEF_ME_LIST({
#        LIST:'ME_DUMMY_2',
#        OPERATOR:'DUMMY_2',
#        IRREP:1,
#        '2MS':0,
#        AB_SYM:+1,
#        DIAG_TYPE:1,
# #       MAX_REC:3,
#        MIN_REC:1,
#        REC:2
#})
#

#ADD_UNITY({
# #       LIST:'ME_DUMMY_2',
#        FAC:1.0,
#        INIT:True,
#        MS_SYM_SIGN:-1
#})
#
#PRINT_MEL({LIST:'ME_DUMMY_2'})
