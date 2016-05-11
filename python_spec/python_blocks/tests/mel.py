from gecco_interface import *



new_target("PrepareDiag")

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
        DIAG_TYPE:1
})

ADD_UNITY({
        LIST:'ME_DUMMY_1',
        FAC:1.0,
        INIT:True,
})


PRINT_MEL({LIST:'ME_DUMMY_1'})

new_target('TEST_SET_BLOCKS',True)
depend("PrepareDiag")

SET_BLOCKS({
        LIST:'ME_DUMMY_1',
        FAC:2.0,
        DESCR:"V,V|VV,VV",
})

PRINT_MEL({LIST:'ME_DUMMY_1'})
