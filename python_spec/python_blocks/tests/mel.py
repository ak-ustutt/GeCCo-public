from gecco_interface import *




new_target('TEST_MODIFY_BLOCKS',True)

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



ADD_UNITY({
        LIST:'ME_DUMMY_1',
        FAC:1.0,
        INIT:True,
})




PRINT_MEL({LIST:'ME_DUMMY_1'})


PRINT({STRING:"Shifting by 0.5"})
MODIFY_BLOCKS({
        LIST:'ME_DUMMY_1',
        FAC:0.5,
        DESCR:"P,P|VV,VV|HH,HH",
        MODE:"SHIFT",
})




PRINT_MEL({LIST:'ME_DUMMY_1'})

PRINT({STRING:"Setting to 2.5"})

MODIFY_BLOCKS({
        LIST:'ME_DUMMY_1',
        FAC:2.5,
        DESCR:"PP,PP|VV,VV|H,H",
        MODE:"SET",
})
PRINT_MEL({LIST:'ME_DUMMY_1'})
PRINT({STRING:"Scaling by 4"})

MODIFY_BLOCKS({
        LIST:'ME_DUMMY_1',
        FAC:4.0,
        DESCR:"P,P|V,V|VV,VV",
        MODE:"SCALE",
})

PRINT_MEL({LIST:'ME_DUMMY_1'})
