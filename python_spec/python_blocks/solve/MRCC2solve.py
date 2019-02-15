from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


def me_list_label(root, sym, spin, ms, msc, spinadapt):
    if spinadapt :
        return "{root}G{sym:1}S{spin:0>2:}M{ms:0>2}".format(root = root,
                                                            sym=sym,
                                                            spin=spin,
                                                            ms=ms,
        )
    else:
        return {0:"{root}G{sym:1}SxxM{ms:0>2}",
                1:"{root}G{sym:1}C+1M{ms:0>2}",
                -1:"{root}G{sym:1}C-1MM{ms:0>2}",}[msc].format(root = root,
                                                               sym=sym,
                                                               ms=ms,)
def ab_sym(ms, mult):
    if (ms == 0) and ((mult-1)%4 == 0) :
        return 1
    elif (ms == 0) and ((mult+1)%4 == 0) :
        return -1
    else :
        return 0

spinadapt = keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0

if spinadapt >= 2 : # and gno >= 0
    S2_val=0
else:
    S2_val=-1  # the default
#------------------------------------------------------------------
#------------------------------------------------------------------
#SOLVE it
#------------------------------------------------------------------
#------------------------------------------------------------------

new_target('SOLVE_MRCC2')
depend('EVAL_E0')
depend('DEF_FORM_PT_LAG2')
depend('BUILD_PRECON')

debug_MEL('PRECON_LST')

#SOLVE_NLEQ uses FOPT_T2_orth for both transformations and needs X_TRM_LIST_DAG also bound to that opperator

ASSIGN_ME2OP({
        LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM'})

debug_FORM('FORM_T1_orth')

debug_FORM('FORM_T2_orth')

PRINT_MEL_INFO({LIST:'ME_T1'})
PRINT_MEL_INFO({LIST:'ME_T2g'})

SOLVE_NLEQ({
        LIST_OPT:['ME_T1','ME_T2g'],
        LIST_RESID:['ME_O1','ME_O2g'],
        LIST_PRC:['ME_PRECON1','ME_PRECON2g'],
        LIST_E:'PT_LAG_LST',
        FORM:'FOPT_PT_LAG2',
        MODE:'TRF TR0',
        FORM_SPC:['FOPT_T1_orth','FOPT_T2_orth'],
        LIST_SPC:['ME_T1_orth','ME_X_TRM','ME_X_TRM_DAG','ME_T2_orth']
        })



PUSH_RESULT({LIST:'PT_LAG_LST',COMMENT:"MRCC2", FORMAT:"SCAL F20.14"})




#--------------
new_target('SOLVE_MRCC2ref')
depend("SOLVE_MRCC2")
depend("FOPT_OMG_C0")
depend("MAKE_D0")

spinadapt=keywords.get('calculate.routes.spinadapt')
spinadapt = int(spinadapt) if spinadapt is not None else 0

ciroot=keywords.get('method.MR_P.ciroot')
ciroot = int(ciroot) if ciroot is not None else 1

maxroots = keywords.get('method.MR_P.maxroot')
maxroots = int(maxroots) if maxroots is not None else ciroot


SOLVE_map={
        LIST_OPT:'ME_C0',
        LIST_PRC:'ME_D0',
        OP_MVP:'A_C0',
        OP_SVP:'C0',
        FORM:'FOPT_OMG_C0',
        MODE:'DIA',
        N_ROOTS:maxroots,
        TARG_ROOT:ciroot
}

if (spinadapt != 0): #and refproj = 0  
    SOLVE_map[MODE]='SPP' #'DIA' will be overwritten
    SOLVE_map[LIST_SPC]='ME_C0_sp'
    SOLVE_map[FORM_SPC]='FOPT_C0_sp'

SOLVE_EVP(SOLVE_map)










#-----------------------
new_target("SOLVE_MRCC2_refopt")
depend('DEF_FORM_PT_LAG2')
depend('BUILD_PRECON')
depend("FOPT_OMG_C0")


ASSIGN_ME2OP({
        LIST:'ME_X_TRM_DAG',
        OPERATOR:'X_TRM'})

PRINT_MEL_INFO({LIST:'ME_T1'})
PRINT_MEL_INFO({LIST:'ME_T2g'})

SOLVE_NLEQ({
        LIST_OPT:['ME_T1','ME_T2g'],
        LIST_RESID:['ME_O1','ME_O2g'],
        LIST_PRC:['ME_PRECON1','ME_PRECON2g'],
        LIST_E:'PT_LAG_LST',
        FORM:'FOPT_PT_LAG2',
        MODE:'TRF TR0',
        FORM_SPC:['FOPT_T1_orth','FOPT_GAM_S','FOPT_T2_orth'],
        LIST_SPC:['ME_T1_orth','ME_X_TRM','ME_X_TRM_DAG',"ME_P_PROJ","ME_GAM_S",
                  "ME_GAM_S_ISQ",'ME_T2_orth']
        })


PUSH_RESULT({LIST:'PT_LAG_LST',COMMENT:"MRCC2", FORMAT:"SCAL F20.14"})




new_target("MAKE_GAM0_HMRCC2")
depend("FOPT_GAM0")
#depend("SOLVE_MRCC2ref")
EVALUATE({
        FORM:'FOPT_GAM0'})


new_target("MAKE_MRCC2_E")
depend('MAKE_MRCC2')
DEF_SCALAR({LABEL:"MRCC2_E"})
DEF_ME_LIST({LIST:"ME_MRCC2_E",
             OPERATOR:"MRCC2_E",
             IRREP:1,
             "2MS":0})
SCALE_COPY({LIST_RES:"ME_MRCC2_E",
            LIST_INP:"PT_LAG_LST",
            FAC:1,})



new_target("MAKE_MRCC2")
if (keywords.is_keyword_set("method.MRCC2.relaxref")):
    depend("SOLVE_MRCC2ref")
    depend("MAKE_GAM0_HMRCC2")
elif (keywords.get("calculate.solve.non_linear.optref") is not None
      and int(keywords.get("calculate.solve.non_linear.optref")) == -3 ):
    depend("SOLVE_MRCC2_refopt")
else:
    depend("SOLVE_MRCC2")
