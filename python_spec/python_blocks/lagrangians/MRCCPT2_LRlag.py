from gecco_modules.NoticeUtil import *
from gecco_interface import *
import gecco_modules.string_to_form as stf

#[(A_c',c   A_c',t1  A_c',t2g  )     (1  0   0  ) ](r_c)    0
#[(A_l1,c   A_l1,t1  A_l1,t2g  ) - w (0  S11 S12) ](r_t1)  =0
#[(A_l2g,c  A_l2g,t1 A_l2g,t2g )     (0  S21 S22) ](r_t2g)  0







new_target("DEF_FORM_AR_RSPNS_q")
depend("DEF_FORM_PT_LAG2")
depend('DEF_RESPONSE_OPs')

DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A1_RAW',
        LABEL_RES:'FORM_AR1_RSPNS_q_INT',
        OP_RES:'O1',
        OP_DERIV:'LAM1'
        })

DERIVATIVE({
        LABEL_IN:'FORM_AR1_RSPNS_q_INT',
        LABEL_RES:'FORM_AR1_RSPNS_q',
        OP_RES:'AR1_rspns_q',
        OP_DERIV:['T1','T2g','C0'],
        OP_MULT:['R1_q','R2g_q','R_mu']})


SUM_TERMS({
        LABEL_IN:"FORM_AR1_RSPNS_q",
        LABEL_RES:"FORM_AR1_RSPNS_q",
       })

debug_FORM("FORM_AR1_RSPNS_q", only_this=True)

DERIVATIVE({
        LABEL_IN:'FORM_PT_LAG_A2_RAW',
        LABEL_RES:'FORM_AR2g_RSPNS_q_INT',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'
        })


DERIVATIVE({
        LABEL_IN:'FORM_AR2g_RSPNS_q_INT',
        LABEL_RES:'FORM_AR2g_RSPNS_q',
        OP_RES:'AR2g_rspns_q',
        OP_DERIV:['T1','T2g','C0'],
        OP_MULT:['R1_q','R2g_q','R_mu']})

SUM_TERMS({
        LABEL_IN:"FORM_AR2g_RSPNS_q",
        LABEL_RES:"FORM_AR2g_RSPNS_q",
       })


debug_FORM("FORM_AR2g_RSPNS_q", only_this=True)


new_target("DEF_FORM_AR_RSPNS_mu")
depend("DEF_FORM_PT_LAG2")
depend("EVAL_E0")
depend('DEF_RESPONSE_OPs')

DERIVATIVE({LABEL_RES:'FORM_AR_RSPNS_mu_INT',
            LABEL_IN:'FORM_PT_LAG_E',
            OP_RES:'AR_rspns_mu',
            OP_DERIV:['C0^+']})


DERIVATIVE({LABEL_RES:'FORM_AR_RSPNS_mu',
            LABEL_IN:'FORM_AR_RSPNS_mu_INT',
            OP_RES:'AR_rspns_mu',
            OP_DERIV:['T1','T2g','C0'],
            OP_MULT:['R1_q','R2g_q','R_mu']})

EXPAND_OP_PRODUCT({
            LABEL:'FORM_AR_RSPNS_mu',
            NEW:False,
            OP_RES:'AR_rspns_mu',
            OPERATORS:['AR_rspns_mu','E0','R_mu','AR_rspns_mu'],
            IDX_SV:[1,2,3,1],
            FAC:-1.0})

debug_FORM("FORM_AR_RSPNS_mu", only_this=True)


new_target("DEF_FORMS_METRIC")
depend('DEF_RESPONSE_OPs')


stf.Formula("FORM_SR1:Scal_S1="\
            "<C0^+*(LAM1)*("\
            "R1_q+R2g_q+"\
            "[R1_q+R2g_q,T1+T2g]+"\
            "1/2*[[R1_q+R2g_q,T1+T2g],T1]+1/2*[[R1_q+R2g_q,T1],T2g]+"\
            "1/2*[[R1_q,T2g],T2g]+"\
            "1/6*[[[R1_q+R2g_q,T1+T2g],T1],T1]+1/6*[[[R1_q+R2g_q,T1],T2g],T1]+1/6*[[[R1_q+R2g_q,T1],T1],T2g]+"\
            "1/6*[[[R1_q,T2g],T2g],T1]+1/6*[[[R1_q,T2g],T1],T2g]+1/6*[[[R1_q,T1],T2g],T2g]+"\
            "1/24*[[[[R1_q+R2g_q,T1+T2g],T1],T1],T1]+1/24*[[[[R1_q+R2g_q,T1],T2g],T1],T1]+1/24*[[[[R1_q+R2g_q,T1],T1],T2g],T1]+"\
            "1/24*[[[[R1_q+R2g_q,T1],T1],T1],T2g]+"
            "1/24*[[[[R1_q,T2g],T2g],T1],T1]+1/24*[[[[R1_q,T2g],T1],T2g],T1]+1/24*[[[[R1_q,T2g],T1],T1],T2g]+"\
            "1/24*[[[[R1_q,T1],T2g],T2g],T1]+1/24*[[[[R1_q,T1],T2g],T1],T2g]+1/24*[[[[R1_q,T1],T1],T2g],T2g]"\
            ")*C0>").set_rule()

#stf.Formula("FORM_SR1:Scal_S1="\
#            "<C0^+*(LAM1)*("\
#            "R1_q+"\
#            "[R1_q,T1+T2g]+"\
#            "1/2*[[R1_q,T1+T2g],T1]+1/2*[[R1_q,T1],T2g]+"\
#            "1/2*[[R1_q,T2g],T2g]+"\
#            "1/6*[[[R1_q,T1+T2g],T1],T1]+1/6*[[[R1_q,T1],T2g],T1]+1/6*[[[R1_q,T1],T1],T2g]+"\
#            "1/6*[[[R1_q,T2g],T2g],T1]+1/6*[[[R1_q,T2g],T1],T2g]+1/6*[[[R1_q,T1],T2g],T2g]+"\
#            "1/24*[[[[R1_q,T1+T2g],T1],T1],T1]+1/24*[[[[R1_q,T1],T2g],T1],T1]+1/24*[[[[R1_q,T1],T1],T2g],T1]+"\
#            "1/24*[[[[R1_q,T1],T1],T1],T2g]+"
#            "1/24*[[[[R1_q,T2g],T2g],T1],T1]+1/24*[[[[R1_q,T2g],T1],T2g],T1]+1/24*[[[[R1_q,T2g],T1],T1],T2g]+"\
#            "1/24*[[[[R1_q,T1],T2g],T2g],T1]+1/24*[[[[R1_q,T1],T2g],T1],T2g]+1/24*[[[[R1_q,T1],T1],T2g],T2g]"\
#            ")*C0>").set_rule()



SELECT_SPECIAL({LABEL_RES:'FORM_DENS1_q',
                LABEL_IN:'FORM_SR1',
                TYPE:'nonzero',         
                MODE:'sum'})

#EXPAND_OP_PRODUCT({OPERATORS:['LAM1','R1_q'],
#                   IDX_SV:[1,2,1],
#                   LABEL:'FORM_DENS1_q',
#                   NEW:True, 
#                   OP_RES:'Scal_S1'})





DERIVATIVE({LABEL_RES:'FORM_S1',
            LABEL_IN:'FORM_DENS1_q',
            OP_RES:'SR1_rspns_q',
            OP_DERIV:['LAM1']})

debug_FORM("FORM_S1")







stf.Formula("FORM_SR2:Scal_S2="\
            "<C0^+*(LAM2g)*("\
            "R1_q+R2g_q+"\
            "[R1_q+R2g_q,T1]+"\
            "[R1_q,T2g]+"\
            "1/2*[[R1_q+R2g_q,T1],T1]"\
            "1/2*[[R1_q,T2g],T1]+1/2*[[R1_q,T1],T2g]+"\
            "1/6*[[[R1_q+R2g_q,T1],T1],T1]"\
            "1/6*[[[R1_q,T2g],T1],T1]+1/6*[[[R1_q,T1],T2g],T1]+1/6*[[[R1_q,T1],T1],T2g]+"\
            "1/24*[[[[R1_q+R2g_q,T1],T1],T1],T1]+"\
            "1/24*[[[[R1_q,T1+T2g],T1],T1],T1]+1/24*[[[[R1_q,T1],T2g],T1],T1]+1/24*[[[[R1_q,T1],T1],T2g],T1]"\
            ")*C0>").set_rule()



SELECT_SPECIAL({LABEL_RES:'FORM_DENS2_q',
                LABEL_IN:'FORM_SR2',
                TYPE:'nonzero',
                MODE:'sum'})

#EXPAND_OP_PRODUCT({OPERATORS:['LAM2g','R2g_q'],
#                   IDX_SV:[1,2,1],
#                   LABEL:'FORM_DENS2_q',
#                   NEW:True, 
#                   OP_RES:'Scal_S2'})


DERIVATIVE({LABEL_RES:'FORM_S2g',
            LABEL_IN:'FORM_DENS2_q',
            OP_RES:'SR2g_rspns_q',
            OP_DERIV:['LAM2g']})




debug_FORM("FORM_S2g",only_this=True)







EXPAND_OP_PRODUCT({OPERATORS:['SR_rspns_mu','R_mu','SR_rspns_mu'],
                   IDX_SV:[1,2,1],
                   LABEL:'FORM_SR_RSPNS_mu',
                   NEW:True, 
                   OP_RES:'SR_rspns_mu'})



debug_FORM("FORM_SR_rspns_mu")

#new_target("DEF_FORM_RSPN_MRCC2")
