from python_interface.gecco_modules.NoticeUtil import *
from python_interface.gecco_interface import *
import python_interface.gecco_modules.string_to_form as stf

#[(A_c',c   A_c',t1  A_c',t2g  )     (1  0   0  ) ](r_c)    0
#[(A_l1,c   A_l1,t1  A_l1,t2g  ) - w (0  S11 S12) ](r_t1)  =0
#[(A_l2g,c  A_l2g,t1 A_l2g,t2g )     (0  S21 S22) ](r_t2g)  0



lagrangian = keywords.get('method.MRCC2.lagrangian')
lag_type = int(lagrangian) if lagrangian is not None else 4 
print("lagrangian:",lag_type,type(lag_type))





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

debug_FORM("FORM_AR1_RSPNS_q")

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


debug_FORM("FORM_AR2g_RSPNS_q")


new_target("DEF_FORM_AR_RSPNS_mu")
depend("DEF_FORM_PT_LAG2")
depend("MAKE_MRCC2_E")
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
            OPERATORS:['AR_rspns_mu',"MRCC2_E",'R_mu','AR_rspns_mu'],
            IDX_SV:[1,2,3,1],
            FAC:-1.0})

debug_FORM("FORM_AR_RSPNS_mu")


new_target("DEF_FORMS_METRIC")
depend('DEF_RESPONSE_OPs')








SR1=stf.GenForm("FORM_SR1t","Scal_S1")
if lag_type <=0:
    raise Exception("unknown lagrangian")

if lag_type >=1:
    SR1+="<C0^+*(LAM1)*(R1_q+R2g_q+1/2[R1_q+R2g_q,T1+T2g])*C0>"
if lag_type >=2:
    SR1+="1/6<C0^+*(LAM1)*([[R1_q+R2g_q,T1+T2g],T1]+[[R1_q+R2g_q,T1],T2g]+[[R1_q,T2g],T2g])*C0>"
if lag_type >=3:
    SR1+="1/24<C0^+*(LAM1)*( [[[R1_q+R2g_q,T1+T2g],T1],T1]+[[[R1_q+R2g_q,T1],T2g],T1]+[[[R1_q+R2g_q,T1],T1],T2g]+"\
        "[[[R1_q,T2g],T2g],T1]+[[[R1_q,T2g],T1],T2g]+[[[R1_q,T1],T2g],T2g])*C0>"
if lag_type >=4:
    SR1+="1/120<C0^+*(LAM1)*( [[[[R1_q+R2g_q,T1+T2g],T1],T1],T1]+[[[[R1_q+R2g_q,T1],T2g],T1],T1]+[[[[R1_q+R2g_q,T1],T1],T2g],T1]+"\
        "[[[[R1_q+R2g_q,T1],T1],T1],T2g]+"\
        "[[[[R1_q,T2g],T2g],T1],T1]+[[[[R1_q,T2g],T1],T2g],T1]+[[[[R1_q,T2g],T1],T1],T2g]+"\
        "[[[[R1_q,T1],T2g],T2g],T1]+[[[[R1_q,T1],T2g],T1],T2g]+[[[[R1_q,T1],T1],T2g],T2g]"\
        ")*C0>"
if lag_type >=5:
    raise Exception("unknown lagrangian")


print "\n".join([str(br) for br in SR1.show()])

SR1.set_rule()


SELECT_SPECIAL({LABEL_RES:'FORM_DENS1_q',
                LABEL_IN:'FORM_SR1t',
                TYPE:'nonzero',         
                MODE:'sum'})





DERIVATIVE({LABEL_RES:'FORM_S1',
            LABEL_IN:'FORM_DENS1_q',
            OP_RES:'SR1_rspns_q',
            OP_DERIV:['LAM1']})

debug_FORM("FORM_S1")




if lag_type <=0:
    raise Exception("unknown lagrangian")

SR2=stf.GenForm("FORM_SR2t","Scal_S2")
if lag_type >=1:
    SR2+="<C0^+*(LAM2g)*(R1_q+R2g_q+1/2[R1_q+R2g_q,T1]+1/2[R1_q,T2g])*C0>"
if lag_type >=2:
    SR2+="1/6<C0^+*(LAM2g)*([[R1_q+R2g_q,T1],T1]+[[R1_q,T2g],T1]+[[R1_q,T1],T2g])*C0>"
if lag_type >=3:
    SR2+="1/24<C0^+*(LAM2g)*([[[R1_q+R2g_q,T1],T1],T1]"\
        "+[[[R1_q,T2g],T1],T1]+[[[R1_q,T1],T2g],T1]+[[[R1_q,T1],T1],T2g])*C0>"
if lag_type >=4:
    SR2+="1/120<C0^+*(LAM2g)*([[[[R1_q+R2g_q,T1],T1],T1],T1]"\
        "[[[[R1_q,T1+T2g],T1],T1],T1]+[[[[R1_q,T1],T2g],T1],T1]+[[[[R1_q,T1],T1],T2g],T1]+[[[[R1_q,T1],T1],T1],T2g])*C0>"
if lag_type >=5:
    raise Exception("unknown lagrangian")

print "\n".join([str(br) for br in SR2.show()])
    
SR2.set_rule()

SELECT_SPECIAL({LABEL_RES:'FORM_DENS2_q',
                LABEL_IN:'FORM_SR2t',
                TYPE:'nonzero',
                MODE:'sum'})

DERIVATIVE({LABEL_RES:'FORM_S2g',
            LABEL_IN:'FORM_DENS2_q',
            OP_RES:'SR2g_rspns_q',
            OP_DERIV:['LAM2g']})




debug_FORM("FORM_S2g")







EXPAND_OP_PRODUCT({OPERATORS:['SR_rspns_mu','R_mu','SR_rspns_mu'],
                   IDX_SV:[1,2,1],
                   LABEL:'FORM_SR_RSPNS_mu',
                   NEW:True, 
                   OP_RES:'SR_rspns_mu'})



debug_FORM("FORM_SR_rspns_mu")

#new_target("DEF_FORM_RSPN_MRCC2")
