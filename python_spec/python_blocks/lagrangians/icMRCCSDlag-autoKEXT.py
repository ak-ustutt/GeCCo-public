"""An implementation of the icMRCCSD Lagrangian, done term by term and separated T1 and T2.

 **** THIS IS A TEST VERSION FOR INITIAL TESTS OF THE ITF TRANSLATOR ****
 **** DO NOT MERGE INTO THE MAIN BRANCH !!                           ****

History:

Yuri august 2017: Creation based on MRCC2lag.py. Implementation up to maxcom=2

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import ref_relaxation

i_am="icMRCCSDlag.py"

new_target('DEF_FORM_MRCC_LAG')
heading('Defining the icMRCC Lagrangian')

depend('DEF_T')
#depend('DEF_T2g')
#depend('DEF_T1')

depend('DEF_LAM')
#depend('DEF_LAM2g')
#depend('DEF_LAM1')

depend('DEF_O')
#depend('DEF_O2g')
#depend('DEF_O1')


depend('GAM0_CALC')

DEF_SCALAR({
        LABEL:'MRCC_LAG'})

DEF_ME_LIST({
        LIST:'MRCC_LAG_LST',
        OPERATOR:'MRCC_LAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

DEF_SCALAR({
        LABEL:'MRCC_LAG_A1'})

DEF_SCALAR({
        LABEL:'MRCC_LAG_A2'})


#T2_shape = 'V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'
#T2_shape = 'PP,HH|VV,HH|PV,VV|PV,HH|PP,VV|PP,HV|VV,VH|PV,HV'
T2_shape = 'PP,HH|VV,HH|PV,VV|PV,HH|PP,VV|PP,HV|PV,HV|VV,VH'
#T2_shape = 'PV,HV'
DEF_OP_FROM_OCC({LABEL:'T2',DESCR:T2_shape})
CLONE_OPERATOR({LABEL:'L2',TEMPLATE:'T2',ADJOINT:True})

T1_shape = 'P,H|P,V|V,H'
#T1_shape = 'P,H'
DEF_OP_FROM_OCC({LABEL:'T1n',DESCR:T1_shape})
CLONE_OPERATOR({LABEL:'L1n',TEMPLATE:'T1n',ADJOINT:True})

T1t_shape = 'P,H|P,V|V,H'
#T1t_shape = 'P,H'
DEF_OP_FROM_OCC({LABEL:'T1t',DESCR:T1_shape})
CLONE_OPERATOR({LABEL:'L1t',TEMPLATE:'T1t',ADJOINT:True})

T2t_shape = 'PV,HH'
#T1t_shape = 'P,H'
DEF_OP_FROM_OCC({LABEL:'T2t',DESCR:T1_shape})
CLONE_OPERATOR({LABEL:'L2t',TEMPLATE:'T2t',ADJOINT:True})


doublet = True
if doublet:
    # Every term in the Lagrangian is enclosed by <C0^+ and C0>
    def _refexp(x):
        return "<C0^+*(" + x + ")*C0>"

    # The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
    def _L1_refexp(x):
        return _refexp("L1n(" + x + ")")

    def _L2_refexp(x):
        return _refexp("L2(" + x + ")")

    LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H-H"))

    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L1_refexp("H"))
    LAG_A2.append(_L2_refexp("H"))

    LAG_E.append(_refexp("[H,T1n]"))
    LAG_E.append(_refexp("[H,T2]"))
    LAG_E.append(_refexp("1/2*[[H,T1n+T2],T1n+T2]"))
#    LAG_E.append(_refexp("1/2*[[H,T1n],T1n]"))

#    LAG_A1.append(_L1_refexp("[H,T1n]"))
#    LAG_A1.append(_L1_refexp("[H,T2]"))

    LAG_A2.append(_L1_refexp("[H,T1n]"))
    LAG_A2.append(_L1_refexp("[H,T2]"))
#    LAG_A2.append(_L1_refexp("1/2*[[H,T2t],T2t]"))    
    LAG_A2.append(_L2_refexp("[H,T1n]"))
    LAG_A2.append(_L2_refexp("[H,T2]"))
#    LAG_A2.append(_L2_refexp("1/2*[[H,T2t],T2t]"))    

    LAG_E.set_rule()
    LAG_A1.set_rule()
    LAG_A2.set_rule()

else:

    # Every term in the Lagrangian is enclosed by <C0^+ and C0>
    def _refexp(x):
        return "<C0^+*(" + x + ")*C0>"

    # The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
    def _L1_refexp(x):
        #return _refexp("LAM1(" + x + ")")
        return _refexp("L1n(" + x + ")")

    def _L2_refexp(x):
        #return _refexp("LAM2g(" + x + ")")
        return _refexp("L2(" + x + ")")

    def _L3_refexp(x):
        return _refexp("-LAM2g(" + x + ")")

    def _L4_refexp(x):
        return _refexp("-LAM1(" + x + ")")

    LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H-H"))
    #LAG_A1.append(_L4_refexp("H"))

    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))
    #LAG_A2.append(_L3_refexp("H"))

    #LAG_E.append(_refexp("[H,T1n]"))
    #LAG_E.append(_refexp("[H,T2]"))

    #LAG_E.append(_refexp("0.5*[[H,T1t],T1t]"))
    #LAG_E.append(_refexp("0.5*[[H,T1n],T2]"))
    #LAG_E.append(_refexp("0.5*[[H,T2],T1n]"))
    #LAG_E.append(_refexp("0.5*[[H,T2],T2]"))

    #LAG_A1.append(_L1_refexp("[H,T1n]"))
    #LAG_A1.append(_L1_refexp("[H,T2]"))
    #LAG_A1.append(_L1_refexp("(1/2)*[[H,T1],T2g]"))

    #LAG_A2.append(_L1_refexp("[H,T1n]"))
    #LAG_A2.append(_L1_refexp("[H,T2]"))
    #LAG_A2.append(_L2_refexp("[H,T1n]"))
    #LAG_A2.append(_L2_refexp("[H,T2]"))
    #LAG_A2.append(_L2_refexp("(1/2)*[[H,T2g],T2g]"))



    LAG_E.set_rule()
    LAG_A1.set_rule()
    LAG_A2.set_rule()



REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T2','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T2','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2','T2g','L2','LAM2g']})

REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T1n','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T1n','T2g', 'L1n', 'LAM2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T1n','T2g', 'L1n', 'LAM2g']})

REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T1t','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T1t','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T1t','T2g']})

REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T2t','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T2t','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2t','T2g']})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_E',
        LABEL_RES:'FORM_MRCC_LAG_E',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'FORM_GAM0'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

#Make the Derivative with respect to LAM
# only dummy
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_Amp1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_LAG_Amp1',LABEL_RES:'FORM_MRCC_LAG_Amp1'})

DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_LAG_Amp2',LABEL_RES:'FORM_MRCC_LAG_Amp2'})


K4E = True
if (K4E):
    # generate the intermediate to K4E -> INTkx
    # has the shape of O2g
    #CLONE_OPERATOR({LABEL:'INTpppp',TEMPLATE:'O2g'})
    DEF_OP_FROM_OCC({LABEL:'INTkx',JOIN:2,DESCR:',VV;PP,|,V;PP,H|,;PP,HH'})
    
    DEF_ME_LIST({LIST:'ME_INTkx',OPERATOR:'INTkx',IRREP:1,'2MS':0,AB_SYM:+1})

    # define a formal Hamiltonian that only contains the blocks of interest
    DEF_OP_FROM_OCC({LABEL:'Hpppp',DESCR:'PP;PP'}) 

    # replace relevant blocks of Hamiltonian with Hpppp
    REPLACE({LABEL_RES:'F_preINTpppp',LABEL_IN:'FORM_MRCC_LAG_Amp2',OP_LIST:['H','Hpppp']})
    # remove all terms that do not include Hpppp
    INVARIANT({LABEL_RES:'F_preINTpppp',LABEL_IN:'F_preINTpppp',OPERATORS:'H',OP_RES:'O2g'})
    # now, take derivative with respct to Hpppp
    DERIVATIVE({LABEL_RES:'F_INTpppp',LABEL_IN:'F_preINTpppp',OP_RES:'INTkx',OP_DERIV:'Hpppp'})
    REORDER_FORMULA({LABEL_IN:'F_INTpppp',LABEL_RES:'F_INTpppp'})

    PRINT_FORMULA({LABEL:'F_INTpppp',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_Amp2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        INTERM:'F_INTpppp'})
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})

    

OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        #LABELS_IN:['CCD_FORM_PP','FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})
        LABELS_IN:['F_INTpppp','FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp1',MODE:'SHORT'}) # only dummy
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})


filename = 'icmrcc_mrccsd_11'
if doublet:
    filename = filename + '_doublet'
filename = filename + '.itfaa'


TRANSLATE_ITF({
        LABEL:'FOPT_MRCC_LAG',
        OUTPUT:filename,
        TITLE:'icmrcc_mrccsd_11_doublet.formulae',
        MULTI:True,
        PROCESS:True,
        KEXT:True,
        TASKS:False,
        INIT_RES:False,
        ITIN:True,
        RENAME:['MRCC_LAG','ECC','T1','T','T2g','T','O1','R','O2g','R','GAM0','Ym<RANK>'],
        CODE:['<Update_INTkx>','INTkx','<Residual>','MRCC_LAG','O1','O2g']})

#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
