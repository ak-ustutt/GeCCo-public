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

tasks=False

nc_en=4
nc_rs=2
select=True     # for nc_rs>2: select terms that contribute in SR case
#linear = True
doublet = True
cas22 = False

if doublet or cas22:
    T2_shape = 'VV,HH|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'  # skipped VVV amps here
else:
    T2_shape = 'V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'

# will later be replaced by "T2g" operators
DEF_OP_FROM_OCC({LABEL:'T2',DESCR:T2_shape})
CLONE_OPERATOR({LABEL:'L2',TEMPLATE:'T2',ADJOINT:True})

# test only
T1_shape = 'P,H|P,V|V,H'
DEF_OP_FROM_OCC({LABEL:'T1n',DESCR:T1_shape})
CLONE_OPERATOR({LABEL:'L1n',TEMPLATE:'T1n',ADJOINT:True})


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("L2(" + x + ")")

LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))

LAG_E.append(_refexp("[H,T1]"))
LAG_E.append(_refexp("[H,T2]"))
if nc_en > 1:
   LAG_E.append(_refexp("1/2*[[H,T1+T2],T1+T2]"))
if nc_en > 2:
   LAG_E.append(_refexp("1/6*[[[H,T1+T2],T1+T2],T1+T2]"))
if nc_en > 3:
   LAG_E.append(_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))

LAG_A1.append(_L1_refexp("[H,T1]"))
LAG_A1.append(_L1_refexp("[H,T2]"))
if nc_rs > 1:
    LAG_A1.append(_L1_refexp("1/2*[[H,T1+T2],T1+T2]"))    
if nc_rs > 2:
    if select:
        LAG_A1.append(_L1_refexp("1/6*[[[H,T1],T1],T1]"))
    else:
        LAG_A1.append(_L1_refexp("1/6*[[[H,T1+T2],T1+T2],T1+T2]"))
if nc_rs > 3:
    if not select:
        LAG_A1.append(_L1_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))

LAG_A2.append(_L2_refexp("[H,T1]"))
LAG_A2.append(_L2_refexp("[H,T2]"))
if nc_rs > 1:
    LAG_A2.append(_L2_refexp("1/2*[[H,T1+T2],T1+T2]"))    
if nc_rs > 2:
    if select:
        LAG_A2.append(_L2_refexp("1/6*[[[H,T1],T1],T1]"))
        LAG_A2.append(_L2_refexp("1/6*[[[H,T2],T1],T1]"))
        LAG_A2.append(_L2_refexp("1/6*[[[H,T1],T2],T1]"))
        LAG_A2.append(_L2_refexp("1/6*[[[H,T1],T1],T2]"))
    else:
        LAG_A2.append(_L2_refexp("1/6*[[[H,T1+T2],T1+T2],T1+T2]"))
if nc_rs > 3:
    if select:
        LAG_A2.append(_L2_refexp("1/24*[[[[H,T1],T1],T1],T1]"))
    else:
        LAG_A2.append(_L2_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))

LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()

REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T2','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T2','T2g']})
REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2','T2g','L2','LAM2g']})

# currently redundant:
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T1n','T2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T1n','T2g', 'L1n', 'LAM2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T1n','T2g', 'L1n', 'LAM2g']})

#REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T1t','T2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T1t','T2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T1t','T2g']})

#REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T2t','T2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T2t','T2g']})
#REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2t','T2g']})

# --- factor out densities ---
FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_E',
        LABEL_RES:'FORM_MRCC_LAG_E',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'FORM_GAM0'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

# define an operator that sums the two singles operators
if doublet or cas22:
    T1s_shape = 'P,H'
    DEF_OP_FROM_OCC({LABEL:'T1s',DESCR:T1s_shape})
else:
    CLONE_OPERATOR({LABEL:'T1s',TEMPLATE:'T1'})

CLONE_OPERATOR({LABEL:'T1s2',TEMPLATE:'T2g'}) # for squared T1s

DEF_ME_LIST({LIST:'ME_T1s',OPERATOR:'T1s',IRREP:1,'2MS':0,AB_SYM:+1})

DEF_FORMULA({LABEL:'F_T1SUM',FORMULA:'T1s=T1+T2g'})

FT1SSQ = stf.Formula("F_T1SSQ:T1s2=<T1s2'*T1*T1*T1s2'>") # it seems that "avoid" is not accepted for Formula
# as a work-around, we added this dummy term ^^^^^^^
FT1SSQ.append("<0.5*T1s2'*T1s''*T1s'''*T1s2'>", avoid=["T1s''","T1s'''"])
FT1SSQ.set_rule()
INVARIANT({LABEL_RES:'F_T1SSQ',LABEL_IN:'F_T1SSQ',OPERATORS:'T1',OP_RES:'T1s2'}) # remove the dummy term
PRINT_FORMULA({LABEL:'F_T1SSQ',MODE:'SHORT'})

EXPAND({LABEL_RES:"F_T1SSQ2",LABEL_IN:"F_T1SSQ",INTERM:"F_T1SUM"})

PRINT_FORMULA({LABEL:'F_T1SSQ2',MODE:'SHORT'})

# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM'})
# One more call for quadratic terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM'})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SSQ2'})
EXPAND({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SSQ'})


# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM'})
# One more call for quadratic terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM'})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SSQ2'})
EXPAND({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SSQ'})


# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM'})
# One more call for quadratic terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM'})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SSQ2'})
EXPAND({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SSQ'})

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

# Replace in singles residual as much as possible by singles part of doubles residual
# in cases where T1 and T2 are treated on the same footing, this replaces everything
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_Amp1',LABEL_RES:'FORM_MRCC_LAG_Amp1',INTERM:'FORM_MRCC_LAG_Amp2'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})


K4E = True
if (K4E):
    # generate the intermediate to K4E -> INTkx
    DEF_OP_FROM_OCC({LABEL:'INTkx',JOIN:2,DESCR:',VV;PP,|,V;PP,H|,;PP,HH'})
    
    DEF_ME_LIST({LIST:'ME_INTkx',OPERATOR:'INTkx',IRREP:1,'2MS':0,AB_SYM:+1})

    # define a formal Hamiltonian that only contains the blocks of interest
    DEF_OP_FROM_OCC({LABEL:'Hpppp',DESCR:'PP;PP'}) 

    # replace relevant blocks of Hamiltonian with Hpppp
    REPLACE({LABEL_RES:'F_preINTkx',LABEL_IN:'FORM_MRCC_LAG_Amp2',OP_LIST:['H','Hpppp']})
    # remove all terms that do not include Hpppp
    INVARIANT({LABEL_RES:'F_preINTkx',LABEL_IN:'F_preINTkx',OPERATORS:'H',OP_RES:'O2g'})
    # now, take derivative with respct to Hpppp
    DERIVATIVE({LABEL_RES:'F_INTkx',LABEL_IN:'F_preINTkx',OP_RES:'INTkx',OP_DERIV:'Hpppp'})
    REORDER_FORMULA({LABEL_IN:'F_INTkx',LABEL_RES:'F_INTkx'})

    PRINT_FORMULA({LABEL:'F_INTkx',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_Amp1',
        LABEL_RES:'FORM_MRCC_LAG_Amp1',
        INTERM:'F_INTkx'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_Amp2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        INTERM:'F_INTkx'})
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})

### Add here more terms to be factored out ...
if (True):
    # try this for 3externals
    # has the shape of O2g
    CLONE_OPERATOR({LABEL:'INT3ext',TEMPLATE:'O2g'})
    DEF_ME_LIST({LIST:'ME_INT3ext',OPERATOR:'INT3ext',IRREP:1,'2MS':0,AB_SYM:+1})

    # define a formal Hamiltonian that only contains the blocks of interest
    DEF_OP_FROM_OCC({LABEL:'Hppph',DESCR:'PP;PH|PP;PV'})

    # replace relevant blocks of Hamiltonian with Hppph
    REPLACE({LABEL_RES:'F_preINT3ext',LABEL_IN:'FORM_MRCC_LAG_Amp2',OP_LIST:['H','Hppph']})
    # remove all terms that do not include Hppph
    INVARIANT({LABEL_RES:'F_preINT3ext',LABEL_IN:'F_preINT3ext',OPERATORS:'H',OP_RES:'O2g'})
    # now, take derivative with respct to Hppph
    DERIVATIVE({LABEL_RES:'F_INT3ext',LABEL_IN:'F_preINT3ext',OP_RES:'INT3ext',OP_DERIV:'Hppph'})
    REORDER_FORMULA({LABEL_IN:'F_INT3ext',LABEL_RES:'F_INT3ext'})

    PRINT_FORMULA({LABEL:'F_INT3ext',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_Amp1',
        LABEL_RES:'FORM_MRCC_LAG_Amp1',
        INTERM:'F_INT3ext'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_Amp2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        INTERM:'F_INT3ext'})

    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})
    

OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        LABELS_IN:['F_T1SUM','F_INTkx','F_INT3ext','FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp1',MODE:'SHORT'}) # only dummy
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_Amp2',MODE:'SHORT'})

if not tasks:
    filename = 'icmrcc_mrccsd_'+str(nc_en)+str(nc_rs)
else:
    filename = 'auicmrcc_mrccsd_'+str(nc_en)+str(nc_rs)    
if nc_rs>2 and select:
    filename = filename + 's'
if doublet:
    filename = filename + '_doublet'
elif cas22:
    filename = filename + '_cas22'
filename2 = filename + '.formulae'
filename = filename + '.itfaa'


TRANSLATE_ITF({
        LABEL:'FOPT_MRCC_LAG',
        OUTPUT:filename,
        TITLE:filename2,
        MULTI:True,
        PROCESS:True,
        KEXT:True,
        TASKS:tasks,
        INIT_RES:False,
        ITIN:True,
        RENAME:['MRCC_LAG','ECC','T1','T1','T2g','T2','O1','R1','O2g','R2','GAM0','Ym<RANK>'],
        CODE:['<Sum_T1>','T1s','<Update_INTkx>','INTkx','<Residual>','INT3ext','MRCC_LAG','O1','O2g']})

#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
