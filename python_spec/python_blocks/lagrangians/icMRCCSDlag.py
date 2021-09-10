"""An implementation of the icMRCCSD Lagrangian, done term by term and separated T1 and T2.

 **** THIS IS A TEST VERSION FOR INITIAL TESTS OF THE ITF TRANSLATOR ****
 **** DO NOT MERGE INTO THE MAIN BRANCH !!                           ****

History:

Yuri august 2017: Creation based on MRCC2lag.py. Implementation up to maxcom=2

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import ref_relaxation


#===================================================================================#
# quick helper routine
#===================================================================================#
def create_plist(n1,sym1,n2,sym2):
    plist = []
    # create an entry with all set to sym2
    entry0 = []
    for ii in range(n1+n2):
        entry0.append(sym2)

    #print('e0:', entry0)
    make_plist_rec(entry0,sym1,0,0,n1,plist)
    
    return plist
    
def make_plist_rec(curentry,sym,curlevel,iimin,maxlevel,plist):
    
    #print('>> entered level ',curlevel,' of ',maxlevel)
    
    newentry = list(curentry)
    
    if (iimin>len(newentry)):
        return

    if (curlevel==maxlevel):
        #print ('===========> writing ',newentry)
        plist.append(newentry)
        return
    
    for ii in range(iimin,len(newentry)):
        
        
        if newentry[ii]!=sym:
            vnewentry = list(newentry)
            vnewentry[ii] = sym

            #print('level: ',curlevel,' ii = ',ii,vnewentry)
            
            make_plist_rec(vnewentry,sym,curlevel+1,ii,maxlevel,plist)
        #else:
        #    print ('level: ',curlevel,'ii = ',ii)

#===================================================================================#


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
nc_rs=4
select=True     # for nc_rs>2: select terms that contribute in SR case
#select = False
#linear = True
doublet = False
cas22 = True

remove_gamma0 = True # remove the scalar part of GAM0 (which is just 1.0)

# make sure that GAM0 is recognized as a Hermitian operator (for transpose):
SET_HERMITIAN({LABEL:'GAM0',CA_SYMMETRY:+1})

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
        for nsingles in range(5):
            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
            for entryT in listT:
                print "Generating: "+_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]")
                LAG_A1.append(_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))

#        LAG_A1.append(_L1_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 4:
    if not select:
        for nsingles in range(6):
            listT = create_plist(nsingles,'T1',5-nsingles,'T2')
            for entryT in listT:
                print "Generating: "+_L1_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]")
                LAG_A1.append(_L1_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]"))
#        LAG_A1.append(_L1_refexp("1/120*[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 5:
    if not select:
        LAG_A1.append(_L1_refexp("1/720*[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
# I think that singles can accomodate at most 6-fold
#if nc_rs > 6:
#    if not select:
#        LAG_A1.append(_L1_refexp("1/5040*[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
#if nc_rs > 7:
#    if not select:
#        LAG_A1.append(_L1_refexp("1/40320*[[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))

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
        for nsingles in range(5):
            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
            for entryT in listT:
                print "Generating: "+_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]")
                LAG_A2.append(_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))

        #LAG_A2.append(_L2_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 4:
    if not select:
        for nsingles in range(6):
            listT = create_plist(nsingles,'T1',5-nsingles,'T2')
            for entryT in listT:
                print "Generating: "+_L2_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]")
                LAG_A2.append(_L2_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]"))

       # LAG_A2.append(_L2_refexp("1/120*[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 5:
    if not select:
        LAG_A2.append(_L2_refexp("1/720*[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 6:
    if not select:
        LAG_A2.append(_L2_refexp("1/5040*[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
if nc_rs > 7:
    if not select:
        LAG_A2.append(_L2_refexp("1/40320*[[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))

PRINT({STRING:"Now expanding energy"})
LAG_E.set_rule()
PRINT({STRING:"Now expanding singles projection"})
LAG_A1.set_rule()
PRINT({STRING:"Now expanding doubles projection"})
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

# currently, this messes up the FACTOR_OUT routines
# remove the scalar GAM0 part only at the end 
# remove fully contracted part of C0^+C0
#ASSUME_CONST({
#        LABEL_IN:'FORM_MRCC_LAG_E',
#        LABEL_RES:'FORM_MRCC_LAG_E',
#        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'FORM_GAM0'})

# remove fully contracted part of C0^+C0
#ASSUME_CONST({
#        LABEL_IN:'FORM_MRCC_LAG_A1',
#        LABEL_RES:'FORM_MRCC_LAG_A1',
#        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'FORM_GAM0'})

# remove fully contracted part of C0^+C0
#ASSUME_CONST({
#        LABEL_IN:'FORM_MRCC_LAG_A2',
#        LABEL_RES:'FORM_MRCC_LAG_A2',
#        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

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

# not required any more:
#FT1SSQ = stf.Formula("F_T1SSQ:T1s2=<T1s2'*T1*T1*T1s2'>") # it seems that "avoid" is not accepted for Formula
## as a work-around, we added this dummy term ^^^^^^^
#FT1SSQ.append("<0.5*T1s2'*T1s''*T1s'''*T1s2'>", avoid=["T1s''","T1s'''"])
#FT1SSQ.set_rule()
#INVARIANT({LABEL_RES:'F_T1SSQ',LABEL_IN:'F_T1SSQ',OPERATORS:'T1',OP_RES:'T1s2'}) # remove the dummy term
#PRINT_FORMULA({LABEL:'F_T1SSQ',MODE:'SHORT'})
#
#EXPAND({LABEL_RES:"F_T1SSQ2",LABEL_IN:"F_T1SSQ",INTERM:"F_T1SUM"})
#
#PRINT_FORMULA({LABEL:'F_T1SSQ2',MODE:'SHORT'})

# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_en > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
if nc_en > 2:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
if nc_en > 3:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
#FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SSQ2'})
#EXPAND({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SSQ'})


# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_rs > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})
if nc_rs > 2:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
#FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SSQ2'})
#EXPAND({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SSQ'})


# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_rs > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
if nc_rs > 2:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
# And a last set of call to correctly replace 0.5(T1+T2g)(T1+T2g)->0.5(T1s)(T1s)
#FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SSQ2'})
#EXPAND({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SSQ'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})




#Make the Derivative with respect to LAM
# only dummy
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_RES1_0',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES1_0',LABEL_RES:'FORM_MRCC_RES1_0'})

DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_RES2_0',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES2_0',LABEL_RES:'FORM_MRCC_RES2_0'})

# Replace in singles residual as much as possible by singles part of doubles residual
# in cases where T1 and T2 are treated on the same footing, this replaces everything
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'FORM_MRCC_RES2_0'})

#PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
#PRINT_FORMULA({LABEL:'FORM_MRCC_RES2_0',MODE:'SHORT'})


K4E = True
if (K4E):
    # generate the intermediate to K4E -> INTkx
    DEF_OP_FROM_OCC({LABEL:'INTkx',JOIN:2,DESCR:',VV;PP,|,V;PP,H|,;PP,HH'})
    
    DEF_ME_LIST({LIST:'ME_INTkx',OPERATOR:'INTkx',IRREP:1,'2MS':0,AB_SYM:+1})

    # define a formal Hamiltonian that only contains the blocks of interest
    DEF_OP_FROM_OCC({LABEL:'Hpppp',DESCR:'PP;PP'}) 

    # replace relevant blocks of Hamiltonian with Hpppp
    REPLACE({LABEL_RES:'F_preINTkx',LABEL_IN:'FORM_MRCC_RES2_0',OP_LIST:['H','Hpppp']})
    # remove all terms that do not include Hpppp
    INVARIANT({LABEL_RES:'F_preINTkx',LABEL_IN:'F_preINTkx',OPERATORS:'H',OP_RES:'O2g'})
    # now, take derivative with respct to Hpppp
    DERIVATIVE({LABEL_RES:'F_INTkx',LABEL_IN:'F_preINTkx',OP_RES:'INTkx',OP_DERIV:'Hpppp'})
    REORDER_FORMULA({LABEL_IN:'F_INTkx',LABEL_RES:'F_INTkx'})

    PRINT_FORMULA({LABEL:'F_INTkx',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'F_INTkx',SPLIT:True})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'F_INTkx',SPLIT:True})
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})


### Add here more terms to be factored out ...
I3ext = True
if (I3ext):
    # try this for 3externals
    # has the shape of O2g
    CLONE_OPERATOR({LABEL:'INT3ext',TEMPLATE:'O2g'})
    DEF_ME_LIST({LIST:'ME_INT3ext',OPERATOR:'INT3ext',IRREP:1,'2MS':0,AB_SYM:+1})

    # define a formal Hamiltonian that only contains the blocks of interest
    DEF_OP_FROM_OCC({LABEL:'Hppph',DESCR:'PP;PH|PP;PV'})

    # replace relevant blocks of Hamiltonian with Hppph
    REPLACE({LABEL_RES:'F_preINT3ext',LABEL_IN:'FORM_MRCC_RES2_0',OP_LIST:['H','Hppph']})
    # remove all terms that do not include Hppph
    INVARIANT({LABEL_RES:'F_preINT3ext',LABEL_IN:'F_preINT3ext',OPERATORS:'H',OP_RES:'O2g'})
    # now, take derivative with respct to Hppph
    DERIVATIVE({LABEL_RES:'F_INT3ext',LABEL_IN:'F_preINT3ext',OP_RES:'INT3ext',OP_DERIV:'Hppph'})
    REORDER_FORMULA({LABEL_IN:'F_INT3ext',LABEL_RES:'F_INT3ext'})

    PRINT_FORMULA({LABEL:'F_INT3ext',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'F_INT3ext',SPLIT:True})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'F_INT3ext',SPLIT:True})

    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

# experimental area:
HGamma = False
if (HGamma):
    # define a number of intermediates where the Hamiltonian and the density are connected
    DEF_OP_FROM_OCC({LABEL:'GM0',JOIN:2,DESCR:',;,'})
    DEF_OP_FROM_OCC({LABEL:'GM1',JOIN:2,DESCR:',V;V,'})
    DEF_OP_FROM_OCC({LABEL:'GM2',JOIN:2,DESCR:',VV;VV,'})
    DEF_HAMILTONIAN({LABEL:'H2',MIN_RANK:2,MAX_RANK:2})
    # A: effective fock operator:
    #DEF_OP_FROM_OCC({LABEL:'INTHG0',JOIN:1,DESCR:'P,[HV];V,[HV];H,H;[HV],P;H,V'})
    DEF_OP_FROM_OCC({LABEL:'INTHG0',JOIN:1,DESCR:'P,[HV];V,[HV];H,H'})
    DEF_ME_LIST({LIST:'ME_INTHG0',OPERATOR:'INTHG0',IRREP:1,'2MS':0,AB_SYM:+1})
    FHG0 = stf.Formula("F_HG0:INTHG0=<INTHG0'*GM0'*H*GM0'*INTHG0'>")
    FHG0.append("<INTHG0'*GM1'*H*GM1'*INTHG0'>",connect=[2,3,3,4])
    FHG0.set_rule()
    REPLACE({LABEL_RES:'F_HG0',LABEL_IN:'F_HG0',OP_LIST:['GM0','GAM0','GM1','GAM0']})
    REORDER_FORMULA({LABEL_IN:'F_HG0',LABEL_RES:'F_HG0'})
    # B: H transformed with one-particle density on creation index
    DEF_OP_FROM_OCC({LABEL:'INTHG1',JOIN:2,DESCR:'[HV],[HPV][HPV];V,|[HPV],[HPV][HV];V,'})
    DEF_ME_LIST({LIST:'ME_INTHG1',OPERATOR:'INTHG1',IRREP:1,'2MS':0,AB_SYM:+1})
    FHG1 = stf.Formula("F_HG1:INTHG1=<INTHG1'*INTHG1'*H*INTHG1'*INTHG1'>") # dummy term
    FHG1.append("<INTHG1'*GM1'*H*INTHG1'*INTHG1'*GM1'*INTHG1'>",connect=[2,3,5,6])
    FHG1.set_rule()
    REPLACE({LABEL_RES:'F_HG1',LABEL_IN:'F_HG1',OP_LIST:['GM1','GAM0']})
    REORDER_FORMULA({LABEL_IN:'F_HG1',LABEL_RES:'F_HG1'})
    # C: H transformed with two-particle density on creation indices
    DEF_OP_FROM_OCC({LABEL:'INTHG2',JOIN:2,DESCR:',[HPV][HPV];VV,'})
    DEF_ME_LIST({LIST:'ME_INTHG2',OPERATOR:'INTHG2',IRREP:1,'2MS':0,AB_SYM:+1})
    FHG2 = stf.Formula("F_HG2:INTHG2=<INTHG2'*INTHG2'*H*INTHG2'*INTHG2'>") # dummy term
    FHG2.append("<INTHG2'*GM2'*H2*INTHG2'*INTHG2'*GM2'*INTHG2'>",connect=[2,3,5,6])
    FHG2.set_rule()
    REPLACE({LABEL_RES:'F_HG2',LABEL_IN:'F_HG2',OP_LIST:['GM2','GAM0','H2','H']})
    REORDER_FORMULA({LABEL_IN:'F_HG2',LABEL_RES:'F_HG2'})
    
#    DEF_OP_FROM_OCC({LABEL:'INTHGAM',JOIN:2,DESCR:',;[HPV],[HPV]|,V;[PV],|,V;,H|,V;[HPV]V,[HPV]|,V;[HPV]P,[HV]|,V;[HPV],[HPV]H|,VV;[PV][PV],|,VV;[PV],H|,VV;,HH'})
#    HGproto = stf.Formula("F_HGprotoL:MRCC_LAG="+_L2_refexp("[H,T1]"))
#    HGproto.set_rule()
#    FACTOR_OUT({
#        LABEL_IN:'F_HGprotoL',
#        LABEL_RES:'F_HGprotoL',
#        INTERM:'FORM_GAM0'})

    PRINT_FORMULA({LABEL:'F_HG0',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HG1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HG2',MODE:'SHORT'})

    
    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:['F_HG0','F_HG0^+','F_HG1','F_HG1^+','F_HG2','F_HG2^+']})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:['F_HG0','F_HG0^+','F_HG1','F_HG1^+','F_HG2','F_HG2^+']})
    
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

    #ABORT({COMMENT:'Development Stop'})


# now creating the actual residuals
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_RES1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES1',LABEL_RES:'FORM_MRCC_RES1_0'})

DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_RES2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES2',LABEL_RES:'FORM_MRCC_RES2_0'})


if (remove_gamma0):
    if (K4E):
        ASSUME_CONST({
            LABEL_IN:'F_INTkx',
            LABEL_RES:'F_INTkx',
            OP_LIST:['GAM0'],VAL_LIST:[1.0]})
    if (I3ext):
        ASSUME_CONST({
            LABEL_IN:'F_INT3ext',
            LABEL_RES:'F_INT3ext',
            OP_LIST:['GAM0'],VAL_LIST:[1.0]})
    if (HGamma):
        ASSUME_CONST({
            LABEL_IN:'F_HG0',
            LABEL_RES:'F_HG0',
            OP_LIST:['GAM0'],VAL_LIST:[1.0]})

    ASSUME_CONST({
        LABEL_IN:'FORM_MRCC_RES2',
        LABEL_RES:'FORM_MRCC_RES2',
        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

    ASSUME_CONST({
        LABEL_IN:'FORM_MRCC_RES1',
        LABEL_RES:'FORM_MRCC_RES1',
        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

    ASSUME_CONST({
        LABEL_IN:'FORM_MRCC_LAG_E',
        LABEL_RES:'FORM_MRCC_LAG_E',
        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

    
if (I3ext and not HGamma):
    _opt_label_list = ['F_T1SUM','F_INTkx','F_INT3ext','FORM_MRCC_RES2','FORM_MRCC_RES1','FORM_MRCC_LAG_E']
    _itf_code_list = ['<Sum_T1>','T1s','<Update_INTkx>','INTkx','<Residual>','INT3ext','MRCC_LAG','O1','O2g']
elif (I3ext and HGamma):
    _opt_label_list = ['F_HG0','F_HG1','F_HG2','F_T1SUM','F_INTkx','F_INT3ext','FORM_MRCC_RES2','FORM_MRCC_RES1','FORM_MRCC_LAG_E']
    _itf_code_list = ['<Make_HGAM>','INTHG0','INTHG1','INTHG2','<Sum_T1>','T1s','<Update_INTkx>','INTkx','<Residual>','INT3ext','MRCC_LAG','O1','O2g']
else:
    _opt_label_list = ['F_T1SUM','F_INTkx','FORM_MRCC_RES2','FORM_MRCC_RES1','FORM_MRCC_LAG_E']
    _itf_code_list = ['<Sum_T1>','T1s','<Update_INTkx>','INTkx','<Residual>','MRCC_LAG','O1','O2g']
    
#    OPTIMIZE({
#        LABEL_OPT:'FOPT_MRCC_LAG',
#        LABELS_IN:['F_T1SUM','F_INTkx','F_INT3ext','FORM_MRCC_RES2','FORM_MRCC_RES1','FORM_MRCC_LAG_E']})

    
OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        LABELS_IN:_opt_label_list})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_RES1',MODE:'SHORT'}) # only dummy
PRINT_FORMULA({LABEL:'FORM_MRCC_RES2',MODE:'SHORT'})

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

#    TRANSLATE_ITF({
#        LABEL:'FOPT_MRCC_LAG',
#        OUTPUT:filename,
#        TITLE:filename2,
#        MULTI:True,
#        PROCESS:True,
#        KEXT:True,
#        TASKS:tasks,
#        INIT_RES:False,
#        ITIN:True,
#        RENAME:['MRCC_LAG','ECC','T1','T1','T2g','T2','O1','R1','O2g','R2','GAM0','Ym<RANK>'],
#        CODE:['<Sum_T1>','T1s','<Update_INTkx>','INTkx','<Residual>','INT3ext','MRCC_LAG','O1','O2g']})

skip_itf = False
if not skip_itf:
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
        CODE:_itf_code_list})

#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
