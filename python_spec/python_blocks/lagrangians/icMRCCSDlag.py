"""An implementation of the icMRCCSD Lagrangian, done term by term and separated T1 and T2.

History:

Yuri august 2017: Creation based on MRCC2lag.py. Implementation up to maxcom=2

"""

from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf
import python_blocks.lagrangians.ref_relaxation as ref_relaxation
import python_blocks.lagrangians.mrcc_methods as mrcc_methods


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

# hybrid approximation?
approx = keywords.get('method.MR_P.hybrid')
hybrid = approx if approx is not None else "none"

word = keywords.get('method.MRCC.maxcom_en')
nc_en = int(word) if word is not None else 4
word = keywords.get('method.MRCC.maxcom_res')
nc_rs = int(word) if word is not None else 2

word = keywords.get('general.print')
verbosity = int(word) if word is not None else 0

word = keywords.get('method.MRCC.select')
#print("word = "+str(word))
if word is None:
    select =  True
else:
    if word == "F":
        select = False
    elif word == "T":
        select = True
    else:
        quit_error('select must be T or F, found: '+word)

itfgen = keywords.is_keyword_set('method.ITF')

doublet=False
cas22=False
orbinfo = Orb_Info()
nactel = orbinfo.get('nactel')
nactorb = orbinfo.get('nactorb')
nocc_el = orbinfo.get('nactt_hpv',1)
if (nactel==1 and nactorb==1):
    doublet=True
elif (nactel==2 and nactorb==2):
    cas22=True
if (itfgen and not doublet and not cas22):
    quit_error('ITFgen called for case not considered yet!')

if (hybrid!="none"):
    no_occ=False
    if keywords.is_keyword_set('method.MR_P.no_occ'):
        print("Use of no_occ is deprecated!")
        string = (keywords.get('method.MR_P.no_occ'))
        if string == 'T':
            no_occ = True
        elif string == 'F':
            no_occ = False
        else:
            print("Didn't recognise no_occ argument; setting to False")
            no_occ = False
        print("Are there no occupied orbitals: ", no_occ)

    if (no_occ and nocc_el>0):
        print("WARNING: You set no_occ, but occupied orbitals are present!")

    no_occ = no_occ or nocc_el==0

    string = keywords.get('method.MR_P.separation')
    separation = string if string is not None else "SY"  # standard: Saitow-Yanai separation

    if keywords.is_keyword_set('method.MR_P.singles'): #remove single keyword or CEPA(0) and PT2 calculations will crash accordingly
       singles = int((keywords.get('method.MR_P.singles')))
       if singles == 0:
            print("All singles operators are in the internal space")
       else:
            print("Singles operators are split between the internal and external spaces")
    else:
       if hybrid in ["CEPT2","CCEPA","TCPT2"]:
           singles = 0
       else:
           singles = 8
       print("All singles operators are in the internal space")

    known_methods=["CEPT2","CCEPA","TCPT2","CEPA0","PT2"]
    if hybrid not in known_methods :
        raise Exception(i_am+": unknown method:"+str(hybrid))
    print("Using the special "+str(hybrid)+" method.")

    if hybrid in ["CEPT2","CCEPA","TCPT2"]:
        if separation == "SY":
            print("Saitow-Yanai separation")
        if separation == "SY-INV":
            print("inverted Saitow-Yanai separation")
            singles = 7 # force this
        # maybe:
        #if separation == "WCSK":
        #    print("Werner-Celani-Stoll-Knowles separation")
    
# Select H0
    ham = ''
    hamiltonian = ''
    if hybrid in ['CEPT2','PT2','TCPT2']:
        known_hamiltonians=["DYALL","REPT","F_EFF"]
        hamiltonian = keywords.get('method.MR_P.H0')
        hamiltonian=str(hamiltonian).strip() if hamiltonian is not None else "DYALL"

        if hamiltonian not in known_hamiltonians :
            raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))


else:
    print("Settings for MRCC generator:")
    print(("maxcom_en  = "+str(nc_en)))
    print(("maxcom_res = "+str(nc_rs)))
    print(("select     = "+str(select)))

if (doublet):
    print("detected CAS(1,1) case")
if (cas22):
    print("detected CAS(2,2) case")


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
remove_gamma0 = True # remove the scalar part of GAM0 (which is just 1.0)

# make sure that GAM0 is recognized as a Hermitian operator (for transpose):
SET_HERMITIAN({LABEL:'GAM0',CA_SYMMETRY:+1})

# set requested method
if (hybrid=="none"):
    maxexc=2
    mrcc_methods.set_mrcc(maxexc,nc_en,nc_rs,select,(doublet or cas22))
else:
    mrcc_methods.set_hybrids(hybrid,separation,hamiltonian,singles,no_occ,(doublet or cas22))

if verbosity >= 100:
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'COUNT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'COUNT'})
PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'COUNT'})

if hybrid in ['CEPT2','CCEPA','CEPA0']:
       # Construct energy operator for use in lagrangian
   DEF_ME_LIST({LIST:'ME_CEPA',
                OPERATOR:'ECEPA',
                IRREP:1,
                '2MS':0,
                AB_SYM:+1})
    

# --- factor out densities ---
# do this later, as we need it for reference relaxation form
#FACTOR_OUT({
#        LABEL_IN:'FORM_MRCC_LAG_E',
#        LABEL_RES:'FORM_MRCC_LAG_E',
#        INTERM:'FORM_GAM0'})

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

if verbosity >= 100:
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

#Make the Derivative with respect to LAM
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_RES2_0',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES2_0',LABEL_RES:'FORM_MRCC_RES2_0'})

# Replace in singles residual as much as possible by singles part of doubles residual
# in cases where T1 and T2 are treated on the same footing, this replaces everything
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'FORM_MRCC_RES2_0'})


# define an operator that sums the two singles operators
if doublet or cas22:
    T1s_shape = 'P,H'
    DEF_OP_FROM_OCC({LABEL:'T1s',DESCR:T1s_shape})
else:
    CLONE_OPERATOR({LABEL:'T1s',TEMPLATE:'T1'})

CLONE_OPERATOR({LABEL:'T1s2',TEMPLATE:'T2g'}) # for squared T1s

DEF_ME_LIST({LIST:'ME_T1s',OPERATOR:'T1s',IRREP:1,'2MS':0,AB_SYM:+1})

DEF_FORMULA({LABEL:'F_T1SUM',FORMULA:'T1s=T1+T2g'})

# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_en > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
if nc_en > 2:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})
if nc_en > 3:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_E',LABEL_RES:'FORM_MRCC_LAG_E',INTERM:'F_T1SUM',SPLIT:True})

# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_rs > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})
if nc_rs > 2 and not select:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A1',LABEL_RES:'FORM_MRCC_LAG_A1',INTERM:'F_T1SUM',SPLIT:True})

# Factor out linear terms
FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
# One more call for quadratic terms
if nc_rs > 1:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
if nc_rs > 2 and not select:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})
if nc_rs > 3 and not select:
    FACTOR_OUT({LABEL_IN:'FORM_MRCC_LAG_A2',LABEL_RES:'FORM_MRCC_LAG_A2',INTERM:'F_T1SUM',SPLIT:True})

if verbosity >= 100:
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})



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

    if verbosity >= 10:
        PRINT_FORMULA({LABEL:'F_INTkx',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'F_INTkx',SPLIT:True})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'F_INTkx',SPLIT:True})

    if verbosity >= 100:
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

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

    if (verbosity >= 10):
        PRINT_FORMULA({LABEL:'F_INT3ext',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'F_INT3ext',SPLIT:True})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'F_INT3ext',SPLIT:True})

    if (verbosity >= 100):
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
    # B-1: H transformed with one-particle hole density on creation index
    if (nc_rs<2):
        DEF_OP_FROM_OCC({LABEL:'INTHE1',JOIN:2,DESCR:'[HPV][HPV],VV;,|[HPV]H,HV;,|[HV]H,PV;,'})
        #DEF_OP_FROM_OCC({LABEL:'INTHE1',JOIN:2,DESCR:'HV,VV;,'})
    else:
        DEF_OP_FROM_OCC({LABEL:'INTHE1',JOIN:2,DESCR:'[HPV][HPV],[HV]V;,|[HPV][HV],[HPV]V;,'})
    DEF_ME_LIST({LIST:'ME_INTHE1',OPERATOR:'INTHE1',IRREP:1,'2MS':0,AB_SYM:+1})
    FHE1 = stf.Formula("F_HE1:INTHE1=<INTHE1'*GM0'*H*GM0'*INTHE1'*INTHE1'*INTHE1'>") # first term
    FHE1.append("<INTHE1'*GM1'*H*GM1'*INTHE1'*INTHE1'*INTHE1'>",connect=[2,5,3,4])
    FHE1.set_rule()
    PRINT_FORMULA({LABEL:'F_HE1',MODE:'SHORT'})
    REPLACE({LABEL_RES:'F_HE1',LABEL_IN:'F_HE1',OP_LIST:['GM0','GAM0','GM1','GAM0']})
    REORDER_FORMULA({LABEL_IN:'F_HE1',LABEL_RES:'F_HE1'})
    # C: H transformed with two-particle density on creation indices
    DEF_OP_FROM_OCC({LABEL:'INTHG2',JOIN:2,DESCR:',[HPV][HPV];VV,'})
    DEF_ME_LIST({LIST:'ME_INTHG2',OPERATOR:'INTHG2',IRREP:1,'2MS':0,AB_SYM:+1})
    FHG2 = stf.Formula("F_HG2:INTHG2=<INTHG2'*INTHG2'*H*INTHG2'*INTHG2'>") # dummy term
    FHG2.append("<INTHG2'*GM2'*H2*INTHG2'*INTHG2'*GM2'*INTHG2'>",connect=[2,3,5,6])
    FHG2.set_rule()
    REPLACE({LABEL_RES:'F_HG2',LABEL_IN:'F_HG2',OP_LIST:['GM2','GAM0','H2','H']})
    REORDER_FORMULA({LABEL_IN:'F_HG2',LABEL_RES:'F_HG2'})
    # C-1: H transformed with two-particle hole-density on creation indices
    DEF_OP_FROM_OCC({LABEL:'INTHE2',JOIN:2,DESCR:'[HPV][HPV],VV;,'})
    DEF_ME_LIST({LIST:'ME_INTHE2',OPERATOR:'INTHE2',IRREP:1,'2MS':0,AB_SYM:+1})
    FHE2 = stf.Formula("F_HE2:INTHE2=<INTHE2'*GM0'*H*GM0'*INTHE2'*INTHE2'*INTHE2'>") # first term
    FHE2.append("<INTHE2'*GM1'*H2*GM1'*INTHE2'*INTHE2'*INTHE2'>",connect=[2,5,3,4])
    FHE2.append("<INTHE2'*GM2'*H2*GM2'*INTHE2'*INTHE2'*INTHE2'>",connect=[2,5,3,4],avoid=[2,4,1,4])
    FHE2.set_rule()
    REPLACE({LABEL_RES:'F_HE2',LABEL_IN:'F_HE2',OP_LIST:['GM0','GAM0','GM1','GAM0','GM2','GAM0','H2','H']})
    REORDER_FORMULA({LABEL_IN:'F_HE2',LABEL_RES:'F_HE2'})
    
    PRINT_FORMULA({LABEL:'F_HG0',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HG1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HE1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HG2',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'F_HE2',MODE:'SHORT'})

    
    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
#        INTERM:['F_HG0','F_HG0^+'],
        INTERM:['F_HG0','F_HG0^+','F_HE2','F_HE2^+']})
#        INTERM:['F_HG0','F_HG0^+','F_HG1','F_HG1^+','F_HG2','F_HG2^+']})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
#        INTERM:['F_HG0','F_HG0^+']})
        INTERM:['F_HG0','F_HG0^+','F_HE2','F_HE2^+']})
#        INTERM:['F_HG0','F_HG0^+','F_HG1','F_HG1^+','F_HG2','F_HG2^+']})
    
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

HTint = False
if (HTint):

    # this was my starting point:
    #DEF_OP_FROM_OCC({LABEL:'INTHT',JOIN:3,DESCR:',V;V,V;V,'})
    # slightly extended
#    DEF_OP_FROM_OCC({LABEL:'INTHT',JOIN:3,DESCR:',V;H,H;V,|,V;V,V;V,|,V;P,P;V,'})  # AA: 108 terms 
#    DEF_OP_FROM_OCC({LABEL:'INTHT',JOIN:3,DESCR:',V;[HPV],[HPV];V,'})  but only HH, PP, VV generate
    DEF_OP_FROM_OCC({LABEL:'INTHT',JOIN:3,DESCR:',;P,VH;V,'}) # BA: this is the only block
#    DEF_OP_FROM_OCC({LABEL:'Ltgt',JOIN:1,DESCR:'HV,PV'})  # A.
    DEF_OP_FROM_OCC({LABEL:'Ltgt',JOIN:1,DESCR:'HH,PP'})  # B.
    DEF_OP_FROM_OCC({LABEL:'Ttgt',JOIN:1,DESCR:'PV,HV'})  # .A

    # mask the target vertices
    REPLACE({LABEL_RES:'F_preHT0',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2g','Ttgt','LAM2g','Ltgt']})

    # IMPORTANT: our starting point is not the formula with all replacements made so far:
    DERIVATIVE({
        LABEL_IN:'F_preHT0',
        LABEL_RES:'F_preHT1', # <<--- presently, this has to be a new name!
        OP_RES:'O2g',
        OP_DERIV:'Ltgt'})

    REORDER_FORMULA({LABEL_IN:'F_preHT1',LABEL_RES:'F_preHT1'})
    
    DERIVATIVE({LABEL_RES:'F_INTHT',LABEL_IN:'F_preHT1',OP_RES:'INTHT',OP_DERIV:'Ttgt'})
    REPLACE({LABEL_RES:'F_INTHT',LABEL_IN:'F_INTHT',OP_LIST:['Ttgt','T2g']})

    # Very important: REORDER
    REORDER_FORMULA({LABEL_IN:'F_INTHT',LABEL_RES:'F_INTHT'})

    PRINT_FORMULA({LABEL:'F_INTHT',MODE:'SHORT'})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_A1',
        INTERM:'F_INTHT',SPLIT:True})

    FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_A2',
        INTERM:'F_INTHT',SPLIT:True})
    
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2',MODE:'SHORT'})

    # test
    Ftest = stf.Formula("F_test:MRCC_LAG=<INTHT'*Ltgt*INTHT'*Ttgt*INTHT'>")
    Ftest.set_rule()
    PRINT_FORMULA({LABEL:'F_test',MODE:'SHORT'})

    EXPAND({LABEL_IN:'F_test',LABEL_RES:'F_test_2',INTERM:'F_INTHT'})
    REPLACE({LABEL_RES:'F_test_2',LABEL_IN:'F_test_2',OP_LIST:['Ttgt','T2g','Ltgt','LAM2g']})
    SUM_TERMS({LABEL_RES:'F_test_2',LABEL_IN:'F_test_2'})
    
    PRINT_FORMULA({LABEL:'F_test_2',MODE:'SHORT'})

    SELECT_TERMS({LABEL_RES:'F_test_3',LABEL_IN:'F_preHT0',OP_RES:'MRCC_LAG',OP_INCL:['Ltgt','Ttgt'],BLK_INCL:[1,1]})
    REPLACE({LABEL_RES:'F_test_3',LABEL_IN:'F_test_3',OP_LIST:['Ttgt','T2g','Ltgt','LAM2g']})

    PRINT_FORMULA({LABEL:'F_test_3',MODE:'SHORT'})

    CONCAT({LABEL_RES:'F_test_4',LABEL_IN:['F_test_3','F_test_2'],FAC:[1.0,-1.0]})
    PRINT_FORMULA({LABEL:'F_test_4',MODE:'SHORT'})

    SUM_TERMS({LABEL_RES:'F_test_4',LABEL_IN:'F_test_4'})

    PRINT_FORMULA({LABEL:'F_test_4',MODE:'SHORT'})
    
    ABORT({COMMENT:'Development Stop'})


# now creating the actual residuals
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_RES1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES1',LABEL_RES:'FORM_MRCC_RES1'})

DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_RES2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

REORDER_FORMULA({LABEL_IN:'FORM_MRCC_RES2',LABEL_RES:'FORM_MRCC_RES2'})

# factor now GAM0 in E but keep copy for reference relaxation
FACTOR_OUT({
        LABEL_IN:'FORM_MRCC_LAG_E',
        LABEL_RES:'FORM_MRCC_LAG_ENGY',
        INTERM:'FORM_GAM0'})

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
            LABEL_IN:'F_HE1',
            LABEL_RES:'F_HE1',
            OP_LIST:['GAM0'],VAL_LIST:[1.0]})
        ASSUME_CONST({
            LABEL_IN:'F_HE2',
            LABEL_RES:'F_HE2',
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
        LABEL_IN:'FORM_MRCC_LAG_ENGY',
        LABEL_RES:'FORM_MRCC_LAG_ENGY',
        OP_LIST:['GAM0'],VAL_LIST:[1.0]})

_opt_label_list = []
_itf_code_list = []

if (HGamma):
    _opt_label_list.append('F_HG0')
    #_opt_label_list.append('F_HG1')
    #_opt_label_list.append('F_HG2')
    #_opt_label_list.append('F_HE1')
    _opt_label_list.append('F_HE2')
    _itf_code_list.append('<Make_HGAM>')
    _itf_code_list.append('INTHG0')
    _itf_code_list.append('INTHG1')
    _itf_code_list.append('INTHG2')
    _itf_code_list.append('INTHE1')
    _itf_code_list.append('INTHE2')

if hybrid in ['CEPT2','CCEPA','CEPA0']:
    _opt_label_list.append('FORM_ECEPA')
    
_opt_label_list.append('F_T1SUM')
_opt_label_list.append('FORM_GAM0')
_opt_label_list.append('F_INTkx')
_itf_code_list.append('<Sum_T1>')
_itf_code_list.append('T1s')
_itf_code_list.append('<Update_INTkx>')
_itf_code_list.append('INTkx')
_itf_code_list.append('<Residual>')

if (I3ext):
    _opt_label_list.append('F_INT3ext')
    _itf_code_list.append('INT3ext')

_opt_label_list.append('FORM_MRCC_RES2')
_opt_label_list.append('FORM_MRCC_RES1')
_opt_label_list.append('FORM_MRCC_LAG_ENGY')
_itf_code_list.append('MRCC_LAG')
_itf_code_list.append('O1')
_itf_code_list.append('O2g')
    
    
OPTIMIZE({
        LABEL_OPT:'FOPT_MRCC_LAG',
        LABELS_IN:_opt_label_list})

if verbosity >= 50:
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_ENGY',MODE:'SHORT'})
    PRINT_FORMULA({LABEL:'FORM_MRCC_RES1',MODE:'SHORT'}) # only dummy
    PRINT_FORMULA({LABEL:'FORM_MRCC_RES2',MODE:'SHORT'})

if verbosity >= 1000:
    PRINT_FORMULA({LABEL:'FOPT_MRCC_LAG',OUTPUT:'FOPT.out'})

PRINT_FORMULA({LABEL:'FOPT_MRCC_LAG',MODE:'COUNT'})

if itfgen:
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
        CODE:_itf_code_list})

#-----
# this must be the equation that still contains C0
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
