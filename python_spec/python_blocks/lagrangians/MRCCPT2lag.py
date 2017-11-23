from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

i_am="MRCCPT2lag.py"

known_hamiltonians=["DYALL","DYALL-X","REPT","F_EFF","F_EFF-D"]
hamiltonian="DYALL"
if keywords.is_keyword_set('method.MRCCPT2.hamiltonian'):
    hamiltonian=str(keywords.get('method.MRCCPT2.hamiltonian')).strip()
print("hamiltonian: ", hamiltonian, type(hamiltonian))

if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

connected=True
if keywords.is_keyword_set('method.MRCCPT2.connected'):
    if (keywords.get('method.MRCCPT2.connected') == "T"):
        connected=True
    elif(keywords.get('method.MRCCPT2.connected') == "F"):
        connected=False
    else :
        raise Exception(i_am+": unrecognised value for option connected (must be T or F)")
print("connected ", connected, type(connected))

connected_O3=connected
if keywords.is_keyword_set('method.MRCCPT2.connected_O3'):
    if (keywords.get('method.MRCCPT2.connected_O3') == "T"):
        connected_O3=True
    elif(keywords.get('method.MRCCPT2.connected_O3') == "F"):
        connected_O3=False
    else :
        raise Exception(i_am+": unrecognised value for option connected_O3 (must be T or F)")
print("connected_O3 ", connected_O3, type(connected_O3))

# Third order term for the energy
third_ord_energy = False
if keywords.is_keyword_set('method.MRCCPT2.3rd_E'):
    if (keywords.get('method.MRCCPT2.3rd_E') == "T"):
        third_ord_energy=True
    elif(keywords.get('method.MRCCPT2.3rd_E') == "F"):
        third_ord_energy=False
    else :
        raise Exception(i_am+": unrecognised value for option 3rd_E (must be T or F)")
print("3rd_E ", third_ord_energy, type(third_ord_energy))

ampl_type = "PT2"
if keywords.is_keyword_set('method.MRCCPT2.ampl_type'):
    ampl_type = keywords.get('method.MRCCPT2.ampl_type')
print("ampl_type: ", ampl_type, type(ampl_type))

spinadapt=0
if keywords.is_keyword_set('calculate.routes.spinadapt'):
    spinadapt=int(keywords.get('calculate.routes.spinadapt'))

#------------------------------------------------------------------------------#
#Define the MRCCPT(2) lagrangian
#------------------------------------------------------------------------------#

new_target('DEF_FORM_PT_LAG')

# note: we use T2g as (T1,T2); this is only valid if we treat T1 and T2 on exactly the same footing
depend('DEF_T2g')
#depend('DEF_T1')

depend('DEF_LAM2g')
#depend('DEF_LAM1')

depend('DEF_O2g')
#depend('DEF_O1')

depend('MakeRefState')
depend('GAM0_CALC')
depend('H0')


if hamiltonian=="DYALL":
    depend('EVAL_HAM_D')
elif hamiltonian=="DYALL-X":
    depend('EVAL_HAM_DX')
elif hamiltonian=="REPT":
    depend('EVAL_REPT_HAM')
elif hamiltonian=="F_EFF":
    depend('EVAL_F_EFF')
elif hamiltonian=="F_EFF-D":
    depend('EVAL_F_EFF_D')

DEF_SCALAR({
        LABEL:'PT_LAG'})
DEF_ME_LIST({
        LIST:'PT_LAG_LST',
        OPERATOR:'PT_LAG',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1})

#dummy for residue of amplitude equation
DEF_SCALAR({
        LABEL:'PT_LAG_A'})
DEF_SCALAR({
        LABEL:'OVL_A'})

#some auxiliary operators and lists for the solver
CLONE_OPERATOR({LABEL:'H0xT',TEMPLATE:'O2g'})
CLONE_OPERATOR({LABEL:'SxT',TEMPLATE:'O2g'})
CLONE_OPERATOR({LABEL:'H1rhs',TEMPLATE:'O2g'})


DEF_ME_LIST({LIST:'H0xT_LST',OPERATOR:'H0xT',
             IRREP:1,'2MS':0,AB_SYM:+1,'S2':0})

DEF_ME_LIST({LIST:'SxT_LST',OPERATOR:'SxT',
             IRREP:1,'2MS':0,AB_SYM:+1,'S2':0})

DEF_ME_LIST({LIST:'H1rhs_LST',OPERATOR:'H1rhs',
             IRREP:1,'2MS':0,AB_SYM:+1,'S2':0})

if hamiltonian=="DYALL":
    _h0_='HAM_D'
    _h0exp_='HAM_D_EXP'
elif hamiltonian=="DYALL-X":
    _h0_='HAM_DX'
    _h0exp_='missing'
elif hamiltonian=="REPT":
    _h0_='REPT_HAM'
    _h0exp_='missing'
elif hamiltonian=="F_EFF":
    _h0_='FOCK_EFF'
    _h0exp_='FOCK_EFF_EXP'
elif hamiltonian=="F_EFF-D":
    _h0_='FOCK_EFF_D'
    _h0exp_='FOCK_EFF_D_EXP'

_h1_ = "(H-" +_h0_ + ")"

#Starting energy and amplitudes equations
LAG_E=stf.Formula("FORM_PT_LAG:PT_LAG="\
                  "<C0^+*H*C0>")
LAG_A=stf.Formula("FORM_PT_LAG_A:PT_LAG="\
                  "<C0^+*(LAM2g*H)*C0>")

LAG_E.append("<C0^+*H*T2g*C0>")

# How the amplitudes equations look like:
if ampl_type == 'PT2':
    if (connected):
        LAG_A.append("<C0^+*(LAM2g)*(["+_h0_+",T2g])*C0>")
    else:
        LAG_A.append("<C0^+*(LAM2g)*("+_h0_+"-"+_h0exp_+")*T2g*C0>")

elif ampl_type == 'LCC':
    LAG_A.append("<C0^+*(LAM2g)*([H,T2g])*C0>")

elif ampl_type == 'MRPT_3O-CC':
    LAG_A.append("<C0^+*(LAM2g)*([H,T2g])*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*([["+_h0_+",T2g],T2g])*C0>")

elif ampl_type == 'MRPT_3O-CI':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")

elif ampl_type == 'CC_2comm':
    LAG_A.append("<C0^+*(LAM2g)*([H,T2g])*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*([[H,T2g],T2g])*C0>")

# testing: H + (HT)c - (TH0)c.
elif ampl_type == 'HTc_TH0c':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>", connect=["H","T2g"])
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>", connect=[_h0_,"T2g"])

# testing: H + HT. Does not work. Does not converge.
elif ampl_type == 'HT_conn':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>", connect=["H","T2g"])

elif ampl_type == 'HT_HTT_conn':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>", connect=["H","T2g"])
    LAG_A.append("<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])

elif ampl_type == 'HT_all':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")

# Amplitudes equations as first order, but with Heff up to second order
# Does not work well...
elif ampl_type == 'PT1_Heff2':
    LAG_A.append("<C0^+*(LAM2g)*(["+_h0_+",T2g])*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*C0*C0^+'*H*T2g'*C0'>", avoid=["C0^+","H",
                                                                   "C0^+","T2g'",
                                                                   "C0^+","C0'",
                                                                   "LAM2g","H",
                                                                   "LAM2g","T2g'",
                                                                   "LAM2g","C0'",
                                                                   "T2g","H",
                                                                   "T2g","T2g'",
                                                                   "T2g","C0'"
                                                                   ])

else:
    quit_error(i_am+': Unknown ampl_type: ' + ampl_type)



# optional penality term
#LAG_E.append("<C0^+*(T2g^+)*O2g*C0>")

#dbg
#for item in LAG_E.show():
#    print item
LAG_E.set_rule()

#print("LAG_E finished")

#dbg
#for item in LAG_A.show():
#    print item
LAG_A.set_rule()

#print("LAG_A finished")
debug_FORM('FORM_PT_LAG_A')

#we need the overlap terms for the LEQ solver
OVL=stf.Formula("FORM_OVL:OVL_A=<C0^+*LAM2g*T2g*C0>")
OVL.set_rule()


# By now, for testing purposes.
# Replaces h0 to H.
# Is it okay for every kind of amplitude equation?
#if True:
#    REPLACE({
#        LABEL_RES:'FORM_PT_LAG_A',
#        LABEL_IN:'FORM_PT_LAG_A',
#        OP_LIST:[_h0_,'H']
#        })

FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG',
        LABEL_IN:'FORM_PT_LAG',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
        LABEL_RES:'FORM_PT_LAG_A',
        LABEL_IN:'FORM_PT_LAG_A',
        INTERM:'FORM_GAM0'})

FACTOR_OUT({
       LABEL_RES:'FORM_OVL',
       LABEL_IN:'FORM_OVL',
       INTERM:'FORM_GAM0'
       })

debug_FORM('FORM_PT_LAG_A')

mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_LAG_A'})


debug_FORM('FORM_PT_LAG_A')



#Make the Derivative with respect to LAM2g.
DERIVATIVE({LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_Amp',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})

DERIVATIVE({LABEL_IN:'FORM_OVL',
        LABEL_RES:'FORM_PT_SxT',
        OP_RES:'SxT',
        OP_DERIV:'LAM2g'})

debug_FORM('FORM_PT_Amp')
debug_FORM('FORM_PT_SxT')

#devide into RHS and transformation part
LEQ_SPLIT({LABEL_RAW:'FORM_PT_Amp',
           LABEL_TRF:'FORM_PT_MVP',
           LABEL_RHS:'FORM_PT_RHS',
           OP_TRF:'H0xT',
           OP_RHS:'H1rhs',
           OP_X:'T2g'})

debug_FORM('FORM_PT_MVP')
debug_FORM('FORM_PT_RHS')

TEX_FORMULA({LABEL:'FORM_PT_LAG_A',OUTPUT:'PT2-LAG-A.tex'})

OPTIMIZE({
        LABEL_OPT:'FOPT_PT_LAG',
        LABELS_IN:['FORM_PT_Amp','FORM_PT_LAG']})

#OPTIMIZE({
#        LABEL_OPT:'FOPT_PT_LAG',
#        LABELS_IN:['FORM_PT_LAG']})

OPTIMIZE({
        LABEL_OPT:'FOPT_PT_EQ',
        LABELS_IN:['FORM_PT_RHS','FORM_PT_MVP','FORM_PT_SxT']})


# Comments from master
# We need to check the third-order energy expression again
# (e.g. replace H by H_N (i.e. without scalar contribution) for formula generation
#  and replace back H_N->H aferwards (using "REPLACE")

# Third order terms
term_CI = "<C0^+ * (T2g^+) * (" + _h1_ + "*T2g) * C0>"
term_CC0 =       "<C0^+ * (T2g^+) * (["  + _h1_ + ",T2g]) * C0>"
term_CCa = "1/2 * <C0^+ *           ([[" + _h1_ + ",T2g],T2g]) * C0>"
term_CCb = "1/2 * <C0^+ * (T2g^+) * ([[" + _h0_ + ",T2g],T2g]) * C0>"

# The CC like third order term, written an another way
term_pureV = "-<C0^+ * (T2g^+) * T2g * ( " +_h1_+ " + [" +_h0_+ ",T2g])  * C0>"
term_highO = "1/2 * <C0^+ * ( " +_h1_+ " + [(T2g^+)," +_h0_+ "]) * T2g * T2g * C0>"

term_2Oa = "<C0^+ * (T2g^+) * H * C0>"
term_2Ob = "<C0^+ * (T2g^+) * [" + _h0_ + ",T2g] * C0>"

if (third_ord_energy):
    new_target('MRCCPT_E_3rd_O', True)
    heading('Third order correction for the energy')
    depend('SOLVE_MRCCPT2')

    for i in ['_CI', '_CC0', '_CCa', '_CCb', '_CC', '_CC_pureV', '_CC_higherO']:
        if (i == '_CI'):
            term = term_CI
            str_i = "Rayleigh-Schrodinger like term: <T^+ W T>"

        elif (i == '_CC0'):
            term = term_CC0
            str_i = "term 0 of CC like: <T^+ [H1,T]>"

        elif (i == '_CCa'):
            term = term_CCa
            str_i = "term a of CC like: 1/2 <[[H1,T],T]>"

        elif (i == '_CCb'):
            term = term_CCb
            str_i = "term b of CC like: 1/2 <T^+ [[H0,T],T]>"

        elif (i == '_CC'):
            term = [term_CC0,
                    term_CCa,
                    term_CCb]
            str_i = "CC like term"

        elif (i == '_CC_pureV'):
            term = term_pureV
            str_i = "CC term with pure virtual excitations"

        elif (i == '_CC_higherO'):
            term = term_highO
            str_i = "CC term with higher order excitations"

        else:
            raise Exception(i_am + ": unrecognised kind of third order correction: " + i)

# This includes the terms from the 2nd order Lagrangian
# that are exactly zero if ampl_type == 'PT2'
#        if isinstance(term, list):
#            term.extend([term_2Oa, term_2Ob])
#        else:
#            term = [term, term_2Oa, term_2Ob]

        DEF_SCALAR({LABEL:'MRCCPT_O3'+i})
        DEF_ME_LIST({LIST:'ME_MRCCPT_O3'+i,
                     OPERATOR:'MRCCPT_O3'+i,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        if (isinstance(term, str)):
            E_O3 = stf.Formula("FORM_MRCCPT_O3"+i+":MRCCPT_O3"+i+"="+term)
        else:
            E_O3 = stf.Formula("FORM_MRCCPT_O3"+i+":MRCCPT_O3"+i+"="+term[0])
            for j in range(1, len(term)):
                E_O3.append(term[j])

        E_O3.set_rule()

        OPTIMIZE({LABEL_OPT:'FOPT_MRCCPT_O3'+i,
                  LABELS_IN:['FORM_MRCCPT_O3'+i]})

        EVALUATE({FORM:'FOPT_MRCCPT_O3'+i})

        DEF_SCALAR({LABEL:'E_MRCCPT2_plus_O3'+i})
        DEF_ME_LIST({LIST:'ME_E_MRCCPT2_plus_O3'+i,
                     OPERATOR:'E_MRCCPT2_plus_O3'+i,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        ASSIGN_ME2OP({LIST:'PT_LAG_LST',OPERATOR:'PT_LAG'})

        E_plus_3rd = stf.Formula("F_MRCCPT2_plus_O3"+i+":E_MRCCPT2_plus_O3"+i+"=<PT_LAG>+<MRCCPT_O3"+i+">")
        E_plus_3rd.set_rule()

        OPTIMIZE({LABEL_OPT:'FOPT_MRCCPT2_plus_O3'+i,
                  LABELS_IN:['F_MRCCPT2_plus_O3'+i]})

        EVALUATE({FORM:'FOPT_MRCCPT2_plus_O3'+i})

        PRINT_MEL({
                LIST:'ME_MRCCPT_O3'+i,
                COMMENT:'3rd order correction          ('+str_i+')',
                FORMAT:'SCAL F24.14'})
        PRINT_MEL({
                LIST:'ME_E_MRCCPT2_plus_O3'+i,
                COMMENT:'MRCCPT3 = MRCCPT2 + 3rd order ('+str_i+')',
                FORMAT:'SCAL F24.14'})
        PUSH_RESULT({
                LIST:'ME_E_MRCCPT2_plus_O3'+i,
                COMMENT:"icMRPT3"+i,
                FORMAT:"SCAL F24.14"})
