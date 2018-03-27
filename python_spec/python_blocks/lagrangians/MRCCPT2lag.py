from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf


i_am="MRCCPT2lag.py"

known_hamiltonians=["DYALL","DYALL-X","REPT","F_EFF","F_EFF-D","EXT_DYALL","SIMP_REPT"]
hamiltonian="DYALL"
if keywords.is_keyword_set('method.MRCCPT2.hamiltonian'):
    hamiltonian=str(keywords.get('method.MRCCPT2.hamiltonian')).strip()
print("hamiltonian: ", hamiltonian, type(hamiltonian))

if hamiltonian not in known_hamiltonians : 
    raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

# Option 'connected':
#
# If connected = False, triggers the equations assuming
# a zeroth order Hamiltonian of the kind
#  H^(0) = P hamiltonian P + Q hamiltonian Q
# where P and Q are the projectors into the reference and
# external spaces. This is as CASPT theory and it is needed if
# the CASCI functions are not eigenfunctions of hamiltonian
#
# If connected = True, just aplies the equations from
# perturbation theory. This will work well
# as long that the CASCI wave functions are eigenfunctions
# of hamiltonian with the same eigenvalues as the CASCI energies
# (and thus E^(1) = 0)
#
# By default, uses the appropriate version according to the
# hamiltonian
if hamiltonian in ["F_EFF","F_EFF-D"]:
    connected = False
else:
    connected = True
if keywords.is_keyword_set('method.MRCCPT2.connected'):
    if (keywords.get('method.MRCCPT2.connected') == "T"):
        connected=True
    elif(keywords.get('method.MRCCPT2.connected') == "F"):
        connected=False
    else :
        raise Exception(i_am+": unrecognised value for option connected (must be T or F)")
print("connected ", connected, type(connected))

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


# Evaluate specific terms
test_terms = False
if keywords.is_keyword_set('method.MRCCPT2.test_terms'):
    if (keywords.get('method.MRCCPT2.test_terms') == "T"):
        test_terms=True
    elif(keywords.get('method.MRCCPT2.test_terms') == "F"):
        test_terms=False
    else :
        raise Exception(i_am+": unrecognised value for option test_terms (must be T or F)")
print("test_terms ", test_terms, type(test_terms))


# Calculate 'stabm' roots of the stability matrix (Jacobian)
stabm = 0
if keywords.is_keyword_set('method.MRCCPT2.stabm'):
    stabm = keywords.get('method.MRCCPT2.stabm')
print("Stability matrix roots: ", stabm, type(stabm))


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
elif hamiltonian=="EXT_DYALL":
    depend('EVAL_HAM_EXT_D')
elif hamiltonian=="SIMP_REPT":
    depend('EVAL_SIMP_REPT_HAM')

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
elif hamiltonian=="EXT_DYALL":
    _h0_='HAM_EXT_D'
    _h0exp_='missing'
elif hamiltonian=="SIMP_REPT":
    _h0_='SIMP_REPT_HAM'
    _h0exp_='missing'

if not(connected) and _h0exp_ == 'missing':
    raise Exception(i_am+": The use of this hamiltonian: " + str(hamiltonian) + " is not compatible to not connected=F.")

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

elif ampl_type == 'CEPA-like':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like-h0':
    # CCSD-like, but with H0
    # Results in terrible energy
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like-2h0':
    # CCSD-like, but with H0 in both terms
    # Results in terrible energy
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like-3h0':
    # CCSD-like, but assum V*T*T is small
    # Better than above, but suffers from Lcc discontinuity
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*"+_h0_+"*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like-4h0':
    # Like above, but using CEPA approx for linear terms
    # Reproduces LCC discontinuity...
    # This implies that some quadratic terms are responsible for this...
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*"+_h0_+"*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-E0)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like-5h0':
    # Like above, but use full E_c
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-E0)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*("+_h0_+"-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like_2':
    # CCSD-like, but only -E_c in quadratic term (call this CCSD-like_2)
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-E0)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g*T2g*C0>")

elif ampl_type == 'CCSD-like_c':
    # CCSD-like, but only connected terms
    # With just single refernce terms - results in terrible energy
    # Try with SR and MR terms...
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>", connect=["H","T2g"])
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>", connect=["H","T2g"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g'*T2g''*C0>", connect=["H","T2g'", "H","T2g''"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])

elif ampl_type == 'CISD-like':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")

elif ampl_type == 'CEPAn-like':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")

elif ampl_type == 'HT_HTTc':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])

elif ampl_type == 'HT_HTT2':
    LAG_A.append("<C0^+*(LAM2g)*[H,T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])

elif ampl_type == '3O-SCF':
    # Third order SCF method
    LAG_E.append("1/2*<C0^+*[[(H-"+_h0_+"),T2g],T2g]*C0>")
    LAG_A.append("<C0^+*(LAM2g)*[H,T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")

elif ampl_type == '3O-SCF-CEPA':
    # Third order SCF method
    LAG_E.append("1/2*<C0^+*[[(H-"+_h0_+"),T2g],T2g]*C0>")
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")

elif ampl_type == 'no_Ec':
    # CCSD-like, but with no -Ec contribution
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g*T2g*C0>")

# Adding quadratic terms to CCSD-like
elif ampl_type == 'CCSD-like+1':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*T2g*C0>")

elif ampl_type == 'CCSD-like+2':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-ECEPA)*T2g*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*T2g*T2g*"+_h0_+"*C0>")

elif ampl_type == 'CCSD-like+3':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.append("<C0^+*H*T2g*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[[(H-ECEPA),T2g],T2g]*C0>")

elif ampl_type == 'CEPA+quad':
    # Promising!
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    
    LAG_A.append("<C0^+*(LAM2g)*(H-ECEPA)*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[[(H-ECEPA),T2g],T2g]*C0>")

elif ampl_type == 'IDEA1_2':
    # Add disconnected terms to connected cepa
    # Reproduces MRPT_3O-CI
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")

elif ampl_type == 'IDEA1_3':
    # Add/subrtract connected terms to connected cepa
    # Poor convergence
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")

elif ampl_type == 'CEPA_2c':
    # CEPA + 2 connected terms from LCC
    # This produces errors like LCC for N2
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")

elif ampl_type == 'LCC_2c':
    # LCC - 2 connected terms from LCC
    # This reproduces IDEA1, ie. connected CEPA
    LAG_A.append("<C0^+*(LAM2g)*[H,T2g]*C0>")

elif ampl_type == 'HT_HTT':
    # Simply HT + 1/2*HTT
    # Doesn't converge
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g*T2g*C0>")


elif ampl_type == 'IDEA1':
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")

elif ampl_type == 'IDEA2':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*(H-"+_h0_+")*T2g*C0>")

elif ampl_type == 'IDEA3':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g*T2g*C0>")
    LAG_A.append("-1/2*<C0^+*(LAM2g)*T2g*H*T2g*C0>")
    LAG_A.append("-1/2*<C0^+*(LAM2g)*T2g*"+_h0_+"*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*T2g*T2g*"+_h0_+"*C0>")

elif ampl_type == 'IDEA4':
    LAG_A.append("<C0^+*(LAM2g)*[H,T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*H*T2g*C0>")

elif ampl_type == 'IDEA5':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*([["+_h0_+",T2g],T2g])*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h1_+"*T2g*C0>")

elif ampl_type == 'IDEA6':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*([["+_h0_+",T2g],T2g])*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h1_+"*T2g*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*"+_h1_+"*T2g*T2g*C0>")

elif ampl_type == 'IDEA7':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)"+_h0_+"*T2g*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*T2g*C0>")

elif ampl_type == 'IDEA8':
    LAG_A.append("<C0^+*(LAM2g)*H*T2g*C0>")
    LAG_A.append("-<C0^+*(LAM2g)*T2g*"+_h0_+"*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g*T2g*C0>")
    LAG_A.append("-1/2*<C0^+*(LAM2g)*T2g*T2g*"+_h0_+"*C0>")

elif ampl_type == 'IDEA9':
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")

elif ampl_type == 'IDEA10':
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")

elif ampl_type == 'IDEA11':
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")

elif ampl_type == 'IDEA12':
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("-1/2*<C0^+*(LAM2g)*"+_h0_+"*T2g'*T2g''*C0>", connect=[_h0_,"T2g'",  _h0_,"T2g''"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])
    LAG_A.append("-1/2*<C0^+*(LAM2g)*"+_h0_+"*T2g'*T2g''*C0>", avoid=[_h0_,"T2g''"], connect=[_h0_,"T2g'",  "T2g'","T2g''"])

elif ampl_type == 'IDEA13':
    LAG_A.append("<C0^+*(LAM2g)*(H-"+_h0_+")*T2g*C0>", connect=["H","T2g",  _h0_,"T2g"])
    LAG_A.append("<C0^+*(LAM2g)*["+_h0_+",T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*[["+_h0_+",T2g],T2g]*C0>")
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("-1/2*<C0^+*(LAM2g)*"+_h0_+"*T2g'*T2g''*C0>", connect=[_h0_,"T2g'",  _h0_,"T2g''"])
    LAG_A.append("1/2*<C0^+*(LAM2g)*H*T2g'*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])
    LAG_A.append("-1/2*<C0^+*(LAM2g)*"+_h0_+"*T2g'*T2g''*C0>", avoid=[_h0_,"T2g''"], connect=[_h0_,"T2g'",  "T2g'","T2g''"])
    LAG_A.append("-<C0^+*(LAM2g)*T2g'*H*T2g''*C0>", connect=["H","T2g'",  "H","T2g''"])
    LAG_A.append("<C0^+*(LAM2g)*T2g'*"+_h0_+"*T2g''*C0>", connect=[_h0_,"T2g'",  _h0_,"T2g''"])
    LAG_A.append("-<C0^+*(LAM2g)*T2g'*H*T2g''*C0>", avoid=["H","T2g''"], connect=["H","T2g'",  "T2g'","T2g''"])
    LAG_A.append("<C0^+*(LAM2g)*T2g'*"+_h0_+"*T2g''*C0>", avoid=[_h0_,"T2g''"], connect=[_h0_,"T2g'",  "T2g'","T2g''"])


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


# Option to add specific terms using EXPAND_OP_PRODUCT to lagragian
if ampl_type == 'IDEA1_2':
    # Add disconnected terms onto FORM_PT_LAG_A
    LAG_A.set_rule()

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,VP,VV"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,VP,VV"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,P,V"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,P,V"],
                       AVOID:[3,4]})

    PRINT_FORMULA({LABEL:'FORM_PT_LAG_A',MODE:'SHORT'})

elif ampl_type == 'IDEA1_3':
    # Add connected terms onto FORM_PT_LAG_A
    LAG_A.set_rule()

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,V,P","3,,VV,VP","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,V,P","3,,V,P","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,VV,PV","3,,VV,VP","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',
                       OPERATORS:['C0^+','LAM2g','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,VV,PV","3,,V,P","4,,PP,VV"],
                       AVOID:[2,3]})

elif ampl_type == 'CEPA_2c':
    # Add connected terms to CEPA from LCC
    # These make CEPA match the LCC behaviour for N2
    LAG_A.set_rule()
    PRINT_FORMULA({LABEL:'FORM_PT_LAG_A',MODE:'SHORT'})

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',FAC:-1,
                       OPERATORS:['C0^+','LAM2g','T2g','H','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[4,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,VP,VV"]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',FAC:-1,
                       OPERATORS:['C0^+','LAM2g','T2g','H','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[4,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,VP,VV"]})
    
    PRINT_FORMULA({LABEL:'FORM_PT_LAG_A',MODE:'SHORT'})

elif ampl_type == 'LCC_2c':
    # Take away 2 connected terms which are hypothesied to cause trouble in LCC compared to CEPA
    # This reproduces connected CEPA method
    LAG_A.set_rule()
    PRINT_FORMULA({LABEL:'FORM_PT_LAG_A',MODE:'SHORT'})

    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',FAC:+1,
                       OPERATORS:['C0^+','LAM2g','T2g','H','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[4,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,VP,VV"]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_PT_LAG_A',NEW:False,OP_RES:'PT_LAG',FAC:+1,
                       OPERATORS:['C0^+','LAM2g','T2g','H','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[4,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,VP,VV"]})
    
    PRINT_FORMULA({LABEL:'FORM_PT_LAG_A',MODE:'SHORT'})

else:
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

if ampl_type in ['CEPA-like','CEPAn-like','CISD-like','CCSD-like','CCSD-like_c','no_Ec','CCSD-like-h0','CCSD-like-2h0','CCSD-like-3h0','CCSD-like-4h0','CCSD-like-5h0','CCSD-like_2','CCSD-like+1','CCSD-like+2','CCSD-like+3','CEPA+quad']:
    # Construct energy operator for use in lagradian
    DEF_ME_LIST({LIST:'ME_CEPA',
                OPERATOR:'ECEPA',
                IRREP:1,
                '2MS':0,
                AB_SYM:+1})
    OPTIMIZE({
            LABEL_OPT:'FOPT_PT_LAG',
            LABELS_IN:['FORM_ECEPA','FORM_PT_Amp','FORM_PT_LAG']})
else:
    OPTIMIZE({
            LABEL_OPT:'FOPT_PT_LAG',
            LABELS_IN:['FORM_PT_Amp','FORM_PT_LAG']})

#OPTIMIZE({
#        LABEL_OPT:'FOPT_PT_LAG',
#        LABELS_IN:['FORM_PT_LAG']})

OPTIMIZE({
        LABEL_OPT:'FOPT_PT_EQ',
        LABELS_IN:['FORM_PT_RHS','FORM_PT_MVP','FORM_PT_SxT']})


if ampl_type == 'LCC':
    # Print out terms used for testing in ITF
    PRINT_FORMULA({LABEL:'FORM_PT_LAG',MODE:'SHORT'})


# Comments from master
# We need to check the third-order energy expression again
# (e.g. replace H by H_N (i.e. without scalar contribution) for formula generation
#  and replace back H_N->H aferwards (using "REPLACE")

# =======================
# Third order terms
# =======================

# RS-based (CI-based) third order energy: <T^+ W T> - <W><T^+ T>
term_CI = "<C0^+ * (T2g^+) * " + _h1_ + " * T2g * C0>"
if not connected:
    # <W><T^+ T> term:
    # This term vanishes if E^(1) = <W> = 0
    # This happens for Fink's and Dyall's H0, for example
    term_CI_E1a = ("-<C0^+ * H * C0 * C0^+' * (T2g^+') * T2g' * C0'>",
                   ["C0^+","T2g^+'",
                    "C0^+","T2g'",
                    "C0^+","C0'",
                    "H","T2g^+'",
                    "H","T2g'",
                    "H","C0'"])
    term_CI_E1b = "<" +_h0exp_+ " * C0^+ * (T2g^+) * T2g * C0>"
    term_CI = [term_CI, term_CI_E1a, term_CI_E1b]


# CC-based third order energy: <T^+ [W,T]> + 1/2 <[[W,T],T]> + 1/2 <T^+ [[H0,T],T]>
#                                CC0    +       CCa       +       CCb
if connected:
    term_CC0 = "<C0^+ * (T2g^+) * ([" +_h1_+ ",T2g]) * C0>"
else:
    term_CC0 = ["<C0^+ * (T2g^+) * [H,T2g] * C0>",
                "<C0^+ * (T2g^+) * (" +_h0exp_+"-"+_h0_+ ") * T2g * C0>"]

# Originally <[[W,T],T]>, but T on the bra is zero and just <WTT> remains
# Also, <WTT> = <HTT>
term_CCa = "1/2 * <C0^+ * H * T2g * T2g * C0>"

if connected:
    term_CCb = "1/2 * <C0^+ * (T2g^+) * ([[" +_h0_+ ",T2g],T2g]) * C0>"
    term_CCb_Q = None
else:
    term_CCb = ["1/2 * <C0^+ * (T2g^+) * " +_h0_+ " * T2g * T2g * C0>",
                "1/2 * <C0^+ * (T2g^+) * T2g * T2g * " +_h0exp_+ " * C0>"]


    # Term -<T^+ T Q H0 T>
    # where Q is the projector into the external space
    #
    # Version 1: Q = 1-|C0><C0^+| is the projector into the full external space
    # (including the part of the reference space orthogonal to the CASCI wave function
    term_TTQH0T_all = {OPERATORS: ['C0^+', 'T2g^+', 'T2g', _h0_ , 'T2g', 'C0'],
                       IDX_SV: [1,2,3,4,5,6],
                       FAC:-1.0,
                       FIX_VTX:True}

    term_TTQH0T_ref = {OPERATORS: ['C0^+', 'T2g^+', 'T2g', 'C0^+', 'C0^+', _h0_ , 'T2g', 'C0'],
                       IDX_SV:    [1     , 2      , 3    , 4     , 5     , 6    , 7    , 8],
                       FAC:1.0,
                       FIX_VTX:True,
                       AVOID: [1,6,  1,7,  1,8,
                               2,6,  2,7,  2,8,
                               3,6,  3,7,  3,8]}

    term_TTQH0T = [term_TTQH0T_all, term_TTQH0T_ref]


    # Version 2: Q = 1-P is the projection into the orthogonal complement
    # of the complete reference space; P = sum_\mu |\mu> <\mu| where the sum
    # go over all the CAS determinants

    # We implement -<T^+ T H0 T> +<T^+ T P H0 T>
    # The first (negative) term is straightforward;
    # For the second, we consider the possible terms in H0 T>
    # such that no inactive lines are left to connect to the
    # (<T^+ T) side.
    # This is restricted to one body H0!!!

    ops_for_TTQH0T =    ['C0^+', 'T2g^+', 'T2g', _h0_ , 'T2g', 'C0']
    idx_sv_for_TTQH0T = [1     , 2      , 3    , 4    , 5    , 6]

    term_TTH0T_full = {OPERATORS:ops_for_TTQH0T,
                       IDX_SV:idx_sv_for_TTQH0T,
                       FAC:-1.0,
                       FIX_VTX:True}


    term_TTQH0T_hp = {OPERATORS:ops_for_TTQH0T,
                      IDX_SV:idx_sv_for_TTQH0T,
                      FAC:1.0,
                      FIX_VTX:True,
                      LABEL_DESCR:['4,,H,P',
                                   '5,,P,H',
                                   '5,,VP,VH',
                                   '4,5,H,P']}

    term_TTQH0T_vp = {OPERATORS:ops_for_TTQH0T,
                      IDX_SV:idx_sv_for_TTQH0T,
                      FAC:1.0,
                      FIX_VTX:True,
                      LABEL_DESCR:['4,,V,P',
                                   '5,,P,V',
                                   '5,,PV,VV',
                                   '4,5,,P']}

    term_TTQH0T_hv = {OPERATORS:ops_for_TTQH0T,
                      IDX_SV:idx_sv_for_TTQH0T,
                      FAC:1.0,
                      FIX_VTX:True,
                      LABEL_DESCR:['4,,H,V',
                                   '5,,V,H',
                                   '5,,VV,VH',
                                   '4,5,H,']}

    term_TTQH0T = [term_TTH0T_full, term_TTQH0T_hp, term_TTQH0T_vp, term_TTQH0T_hv]


# The difference between CC-based and RS-based
# third order terms, written as:
# 1/2 <(W + [T^+, H0])TT> - <T^+ T (W-<W> + [H0,T])>
# term_highO             + term_pureV
if connected:
    term_highO = "1/2 * <C0^+ * ( H + [(T2g^+)," +_h0_+ "]) * T2g * T2g * C0>"
    term_pureV = "-<C0^+ * (T2g^+) * T2g * ( " +_h1_+ " + [" +_h0_+ ",T2g]) * C0>"
else:
    term_highO = "1/2 * <C0^+ * ( H +  (T2g^+)*(" +_h0_+"-"+_h0exp_+ ") ) * T2g * T2g * C0>"
    # term_pureV is as below, plus term_TTQH0T
    term_pureV = ["-<C0^+ * (T2g^+) * T2g * H * C0>",
                  ("<C0^+ * (T2g^+) * T2g * C0 * C0^+' * H * C0' >",
                   ["C0^+","H",
                    "C0^+","C0'",
                    "T2g^+","H",
                    "T2g^+","C0'",
                    "C0","H",
                    "C0","C0'"]),
                  "<C0^+ * (T2g^+) * T2g * T2g * " +_h0exp_+ " * C0>"
                  ]

third_ord_terms = ['_CI', '_CC', '_CC0', '_CCa', '_CCb', '_CC_higherO', '_CC_pureV']
#third_ord_terms = ['_CI', '_CC']
#third_ord_terms = ['_CCb']


if third_ord_energy:
    new_target('MRCCPT_E_3rd_O', True)
    heading('Third order correction for the energy')
    depend('SOLVE_MRCCPT2')

    for i in third_ord_terms:
        terms_extra = None
        if (i == '_CI'):
            str_i = "RS-based term"
            term = term_CI

        elif (i == '_CC'):
            str_i = "CC-based term"
            term = []
            for t in (term_CC0, term_CCa, term_CCb):
                if isinstance(t,str):
                    term.append(t)
                else:
                    term.extend(t)
            if not(connected):
                terms_extra = term_TTQH0T

        elif (i == '_CC0'):
            str_i = "term 0 of CC-based term: <T^+ [H1,T]>"
            term = term_CC0

        elif (i == '_CCa'):
            str_i = "term a of CC-based term: 1/2 <[[H1,T],T]>"
            term = term_CCa

        elif (i == '_CCb'):
            str_i = "term b of CC-based term: 1/2 <T^+ [[H0,T],T]>"
            term = term_CCb
            if not(connected):
                terms_extra = term_TTQH0T

        elif (i == '_CC_higherO'):
            str_i = "CC term with higher order excitations"
            if isinstance(term_CI, str):
                term = [term_CI]
            else:
                term = list(term_CI)
            if isinstance(term_highO, str):
                term.append(term_highO)
            else:
                term.extend(term_highO)

        elif (i == '_CC_pureV'):
            str_i = "CC term with pure virtual excitations"
            if isinstance(term_CI, str):
                term = [term_CI]
            else:
                term = list(term_CI)
            if isinstance(term_pureV, str):
                term.append(term_pureV)
            else:
                term.extend(term_pureV)
            if not(connected):
                terms_extra = term_TTQH0T

        else:
            raise Exception(i_am + ": unrecognised kind of third order correction: " + i)

        DEF_SCALAR({LABEL:'MRCCPT_O3'+i})
        DEF_ME_LIST({LIST:'ME_MRCCPT_O3'+i,
                     OPERATOR:'MRCCPT_O3'+i,
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        # dbg
        # print "Third order terms for i = " + i + ":"
        # if (isinstance(term, str)):
        #     print term
        # else:
        #     for t in term:
        #         print str(t)
        # if terms_extra is not None:
        #     print "Extra terms:"
        #     for t in terms_extra:
        #         print t
        # end dbg

        if (isinstance(term, str)):
            E_O3 = stf.Formula("FORM_MRCCPT_O3"+i+":MRCCPT_O3"+i+"="+term)
        else:
            E_O3 = stf.Formula("FORM_MRCCPT_O3"+i+":MRCCPT_O3"+i+"="+term[0])
            for t in term[1:]:
                if isinstance(t,str):
                    E_O3.append(t)
                elif isinstance(t, tuple):
                    E_O3.append(t[0], avoid=t[1])
                else:
                    raise Exception(i_am+": I don't know what to do with this: " + str(t))

        E_O3.set_rule()

        # Extra terms, that need the explicit EXPAND_OP_PRODUCT
        # because the stf does not allow LABEL_DESCR
        if terms_extra is not None:
            for t in terms_extra:
                t_actual = dict(t)
                t_actual.update(
                    {LABEL:"FORM_MRCCPT_O3"+i,
                     NEW:False,
                     OP_RES:"MRCCPT_O3"+i})
                EXPAND_OP_PRODUCT(t_actual)
            SUM_TERMS({LABEL_RES:"FORM_MRCCPT_O3"+i,
                       LABEL_IN:"FORM_MRCCPT_O3"+i})

        # if i == '_CCb':
        #     debug_FORM("FORM_MRCCPT_O3"+i,True)
        #     ABORT({})

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

if test_terms:
    new_target('testing_terms', True)
    heading('Testing energy terms')
    depend('SOLVE_MRCCPT2')
    depend(('DEF_FORM_PT_LAG'))


#    DEF_SCALAR({LABEL:'E_J1'})
#    DEF_ME_LIST({LIST:'ME_J1',
#                 OPERATOR:'E_J1',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_1=stf.Formula("FORM_J1:E_J1=<C0^+*T2g^+*H*C0>")
#    E_LAG_1.set_rule()
#    SUM_TERMS({LABEL_RES:'FORM_J1',
#               LABEL_IN:'FORM_J1'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J1',
#              LABELS_IN:['FORM_J1']})
#    PRINT_FORMULA({LABEL:'FORM_J1',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J1'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J2'})
#    DEF_ME_LIST({LIST:'ME_J2',
#                 OPERATOR:'E_J2',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_2=stf.Formula("FORM_J2:E_J2=<C0^+*T2g^+*"+_h1_+"*T2g*C0>")
#    E_LAG_2.set_rule()
#    REPLACE({LABEL_RES:'FORM_J2',
#             LABEL_IN:'FORM_J2',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J2',
#               LABEL_IN:'FORM_J2'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J2',
#              LABELS_IN:['FORM_J2']})
#    PRINT_FORMULA({LABEL:'FORM_J2',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J2'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J3'})
#    DEF_ME_LIST({LIST:'ME_J3',
#                 OPERATOR:'E_J3',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_3=stf.Formula("FORM_J3:E_J3=<C0^+*T2g^+*T2g*"+_h1_+"*C0>")
#    E_LAG_3.set_rule()
#    REPLACE({LABEL_RES:'FORM_J3',
#             LABEL_IN:'FORM_J3',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J3',
#               LABEL_IN:'FORM_J3'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J3',
#              LABELS_IN:['FORM_J3']})
#    PRINT_FORMULA({LABEL:'FORM_J3',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J3'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J4'})
#    DEF_ME_LIST({LIST:'ME_J4',
#                 OPERATOR:'E_J4',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_4=stf.Formula("FORM_J4:E_J4=<C0^+*T2g^+*"+_h0_+"*T2g*C0>")
#    E_LAG_4.set_rule()
#    REPLACE({LABEL_RES:'FORM_J4',
#             LABEL_IN:'FORM_J4',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J4',
#               LABEL_IN:'FORM_J4'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J4',
#              LABELS_IN:['FORM_J4']})
#    PRINT_FORMULA({LABEL:'FORM_J4',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J4'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J5'})
#    DEF_ME_LIST({LIST:'ME_J5',
#                 OPERATOR:'E_J5',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_5=stf.Formula("FORM_J5:E_J5=<C0^+*T2g^+*T2g*"+_h0_+"*C0>")
#    E_LAG_5.set_rule()
#    REPLACE({LABEL_RES:'FORM_J5',
#             LABEL_IN:'FORM_J5',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J5',
#               LABEL_IN:'FORM_J5'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J5',
#              LABELS_IN:['FORM_J5']})
#    PRINT_FORMULA({LABEL:'FORM_J5',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J5'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J6'})
#    DEF_ME_LIST({LIST:'ME_J6',
#                 OPERATOR:'E_J6',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_6=stf.Formula("FORM_J6:E_J6=1/2*<C0^+*T2g^+*"+_h1_+"*T2g*T2g*C0>")
#    E_LAG_6.set_rule()
#    REPLACE({LABEL_RES:'FORM_J6',
#             LABEL_IN:'FORM_J6',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J6',
#               LABEL_IN:'FORM_J6'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J6',
#              LABELS_IN:['FORM_J6']})
#    PRINT_FORMULA({LABEL:'FORM_J6',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J6'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J7'})
#    DEF_ME_LIST({LIST:'ME_J7',
#                 OPERATOR:'E_J7',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_7=stf.Formula("FORM_J7:E_J7=-<C0^+*T2g^+*T2g*"+_h1_+"*T2g*C0>")
#    E_LAG_7.set_rule()
#    REPLACE({LABEL_RES:'FORM_J7',
#             LABEL_IN:'FORM_J7',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J7',
#               LABEL_IN:'FORM_J7'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J7',
#              LABELS_IN:['FORM_J7']})
#    PRINT_FORMULA({LABEL:'FORM_J7',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J7'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J8'})
#    DEF_ME_LIST({LIST:'ME_J8',
#                 OPERATOR:'E_J8',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_8=stf.Formula("FORM_J8:E_J8=1/2*<C0^+*T2g^+*T2g*T2g*"+_h1_+"*C0>")
#    E_LAG_8.set_rule()
#    REPLACE({LABEL_RES:'FORM_J8',
#             LABEL_IN:'FORM_J8',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J8',
#               LABEL_IN:'FORM_J8'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J8',
#              LABELS_IN:['FORM_J8']})
#    PRINT_FORMULA({LABEL:'FORM_J8',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J8'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J9'})
#    DEF_ME_LIST({LIST:'ME_J9',
#                 OPERATOR:'E_J9',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_9=stf.Formula("FORM_J9:E_J9=1/2*<C0^+*T2g^+*"+_h0_+"*T2g*T2g*C0>")
#    E_LAG_9.set_rule()
#    REPLACE({LABEL_RES:'FORM_J9',
#             LABEL_IN:'FORM_J9',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J9',
#               LABEL_IN:'FORM_J9'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J9',
#              LABELS_IN:['FORM_J9']})
#    PRINT_FORMULA({LABEL:'FORM_J9',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J9'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J10'})
#    DEF_ME_LIST({LIST:'ME_J10',
#                 OPERATOR:'E_J10',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_10=stf.Formula("FORM_J10:E_J10=-<C0^+*T2g^+*T2g*"+_h0_+"*T2g*C0>")
#    E_LAG_10.set_rule()
#    REPLACE({LABEL_RES:'FORM_J10',
#             LABEL_IN:'FORM_J10',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J10',
#               LABEL_IN:'FORM_J10'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J10',
#              LABELS_IN:['FORM_J10']})
#    PRINT_FORMULA({LABEL:'FORM_J10',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J10'})
#    
#    
#    DEF_SCALAR({LABEL:'E_J11'})
#    DEF_ME_LIST({LIST:'ME_J11',
#                 OPERATOR:'E_J11',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_11=stf.Formula("FORM_J11:E_J11=1/2*<C0^+*T2g^+*T2g*T2g*"+_h0_+"*C0>")
#    E_LAG_11.set_rule()
#    REPLACE({LABEL_RES:'FORM_J11',
#             LABEL_IN:'FORM_J11',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J11',
#               LABEL_IN:'FORM_J11'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J11',
#              LABELS_IN:['FORM_J11']})
#    PRINT_FORMULA({LABEL:'FORM_J11',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J11'})
   


    
#    DEF_SCALAR({LABEL:'E_J1'})
#    DEF_ME_LIST({LIST:'ME_J1',
#                 OPERATOR:'E_J1',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_1=stf.Formula("FORM_J1:E_J1=-<C0^+*T2g^+*(T2g*("+_h1_+"))*C0>")
#    E_LAG_1.set_rule()
#    REPLACE({LABEL_RES:'FORM_J1',
#             LABEL_IN:'FORM_J1',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J1',
#               LABEL_IN:'FORM_J1'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J1',
#              LABELS_IN:['FORM_J1']})
#    PRINT_FORMULA({LABEL:'FORM_J1',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J1'})
#   
#
#    # CEPA - IDEA1
#    DEF_SCALAR({LABEL:'E_J2'})
#    DEF_ME_LIST({LIST:'ME_J2',
#                 OPERATOR:'E_J2',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_2=stf.Formula("FORM_J2:E_J2=<C0^+*T2g^+*H*T2g*C0>")
#    E_LAG_2.append("-<C0^+*T2g^+*T2g*"+_h0_+"*C0>")
#    E_LAG_2.append("-<C0^+*T2g^+*H*T2g*C0>", connect=["H","T2g"])
#    E_LAG_2.append("<C0^+*T2g^+*"+_h0_+"*T2g*C0>", connect=[_h0_,"T2g"])
#    E_LAG_2.append("-<C0^+*T2g^+*"+_h0_+"*T2g*C0>")
#    E_LAG_2.append("<C0^+*T2g^+*T2g*"+_h0_+"*C0>")
#    E_LAG_2.set_rule()
#    REPLACE({LABEL_RES:'FORM_J2',
#             LABEL_IN:'FORM_J2',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J2',
#               LABEL_IN:'FORM_J2'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J2',
#              LABELS_IN:['FORM_J2']})
#    PRINT_FORMULA({LABEL:'FORM_J2',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J2'})
#    
#
#    # HT
#    DEF_SCALAR({LABEL:'E_J3'})
#    DEF_ME_LIST({LIST:'ME_J3',
#                 OPERATOR:'E_J3',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_3=stf.Formula("FORM_J3:E_J3=<C0^+*T2g^+*H*T2g*C0>")
#    E_LAG_3.set_rule()
#    SUM_TERMS({LABEL_RES:'FORM_J3',
#               LABEL_IN:'FORM_J3'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J3',
#              LABELS_IN:['FORM_J3']})
#    PRINT_FORMULA({LABEL:'FORM_J3',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J3'})
#    
#    
#    # H_0T
#    DEF_SCALAR({LABEL:'E_J4'})
#    DEF_ME_LIST({LIST:'ME_J4',
#                 OPERATOR:'E_J4',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_4=stf.Formula("FORM_J4:E_J4=<C0^+*T2g^+*"+_h0_+"*T2g*C0>")
#    E_LAG_4.set_rule()
#    SUM_TERMS({LABEL_RES:'FORM_J4',
#               LABEL_IN:'FORM_J4'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J4',
#              LABELS_IN:['FORM_J4']})
#    PRINT_FORMULA({LABEL:'FORM_J4',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J4'})
#    
#    
#    # -TH_0
#    DEF_SCALAR({LABEL:'E_J5'})
#    DEF_ME_LIST({LIST:'ME_J5',
#                 OPERATOR:'E_J5',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_5=stf.Formula("FORM_J5:E_J5=-<C0^+*T2g^+*T2g*"+_h0_+"*C0>")
#    E_LAG_5.set_rule()
#    SUM_TERMS({LABEL_RES:'FORM_J5',
#               LABEL_IN:'FORM_J5'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J5',
#              LABELS_IN:['FORM_J5']})
#    PRINT_FORMULA({LABEL:'FORM_J5',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J5'})
#
#
#    # (WT)_c
#    DEF_SCALAR({LABEL:'E_J6'})
#    DEF_ME_LIST({LIST:'ME_J6',
#                 OPERATOR:'E_J6',
#                 IRREP:1,
#                 '2MS':0,
#                 AB_SYM:+1})
#
#    E_LAG_6=stf.Formula("FORM_J6:E_J6=<C0^+*H*C0>")
#    E_LAG_6.append("-<C0^+*H*C0>")
#    E_LAG_6.append("<C0^+*T2g^+*H*T2g*C0>", connect=["H","T2g"])
#    E_LAG_6.append("-<C0^+*T2g^+*"+_h0_+"*T2g*C0>", connect=[_h0_,"T2g"])
#    E_LAG_6.set_rule()
#    REPLACE({LABEL_RES:'FORM_J6',
#             LABEL_IN:'FORM_J6',
#             OP_LIST:[_h0_,'H']})
#    SUM_TERMS({LABEL_RES:'FORM_J6',
#               LABEL_IN:'FORM_J6'})
#    OPTIMIZE({LABEL_OPT:'FOPT_J6',
#              LABELS_IN:['FORM_J6']})
#    PRINT_FORMULA({LABEL:'FORM_J6',MODE:'SHORT'})
#    EVALUATE({FORM:'FOPT_J6'})

    # Disconnected terms
    DEF_SCALAR({LABEL:'E_J6'})
    DEF_ME_LIST({LIST:'ME_J6',
                 OPERATOR:'E_J6',
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:+1})
    EXPAND_OP_PRODUCT({LABEL:'FORM_J6',NEW:True,OP_RES:'E_J6',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,VP,VV"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_J6',NEW:False,OP_RES:'E_J6',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,VP,VV"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_J6',NEW:False,OP_RES:'E_J6',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,P,V"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_J6',NEW:False,OP_RES:'E_J6',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,P,V"],
                       AVOID:[3,4]})
    OPTIMIZE({LABEL_OPT:'FOPT_J6',
              LABELS_IN:['FORM_J6']})
    PRINT_FORMULA({LABEL:'FORM_J6',MODE:'SHORT'})
    EVALUATE({FORM:'FOPT_J6'})
   
   
    # connect counterparts to (WT)d
    DEF_SCALAR({LABEL:'E_J7'})
    DEF_ME_LIST({LIST:'ME_J7',
                 OPERATOR:'E_J7',
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:+1})

    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:True,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,VP,VV"],
                       AVOID:[3,4]})

    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,VP,VV"],
                       AVOID:[3,4]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,PV,VV","4,,P,V"],
                       AVOID:[3,4]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,3, 4,2],
                       LABEL_DESCR:["2,,VV,PP","3,,P,V","4,,P,V"],
                       AVOID:[3,4]})
   
    # Connected terms
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[2,4,  3,4],
                       LABEL_DESCR:["2,,V,P","3,,VV,VP","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,V,P","3,,V,P","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,VV,PV","3,,VV,VP","4,,PP,VV"],
                       AVOID:[2,3]})
    
    EXPAND_OP_PRODUCT({LABEL:'FORM_J7',NEW:False,OP_RES:'E_J7',
                       OPERATORS:['C0^+','T2g^+','H','T2g','C0'],
                       IDX_SV   :[1   ,2   ,3, 4, 5],
                       CONNECT:[3,4, 4,2],
                       LABEL_DESCR:["2,,VV,PV","3,,V,P","4,,PP,VV"],
                       AVOID:[2,3]})

    OPTIMIZE({LABEL_OPT:'FOPT_J7',
              LABELS_IN:['FORM_J7']})
    PRINT_FORMULA({LABEL:'FORM_J7',MODE:'SHORT'})
    EVALUATE({FORM:'FOPT_J7'})


    # Linear-CEPA
    DEF_SCALAR({LABEL:'E_J8'})
    DEF_ME_LIST({LIST:'ME_J8',
                 OPERATOR:'E_J8',
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:+1})

    E_LAG_6=stf.Formula("FORM_J8:E_J8=<C0^+*T2g^+*[H,T2g]*C0>")
    E_LAG_6.append("-<C0^+*T2g^+*H*T2g*C0>")
    E_LAG_6.append("<C0^+*T2g^+*T2g*"+_h0_+"*C0>")
    E_LAG_6.set_rule()
    REPLACE({LABEL_RES:'FORM_J8',
             LABEL_IN:'FORM_J8',
             OP_LIST:[_h0_,'H']})
    SUM_TERMS({LABEL_RES:'FORM_J8',
               LABEL_IN:'FORM_J8'})
    OPTIMIZE({LABEL_OPT:'FOPT_J8',
              LABELS_IN:['FORM_J8']})
    PRINT_FORMULA({LABEL:'FORM_J8',MODE:'SHORT'})
    EVALUATE({FORM:'FOPT_J8'})

    if ampl_type in ['CEPA-like','CEPAn-like','CISD-like','CCSD-like','CCSD-like_c','no_Ec']:
        # <T2g^+ E T>
        DEF_SCALAR({LABEL:'E_quad1'})
        DEF_ME_LIST({LIST:'ME_quad1',
                     OPERATOR:'E_quad1',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})
    
        E_LAG_quad1=stf.Formula("FORM_quad1:E_quad1=<C0^+*T2g^+*(H-ECEPA)*T2g*C0>")
        E_LAG_quad1.set_rule()
        SUM_TERMS({LABEL_RES:'FORM_quad1',
                   LABEL_IN:'FORM_quad1'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad1',
                  LABELS_IN:['FORM_quad1']})
        PRINT_FORMULA({LABEL:'FORM_quad1',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad1'})

        # 1/2<T2g^+ HTT>
        DEF_SCALAR({LABEL:'E_quad2'})
        DEF_ME_LIST({LIST:'ME_quad2',
                     OPERATOR:'E_quad2',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})
    
        E_LAG_quad2=stf.Formula("FORM_quad2:E_quad2=<C0^+*T2g^+*H*T2g*T2g*C0>")
        E_LAG_quad2.set_rule()
        SUM_TERMS({LABEL_RES:'FORM_quad2',
                   LABEL_IN:'FORM_quad2'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad2',
                  LABELS_IN:['FORM_quad2']})
        PRINT_FORMULA({LABEL:'FORM_quad2',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad2'})
        
        
        # MRCCSD- 3O-SCF (Don't care about energy expression)
        DEF_SCALAR({LABEL:'E_quad3'})
        DEF_ME_LIST({LIST:'ME_quad3',
                     OPERATOR:'E_quad3',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})
    
        E_LAG_quad3=stf.Formula("FORM_quad3:E_quad3=<C0^+*T2g^+*[H,T2g]*C0>")
        E_LAG_quad3.append("1/2*<C0^+*T2g^+*[[H,T2g],T2g]*C0>")
        E_LAG_quad3.append("-<C0^+*T2g^+*[H,T2g]*C0>")
        E_LAG_quad3.append("-1/2*<C0^+*T2g^+*[["+_h0_+",T2g],T2g]*C0>")
        E_LAG_quad3.set_rule()
        REPLACE({LABEL_RES:'FORM_quad3',
             LABEL_IN:'FORM_quad3',
             OP_LIST:[_h0_,'H']})
        SUM_TERMS({LABEL_RES:'FORM_quad3',
                   LABEL_IN:'FORM_quad3'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad3',
                  LABELS_IN:['FORM_quad3']})
        PRINT_FORMULA({LABEL:'FORM_quad3',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad3'})
        

        # MRCCSD- CCSD-like (Don't care about energy expression)
        DEF_SCALAR({LABEL:'E_quad4'})
        DEF_ME_LIST({LIST:'ME_quad4',
                     OPERATOR:'E_quad4',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        E_LAG_quad4=stf.Formula("FORM_quad4:E_quad4=1/2*<C0^+*T2g^+*H*T2g*T2g*C0>")
        E_LAG_quad4.append("-1/2*<C0^+*T2g^+*(H-ECEPA)*T2g*T2g*C0>")
        E_LAG_quad4.set_rule()
        REPLACE({LABEL_RES:'FORM_quad4',
             LABEL_IN:'FORM_quad4',
             OP_LIST:[_h0_,'H']})
        SUM_TERMS({LABEL_RES:'FORM_quad4',
                   LABEL_IN:'FORM_quad4'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad4',
                  LABELS_IN:['FORM_quad4']})
        PRINT_FORMULA({LABEL:'FORM_quad4',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad4'})


        DEF_SCALAR({LABEL:'E_quad5'})
        DEF_ME_LIST({LIST:'ME_quad5',
                     OPERATOR:'E_quad5',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        E_LAG_quad5=stf.Formula("FORM_quad5:E_quad5=1/2*<C0^+*T2g^+*H*T2g*T2g*C0>")
        E_LAG_quad5.append("-1/2*<C0^+*T2g^+*("+_h0_+"-ECEPA)*T2g*T2g*C0>")
        E_LAG_quad5.set_rule()
        REPLACE({LABEL_RES:'FORM_quad5',
             LABEL_IN:'FORM_quad5',
             OP_LIST:[_h0_,'H']})
        SUM_TERMS({LABEL_RES:'FORM_quad5',
                   LABEL_IN:'FORM_quad5'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad5',
                  LABELS_IN:['FORM_quad5']})
        PRINT_FORMULA({LABEL:'FORM_quad5',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad5'})
        
    
    if ampl_type in ['CC_2comm']:
        DEF_SCALAR({LABEL:'E_quad6'})
        DEF_ME_LIST({LIST:'ME_quad6',
                     OPERATOR:'E_quad6',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:+1})

        E_LAG_quad6=stf.Formula("FORM_quad6:E_quad6=<C0^+*T2g^+*[H,T2g]*C0>")
        E_LAG_quad6.append("1/2*<C0^+*T2g^+*[[H,T2g],T2g]*C0>")
        E_LAG_quad6.set_rule()
        REPLACE({LABEL_RES:'FORM_quad6',
             LABEL_IN:'FORM_quad6',
             OP_LIST:[_h0_,'H']})
        SUM_TERMS({LABEL_RES:'FORM_quad6',
                   LABEL_IN:'FORM_quad6'})
        OPTIMIZE({LABEL_OPT:'FOPT_quad6',
                  LABELS_IN:['FORM_quad6']})
        PRINT_FORMULA({LABEL:'FORM_quad6',MODE:'SHORT'})
        EVALUATE({FORM:'FOPT_quad6'})



# Calculate stability matrix (Jacobian)
# stabm=number of roots
if stabm > 0:
    new_target('stability_matrix', True)
    heading('Calculating eigenvalues of stability matrix (Jacobian)')
    depend('SOLVE_MRCCPT2')
    depend(('DEF_FORM_PT_LAG'))
    depend('BUILD_PRECON')

    CLONE_OPERATOR({LABEL:'R2g',TEMPLATE:'T2g'})
    DEF_ME_LIST({LIST:'ME_R2g',
            OPERATOR:'R2g',
            IRREP:1,
            '2MS':0,
            AB_SYM:+1,MAX_REC:stabm})
    CLONE_OPERATOR({LABEL:'R2g_prime',TEMPLATE:'T2g'})
    DEF_ME_LIST({LIST:'ME_R2g_prime',
            OPERATOR:'R2g_prime',
            IRREP:1,
            '2MS':0,
            AB_SYM:0,MAX_REC:stabm})

    # Construct overlap matrix
    CLONE_OPERATOR({LABEL:'S2g',TEMPLATE:'O2g'})
    DEF_ME_LIST({LIST:'ME_S2g',
            OPERATOR:'S2g',
            IRREP:1,
            '2MS':0,
            AB_SYM:+1,MAX_REC:stabm})
    OVRLAP=stf.Formula("FORM_OVRLAP:S2g=<C0^+*LAM2g*R2g*C0>")
    OVRLAP.set_rule()


    DERIVATIVE({LABEL_IN:'FORM_OVRLAP',
            LABEL_RES:'FORM_S_R',
            OP_RES:'S2g',
            OP_DERIV:'LAM2g'})


    CLONE_OPERATOR({LABEL:'JAC_R',TEMPLATE:'O2g'})
    DEF_ME_LIST({LIST:'ME_JAC_R',
            OPERATOR:'JAC_R',
            IRREP:1,
            '2MS':0,
            AB_SYM:+1,MAX_REC:stabm})

    
    DERIVATIVE({LABEL_IN:'FORM_PT_Amp',
            LABEL_RES:'FORM_PT_Stabm',
            OP_RES:'JAC_R',
            OP_DERIV:'T2g',
            OP_MULT:'R2g'})

    # for transformation:
    DEF_SCALAR({LABEL:'DUMM_Y'})
    INVARIANT({LABEL_RES:'FORM_R2g',
           LABEL_IN:'FORM_T2_orth',
           OP_RES:'R2g',
           OPERATORS:'DUMM_Y'}) # currently the only way to change the OP_RES

    REPLACE({LABEL_RES:'FORM_R2g',
         LABEL_IN:'FORM_R2g',
         OP_LIST:['T2_orth','R2g_prime']})

    OPTIMIZE({LABEL_OPT:'FOPT_JAC_R',
              LABELS_IN:['FORM_PT_Stabm','FORM_S_R','FORM_R2g']})

    ASSIGN_ME2OP({LIST:'ME_X_TRM_DAG',
                OPERATOR:'X_TRM'})


    SOLVE_EVP({LIST_OPT:'ME_R2g',
               LIST_PRC:'ME_PRECON2g',
               #SOLVER:'NEW',
               SOLVER:'OLD',
               OP_MVP:'JAC_R',
               OP_SVP:'S2g',
               FORM:'FOPT_JAC_R',
               MODE:'TRF',
               #FORM_SPC:['FOPT_T2_orth'],
               LIST_SPC:['ME_R2g_prime','ME_X_TRM','ME_X_TRM_DAG'],
               N_ROOTS:stabm})
