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

third_ord_energy=False
if keywords.is_keyword_set('method.MRCCPT2.3rd_E'):
    if (keywords.get('method.MRCCPT2.3rd_E') == "T"):
        third_ord_energy=True
    elif(keywords.get('method.MRCCPT2.3rd_E') == "F"):
        third_ord_energy=False
    else :
        raise Exception(i_am+": unrecognised value for option 3rd_E (must be T or F)")
print("3rd_E ", third_ord_energy, type(third_ord_energy))

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

#Energy equation no 
LAG_E=stf.Formula("FORM_PT_LAG:PT_LAG="\
                  "<C0^+*(H+H*T2g)*C0>")

LAG_A=stf.Formula("FORM_PT_LAG_A:PT_LAG="\
                  "<C0^+*(LAM2g*H)*C0>")

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

if (connected):
      LAG_A.append("<C0^+*(LAM2g)*(["+_h0_+",T2g])*C0>")
else:
      LAG_A.append("<C0^+*(LAM2g)*("+_h0_+"-"+_h0exp_+")*T2g*C0>")


# We need to check the third-order energy expression again
# (e.g. replace H by H_N (i.e. without scalar contribution) for formula generation
#  and replace back H_N->H aferwards (using "REPLACE")
if (third_ord_energy):
    if (connected):
        LAG_E.append("<C0^+*(T2g^+)*[H,T2g]*C0>")
    else:
        LAG_E.append("<C0^+*(T2g^+)*((H-"+_h0_+")*T2g)*C0>")


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

#we need the overlap terms for the LEQ solver
OVL=stf.Formula("FORM_OVL:OVL_A=<C0^+*LAM2g*T2g*C0>")
OVL.set_rule()


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


mark("PT-LAGRANGIAN")

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG',
        LABEL_RES:'FORM_PT_LAG'})

SUM_TERMS({
        LABEL_IN:'FORM_PT_LAG_A',
        LABEL_RES:'FORM_PT_LAG_A'})



debug_FORM('FORM_PT_LAG')


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
debug_FORM('FORM_PT_SxT',True)

#devide into RHS and transformation part
LEQ_SPLIT({LABEL_RAW:'FORM_PT_Amp',
           LABEL_TRF:'FORM_PT_MVP',
           LABEL_RHS:'FORM_PT_RHS',
           OP_TRF:'H0xT',
           OP_RHS:'H1rhs',
           OP_X:'T2g'})

debug_FORM('FORM_PT_MVP',True)
debug_FORM('FORM_PT_RHS',True)

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

