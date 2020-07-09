"""An implementation of the icMRCCSD Lagrangian, done term by term and separated T1 and T2.



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


# Select a method to use
approx = keywords.get('method.MR_P.method')
method = approx if approx is not None else "Yuri"

known_methods=["Yuri", "MRCCSD22", "CEPT2","CCEPA","BETACEPT2","CEPA0","CEPA0TeTi","PT2"]
if method not in known_methods :
    raise Exception(i_am+": unknown method:"+str(method))
print("Using the "+method+" method.")


# Select H0
if method in ['CEPT2','BETACEPT2','PT2']:
    known_hamiltonians=["DYALL","REPT","F_EFF"]
    hamiltonian = keywords.get('method.MR_P.hamiltonian')
    hamiltonian=str(hamiltonian).strip() if hamiltonian is not None else "DYALL"

    if hamiltonian not in known_hamiltonians :
        raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

    depend('H0')

    if hamiltonian=="DYALL":
        depend('EVAL_HAM_D')
        _h0_='HAM_D'
    elif hamiltonian=="REPT":
        depend('EVAL_REPT_HAM')
        _h0_='REPT_HAM'
    elif hamiltonian=="F_EFF":
        depend('EVAL_F_EFF')
        _h0_='FOCK_EFF'
        _h0exp_='FOCK_EFF_EXP'

#----------
#Define the Te and Ti operators
#with_singles_2
#te_shape='PP,HV|PP,VV|PP,HH'
#--------------
#with_singles_3
#te_shape='PP,HV|PP,VV|PP,HH|P,V'
#--------------
#with_singles_4
#te_shape='PP,HV|PP,VV|PP,HH|P,V|P,H'
#--------------
#with_singles_5
#te_shape='PP,HV|PP,VV|PP,HH|P,H'
#--------------
#with_singles_6
te_shape='PP,HV|PP,VV|PP,HH|V,H|P,V|P,H'

#with_singles_2   
ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'
#--------------
#with_singles_3
#ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,H'
#--------------
#with_singles_4
#ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H'
#--------------
#with_singles_5
#ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V'
#--------------
#with_singles_6
#ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH'
#--------------

DEF_OP_FROM_OCC({LABEL:'Te',
                 DESCR:te_shape})

#def LAMe
CLONE_OPERATOR({LABEL:'LAMe',
                TEMPLATE:'Te',
                ADJOINT:True})

DEF_OP_FROM_OCC({LABEL:'Ti',
                 DESCR:ti_shape})
#def LAMi
CLONE_OPERATOR({LABEL:'LAMi',
                TEMPLATE:'Ti',
                ADJOINT:True})


#-----------------TEST___TEST_for seperating singles
ts_shape='V,H|P,V|P,H'

DEF_OP_FROM_OCC({LABEL:'Ts',
                 DESCR:ts_shape})
CLONE_OPERATOR({LABEL:'LAMs',
                TEMPLATE:'Ts',
                ADJOINT:True})
#----------------END_TEST


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

if method in ['CEPT2','CCEPA','BETACEPT2','CEPA0TeTi']: #TeTistuff
   def _Le_refexp(x):
       return _refexp("LAMe(" + x + ")")
   
   def _Li_refexp(x):
       return _refexp("LAMi(" + x + ")")

if method in ['BETACEPT2','CEPA0TeTi']:
   def _Ls_refexp(x):
       return _refexp("LAMs(" + x + ")")


# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")


LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))

if method in ['BETACEPT2','CEPA0TeTi']: #Ts -> T1     Te,Ti -> T2g  --note:testidea1, check idea2

    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Ls_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Li_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Le_refexp("H"))
else:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))

if method=='MRCCSD22':
    # Put MRCCSD(2,2) method here...
    LAG_E.append(_refexp("H*T2g"))
    LAG_E.append(_refexp("(1/2)*H*T2g*T2g"))
    
    LAG_A1.append(_L1_refexp("[H,T2g]"))  #dunno
    LAG_A1.append(_L1_refexp("(1/2)*[[H,T2g],T2g]")) #

    LAG_A2.append(_L2_refexp("[H,T2g]"))
    LAG_A2.append(_L2_refexp("(1/2)*[[H,T2g],T2g]"))
    LAG_A2.append(_L1_refexp("-[H,T2g]"))  #dunno
    LAG_A2.append(_L1_refexp("-(1/2)*[[H,T2g],T2g]")) #
        
elif method=='PT2':

    LAG_E.append(_refexp("(H*T1)+(H*T2g)"))

    LAG_A1.append(_L1_refexp("("+_h0_+"-"+_h0exp_+")*T1"))
    LAG_A1.append(_L1_refexp("("+_h0_+"-"+_h0exp_+")*T2g"))
    
    LAG_A2.append(_L2_refexp("("+_h0_+"-"+_h0exp_+")*T1"))
    LAG_A2.append(_L2_refexp("("+_h0_+"-"+_h0exp_+")*T2g"))
    
elif method=='Yuri':

    #LAG_E.append(_refexp("[H,T1]"))
    #LAG_E.append(_refexp("[H,T2g]"))
    #
    #LAG_A1.append(_L1_refexp("[H,T1]"))
    #LAG_A1.append(_L1_refexp("[H,T2g]"))
    #LAG_A1.append(_L1_refexp("(1/2)*[[H,T1],T2g]"))
    #
    #LAG_A2.append(_L2_refexp("[H,T1]"))
    #LAG_A2.append(_L2_refexp("[H,T2g]"))
    #LAG_A2.append(_L2_refexp("(1/2)*[[H,T2g],T2g]"))


    LAG_E.append(_refexp("(H*T1)+(H*T2g)"))
    LAG_E.append(_refexp("-(T1*H)-(T2g*H)"))

    LAG_E.append(_refexp("(1/2)*(H*T1 *T1 )"))
    LAG_E.append(_refexp("(1/2)*(H*T1 *T2g)"))
    LAG_E.append(_refexp("(1/2)*(H*T2g*T1 )"))
    LAG_E.append(_refexp("(1/2)*(H*T2g*T2g)"))

    #LAG_E.append(_refexp("(1/2)*(H*T2g'*T2g'')"), avoid=["T2g'","T2g''"])

    LAG_E.append(_refexp("-(T1 *H*T1 )"))
    LAG_E.append(_refexp("-(T1 *H*T2g)"))
    LAG_E.append(_refexp("-(T2g*H*T1 )"))
    LAG_E.append(_refexp("-(T2g*H*T2g)"))

    LAG_E.append(_refexp("(1/2)*(T1 *T1 *H)"))
    LAG_E.append(_refexp("(1/2)*(T1 *T2g*H)"))
    LAG_E.append(_refexp("(1/2)*(T2g*T1 *H)"))
    LAG_E.append(_refexp("(1/2)*(T2g*T2g*H)"))



    LAG_A1.append(_L1_refexp("(H*T1)+(H*T2g)"))
    LAG_A1.append(_L1_refexp("-(T1*H)-(T2g*H)"))

    LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T1 )"))
    LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T2g)"))
    LAG_A1.append(_L1_refexp("(1/2)*(H*T2g*T1 )"))
    LAG_A1.append(_L1_refexp("(1/2)*(H*T2g*T2g)"))

    LAG_A1.append(_L1_refexp("-(T1 *H*T1 )"))
    LAG_A1.append(_L1_refexp("-(T1 *H*T2g)"))
    LAG_A1.append(_L1_refexp("-(T2g*H*T1 )"))
    LAG_A1.append(_L1_refexp("-(T2g*H*T2g)"))

    LAG_A1.append(_L1_refexp("(1/2)*(T1 *T1 *H)"))
    LAG_A1.append(_L1_refexp("(1/2)*(T1 *T2g*H)"))
    LAG_A1.append(_L1_refexp("(1/2)*(T2g*T1 *H)"))
    LAG_A1.append(_L1_refexp("(1/2)*(T2g*T2g*H)"))



    LAG_A2.append(_L2_refexp("(H*T1)+(H*T2g)"))
    LAG_A2.append(_L2_refexp("-(T1*H)-(T2g*H)"))

    LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T1 )"))
    LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T2g)"))
    LAG_A2.append(_L2_refexp("(1/2)*(H*T2g*T1 )"))
    LAG_A2.append(_L2_refexp("(1/2)*(H*T2g*T2g)"))

    LAG_A2.append(_L2_refexp("-(T1 *H*T1 )"))
    LAG_A2.append(_L2_refexp("-(T1 *H*T2g)"))
    LAG_A2.append(_L2_refexp("-(T2g*H*T1 )"))
    LAG_A2.append(_L2_refexp("-(T2g*H*T2g)"))

    LAG_A2.append(_L2_refexp("(1/2)*(T1 *T1 *H)"))
    LAG_A2.append(_L2_refexp("(1/2)*(T1 *T2g*H)"))
    LAG_A2.append(_L2_refexp("(1/2)*(T2g*T1 *H)"))
    LAG_A2.append(_L2_refexp("(1/2)*(T2g*T2g*H)"))

elif method == 'CEPT2':
    #here goes cept2
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    
    

elif method == 'BETACEPT2':
    #here goes cept2
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

elif method == 'CEPA0TeTi':
    #here goes cept2
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

    LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))

    LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
    LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ti"))
    LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))

    LAG_A2.append(_Li_refexp("(H-ECEPA)*Ts"))
    LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
    LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))
elif method == 'CEPA0':
    #here goes CEPA0
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

    LAG_E.append(_refexp("(H*T1)+(H*T2g)"))

    LAG_A1.append(_L1_refexp("(H-ECEPA)*T1"))
    LAG_A1.append(_L1_refexp("(H-ECEPA)*T2g"))

    LAG_A2.append(_L2_refexp("(H-ECEPA)*T1"))
    LAG_A2.append(_L2_refexp("(H-ECEPA)*T2g"))
elif method == 'CCEPA':
    #here goes the ccepa
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    LAG_E.append("<C0^+*1/2*H*Te*Te*C0>")
    LAG_A.append("<C0^+*(LAMi)*(H-ECEPA)*Ti*C0>")
    LAG_A.append("<C0^+*(LAMi)*(H-ECEPA)*Te*C0>")
    LAG_A.append("<C0^+*(LAMe)*(H-ECEPA)*Ti*C0>")
    LAG_A.append("<C0^+*(LAMe)*[H,Te]*C0>")
    LAG_A.append("1/2*<C0^+*(LAMe)*[[H,Te],Te]*C0>")





LAG_E.set_rule()


if method in ['CEPT2','CCEPA']:
#Replaces the Ti and Te strings, so that the equations can be parsed to the solver
#may need to introduce a new Tis operator that will take the place of this Ti
#!!!NOT FINISHED!!!
    LAG_A1.set_rule()
    LAG_A2.set_rule()
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',
         LABEL_IN:'FORM_MRCC_LAG_E',
         OP_LIST:['Ti','T2g','Te','T2g']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
         LABEL_IN:'FORM_MRCC_LAG_A1',
         OP_LIST:['Ti','T2g','Te','T2g']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
         LABEL_IN:'FORM_MRCC_LAG_A1',
         OP_LIST:['LAMi','LAM2g','LAMe','LAM2g']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
         LABEL_IN:'FORM_MRCC_LAG_A2',
         OP_LIST:['Ti','T2g','Te','T2g']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
         LABEL_IN:'FORM_MRCC_LAG_A2',
         OP_LIST:['LAMi','LAM2g','LAMe','LAM2g']})

    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

elif method in ['BETACEPT2','CEPA0TeTi']:

    LAG_A1.set_rule()
    LAG_A2.set_rule()
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',
         LABEL_IN:'FORM_MRCC_LAG_E',
         OP_LIST:['Ti','T2g','Te','T2g','Ts','T1']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
         LABEL_IN:'FORM_MRCC_LAG_A1',
         OP_LIST:['Ti','T2g','Te','T2g','Ts','T1']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
         LABEL_IN:'FORM_MRCC_LAG_A1',
         OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMs','LAM1']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
         LABEL_IN:'FORM_MRCC_LAG_A2',
         OP_LIST:['Ti','T2g','Te','T2g','Ts','T1']})

    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
         LABEL_IN:'FORM_MRCC_LAG_A2',
         OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMs','LAM1']})

    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

else:
    LAG_A1.set_rule()
    LAG_A2.set_rule()

    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

#Replaces the Ti and Te strings, so that the equations can be parsed to the solver

#Make the Derivative with respect to LAM
DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A1',
        LABEL_RES:'FORM_MRCC_LAG_Amp1',
        OP_RES:'O1',
        OP_DERIV:'LAM1'})

DERIVATIVE({
        LABEL_IN:'FORM_MRCC_LAG_A2',
        LABEL_RES:'FORM_MRCC_LAG_Amp2',
        OP_RES:'O2g',
        OP_DERIV:'LAM2g'})


if method in ['CEPT2','CCEPA','CEPA0','CEPA0TeTi','BETACEPT2']:
       # Construct energy operator for use in lagrangian
   DEF_ME_LIST({LIST:'ME_CEPA',
                OPERATOR:'ECEPA',
                IRREP:1,
                '2MS':0,
                AB_SYM:+1})
   OPTIMIZE({
           LABEL_OPT:'FOPT_MRCC_LAG',
           LABELS_IN:['FORM_ECEPA','FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})
else:

   OPTIMIZE({
           LABEL_OPT:'FOPT_MRCC_LAG',
           LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})


#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
