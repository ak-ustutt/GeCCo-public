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


# No occupied orbitals? Don't need the To operator
no_occ=False
if keywords.is_keyword_set('method.MR_P.no_occ'):
    string = (keywords.get('method.MR_P.no_occ'))
    if string == 'T':
        no_occ = True
    elif string == 'F':
        no_occ = False
    else:
        print("Didn't recognise no_occ argument; setting to False")
        no_occ = False
    print("Are there no occupied orbitals: ", no_occ)
else:
    print("This calculation has occupied orbitals")


#Workaround for the single-seperation of CEPT2,CCEPA and TCPT2
#   print("Setting singles to: 0")
#if method in ["TCPT20","TCPT21","CEPT20","CEPT21","CEPT22","CEPT23","CEPT24","CEPT25","CEPT26","CEPT27","CCEPA0","CCEPA1","CCEPA2","CCEPA3","CCEPA4","CCEPA5","CCEPA6","CCEPA7"]:
#    singles = int(method[5])
#    print("setting singles to:"+str(singles))
#    method_old = method
#    method = method[0:5]

if keywords.is_keyword_set('method.MR_P.singles'): #remove single keyword or CEPA(0) and PT2 calculations will crash accordingly
    singles = int((keywords.get('method.MR_P.singles')))
    if singles == 0:
        print("All singles operators are in the internal space")
    else:
        print("Singles operators are split between the internal and external spaces")
else:
    if method in ["CEPT2","CCEPA","TCPT2"]:
        singles = 0
    else:
        singles = 8
    print("All singles operators are in the internal space")


known_methods=["Yuri", "MRCCSD22", "CEPT2","CCEPA","TCPT2","CEPA0","PT2"]
if method not in known_methods :
    raise Exception(i_am+": unknown method:"+str(method))
print("Using the "+method+" method.")


# Select H0
ham = ''
if method in ['CEPT2','PT2','TCPT2']:
    known_hamiltonians=["DYALL","REPT","F_EFF"]
    hamiltonian = keywords.get('method.MR_P.hamiltonian')
    hamiltonian=str(hamiltonian).strip() if hamiltonian is not None else "DYALL"

    if hamiltonian not in known_hamiltonians :
        raise Exception(i_am+": unknown hamiltonian type:"+str(hamiltonian))

    depend('H0')

    if hamiltonian=="DYALL":
        depend('EVAL_HAM_D')
        _h0_='HAM_D'
        ham = _h0_
    elif hamiltonian=="REPT":
        depend('EVAL_REPT_HAM')
        _h0_='REPT_HAM'
        ham = _h0_
    elif hamiltonian=="F_EFF":
        depend('EVAL_F_EFF')
        _h0_='FOCK_EFF'
        _h0exp_='FOCK_EFF_EXP'
        ham = _h0_+"-"+_h0exp_


#---Test-for-new-keyword
if method in ["CEPT2","CCEPA","TCPT2"]:

#To -> Internal_Single_EXT
#Tr -> External_Single_EXT
    if singles == 1:

        to_shape='V,H'      #1
        tr_shape='P,V|P,H'  #1

    elif singles == 2:

        to_shape='P,V'      #2
        tr_shape='V,H|P,H'  #2

    elif singles == 3:

        to_shape='P,H'      #3
        tr_shape='V,H|P,V'  #3

    elif singles == 4:

        to_shape='V,H|P,V'  #4
        tr_shape='P,H'      #4

    elif singles == 5:

        to_shape='V,H|P,H'  #5
        tr_shape='P,V'      #5

    elif singles == 6:

        to_shape='P,V|P,H'  #6
        tr_shape='V,H'      #6

    if singles in [1,2,3,4,5,6]:
        if (no_occ and singles==1):
            DEF_OP_FROM_OCC({LABEL:'Tr',
                             DESCR:tr_shape})

            CLONE_OPERATOR({LABEL:'LAMr',
                            TEMPLATE:'Tr',
                            ADJOINT:True})
        else:
            DEF_OP_FROM_OCC({LABEL:'To',
                             DESCR:to_shape})

            CLONE_OPERATOR({LABEL:'LAMo',
                            TEMPLATE:'To',
                            ADJOINT:True})

            DEF_OP_FROM_OCC({LABEL:'Tr',
                             DESCR:tr_shape})

            CLONE_OPERATOR({LABEL:'LAMr',
                            TEMPLATE:'Tr',
                            ADJOINT:True})
#--------------
#CIPT2 seperation ? "all excitations from the active space"
#te_shape='PP,VV|PV,VV|P,V'
#--------------

#Define the Te and Ti operators and the corresponding LAMe and LAMi
te_shape='PP,HV|PP,VV|PP,HH'
ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'

DEF_OP_FROM_OCC({LABEL:'Te',
                 DESCR:te_shape})
CLONE_OPERATOR({LABEL:'LAMe',
                TEMPLATE:'Te',
                ADJOINT:True})

DEF_OP_FROM_OCC({LABEL:'Ti',
                 DESCR:ti_shape})
CLONE_OPERATOR({LABEL:'LAMi',
                TEMPLATE:'Ti',
                ADJOINT:True})


#-----------------TEST___TEST_for seperating singles (delete it and use T1 instead?)
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

if method in ['TCPT2','CEPT2','CCEPA']: #TeTistuff
   def _Le_refexp(x):
       return _refexp("LAMe(" + x + ")")

   def _Li_refexp(x):
       return _refexp("LAMi(" + x + ")")
   print("LAMe definiert")

if singles in [0,7]:
   def _Ls_refexp(x):
       return _refexp("LAMs(" + x + ")")

if singles in [1,2,3,4,5,6]:
   def _Lo_refexp(x):
       return _refexp("LAMo(" + x + ")")

   def _Lr_refexp(x):
       return _refexp("LAMr(" + x + ")")

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")


LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))

if singles in [0,7]:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Ls_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
    LAG_A2.append(_Le_refexp("H"))

elif singles in [1,2,3,4,5,6]:
    if no_occ:
        LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Lr_refexp("H"))
        LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
        LAG_A2.append(_Le_refexp("H"))
    else:
        LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Lo_refexp("H"))
        LAG_A1.append(_Lr_refexp("H"))
        LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
        LAG_A2.append(_Le_refexp("H"))

else:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2_refexp("H"))

if method=='PT2':

    LAG_E.append(_refexp("(H*T1)+(H*T2g)"))

    if hamiltonian == 'DYALL':
    
        LAG_A1.append(_Ls_refexp("(["+_h0_+",T1])"))
        LAG_A1.append(_Ls_refexp("(["+_h0_+",T2g])"))

        LAG_A2.append(_Ls_refexp("(["+_h0_+",T1])"))
        LAG_A2.append(_Ls_refexp("(["+_h0_+",T2g])"))

    else:

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


    if singles == 0:
        LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
        LAG_A2.append(_Le_refexp("(H-ECEPA)*Te")) #this term stays the same

#        LAG_A1.append(_Ls_refexp("(["+ham+",Ts])"))
#        LAG_A1.append(_Ls_refexp("(["+ham+",Ti])"))
#        LAG_A1.append(_Ls_refexp("(["+ham+",Te])"))
#
#        LAG_A2.append(_Li_refexp("(["+ham+",Ts])"))
#        LAG_A2.append(_Li_refexp("(["+ham+",Ti])"))
#        LAG_A2.append(_Li_refexp("(["+ham+",Te])"))
#
#        LAG_A2.append(_Le_refexp("(["+ham+",Ts])"))
#        LAG_A2.append(_Le_refexp("(["+ham+",Ti])"))


        if hamiltonian == "DYALL":

            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
            LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))

            LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])"))
            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

        else:

            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

    elif singles in [1,2,3,4,5,6]:
        LAG_E.append(_refexp("(H*(Te+Ti))+(H*(Tr))"))

        LAG_A1.append(_Lr_refexp("(H-ECEPA)*Te"))
        LAG_A1.append(_Lr_refexp("(H-ECEPA)*Tr"))

        LAG_A2.append(_Le_refexp("(H-ECEPA)*Tr"))
        LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

        if hamiltonian == "DYALL":

            LAG_A1.append(_Lr_refexp("(["+_h0_+",Ti])"))

            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Tr])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

            if not no_occ:

                LAG_E.append(_refexp("(H*(To))"))

                LAG_A1.append(_Lo_refexp("(["+_h0_+",To])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Ti])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Te])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Tr])"))

                LAG_A1.append(_Lr_refexp("(["+_h0_+",To])"))

                LAG_A2.append(_Li_refexp("(["+_h0_+",To])"))

                LAG_A2.append(_Le_refexp("(["+_h0_+",To])"))


        else:

            LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

            if not no_occ:

                LAG_E.append(_refexp("(H*(To))"))

                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*To"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

                LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*To"))

                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*To"))

                LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test


    elif singles == 7:
        LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))

        LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
        LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))

        LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))
        LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

        if hamiltonian == "DYALL":
            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))

            LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

        else:

            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

elif method == 'TCPT2':
    if singles == 0:
        LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
        LAG_E.append(_refexp("1/2*H*Te*Te"))

        LAG_A2.append(_Le_refexp("H*Te"))
        LAG_A2.append(_Le_refexp("-Te*H"))
        LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
        LAG_A2.append(_Le_refexp("-Te*H*Te"))
        LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))

        if hamiltonian == "DYALL": #exchanges the "PT2" terms for the DYALL Ham

            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
            LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))

            LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])"))
            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

        else:

            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

    if singles in [1,2,3,4,5,6]:

        LAG_E.append(_refexp("(H*(Tr))+(H*(Te+Ti))"))
        LAG_E.append(_refexp("1/2*H*Te*Te"))
        LAG_E.append(_refexp("1/2*H*Tr*Tr"))

        LAG_A1.append(_Lr_refexp("H*Te"))
        LAG_A1.append(_Lr_refexp("H*Tr"))
        LAG_A1.append(_Lr_refexp("-Te*H"))
        LAG_A1.append(_Lr_refexp("-Tr*H"))
        LAG_A1.append(_Lr_refexp("1/2*H*Te*Te"))
        LAG_A1.append(_Lr_refexp("1/2*H*Tr*Tr"))
        LAG_A1.append(_Lr_refexp("-Te*H*Te"))
        LAG_A1.append(_Lr_refexp("-Tr*H*Tr"))
        LAG_A1.append(_Lr_refexp("1/2*Te*Te*H"))
        LAG_A1.append(_Lr_refexp("1/2*Tr*Tr*H"))

        LAG_A2.append(_Le_refexp("H*Te"))
        LAG_A2.append(_Le_refexp("H*Tr"))
        LAG_A2.append(_Le_refexp("-Te*H"))
        LAG_A2.append(_Le_refexp("-Tr*H"))
        LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
        LAG_A2.append(_Le_refexp("1/2*H*Tr*Tr"))
        LAG_A2.append(_Le_refexp("-Te*H*Te"))
        LAG_A2.append(_Le_refexp("-Tr*H*Tr"))
        LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
        LAG_A2.append(_Le_refexp("1/2*Tr*Tr*H"))

        if hamiltonian == "DYALL":

            if not no_occ: #excludes To if there are no occupied orbitals, relevant for to_shape='V,H'

                LAG_E.append(_refexp("(H*To)"))

                LAG_A1.append(_Lo_refexp("(["+_h0_+",To])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Ti])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Te])"))
                LAG_A1.append(_Lo_refexp("(["+_h0_+",Tr])"))
                LAG_A1.append(_Lr_refexp("(["+_h0_+",To])"))

                LAG_A2.append(_Li_refexp("(["+_h0_+",To])"))

                LAG_A2.append(_Le_refexp("(["+_h0_+",To])"))


            LAG_A1.append(_Lr_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Tr])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

        else:

            if not no_occ:

                LAG_E.append(_refexp("(H*To)"))

                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*To"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
                LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

                LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test

                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*To"))

                LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*To"))


            LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))


    if singles == 7:

        LAG_E.append(_refexp("(H*(Ts))+(H*(Te+Ti))"))

        if hamiltonian == "DYALL":

            LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))

            LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
            LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))

            LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
        else:

            LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
            LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

            LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))

        LAG_A1.append(_Ls_refexp("H*Te"))
        LAG_A1.append(_Ls_refexp("H*Ts"))
        LAG_A1.append(_Ls_refexp("-Te*H"))
        LAG_A1.append(_Ls_refexp("-Ts*H"))
        LAG_A1.append(_Ls_refexp("1/2*H*Te*Te"))
        LAG_A1.append(_Ls_refexp("1/2*H*Ts*Ts"))
        LAG_A1.append(_Ls_refexp("-Te*H*Te"))
        LAG_A1.append(_Ls_refexp("-Ts*H*Ts"))
        LAG_A1.append(_Ls_refexp("1/2*Te*Te*H"))
        LAG_A1.append(_Ls_refexp("1/2*Ts*Ts*H"))

        LAG_A2.append(_Le_refexp("H*Te"))
        LAG_A2.append(_Le_refexp("H*Ts"))
        LAG_A2.append(_Le_refexp("-Te*H"))
        LAG_A2.append(_Le_refexp("-Ts*H"))
        LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
        LAG_A2.append(_Le_refexp("1/2*H*Ts*Ts"))
        LAG_A2.append(_Le_refexp("-Te*H*Te"))
        LAG_A2.append(_Le_refexp("-Ts*H*Ts"))
        LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))

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

    if singles == 0:
       LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
       LAG_E.append(_refexp("1/2*H*Te*Te"))

       LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
       LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ti"))
       LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))

       LAG_A2.append(_Li_refexp("(H-ECEPA)*Ts"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))

       LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))
       LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
       LAG_A2.append(_Le_refexp("H*Te"))
       LAG_A2.append(_Le_refexp("-Te*H"))
       LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
       LAG_A2.append(_Le_refexp("-Te*H*Te"))
       LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))

    elif singles in [1,2,3,4,5,6]:

        LAG_E.append(_refexp("(H*(Tr))+(H*(Te+Ti))"))
        LAG_E.append(_refexp("1/2*H*Te*Te"))
        LAG_E.append(_refexp("1/2*H*Tr*Tr"))

        LAG_A1.append(_Lr_refexp("(H-ECEPA)*Ti"))
        LAG_A1.append(_Lr_refexp("H*Te"))
        LAG_A1.append(_Lr_refexp("H*Tr"))
        LAG_A1.append(_Lr_refexp("-Te*H"))
        LAG_A1.append(_Lr_refexp("-Tr*H"))
        LAG_A1.append(_Lr_refexp("1/2*H*Te*Te"))
        LAG_A1.append(_Lr_refexp("1/2*H*Tr*Tr"))
        LAG_A1.append(_Lr_refexp("-Te*H*Te"))
        LAG_A1.append(_Lr_refexp("-Tr*H*Tr"))
        LAG_A1.append(_Lr_refexp("1/2*Te*Te*H"))
        LAG_A1.append(_Lr_refexp("1/2*Tr*Tr*H"))

        LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
        LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))
        LAG_A2.append(_Li_refexp("(H-ECEPA)*Tr"))

        LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
        LAG_A2.append(_Le_refexp("H*Te"))
        LAG_A2.append(_Le_refexp("H*Tr"))
        LAG_A2.append(_Le_refexp("-Te*H"))
        LAG_A2.append(_Le_refexp("-Tr*H"))
        LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
        LAG_A2.append(_Le_refexp("1/2*H*Tr*Tr"))
        LAG_A2.append(_Le_refexp("-Te*H*Te"))
        LAG_A2.append(_Le_refexp("-Tr*H*Tr"))
        LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
        LAG_A2.append(_Le_refexp("1/2*Tr*Tr*H"))
        if not no_occ:
            LAG_E.append(_refexp("(H*(To))"))

            LAG_A1.append(_Lo_refexp("(H-ECEPA)*To"))
            LAG_A1.append(_Lo_refexp("(H-ECEPA)*Ti"))
            LAG_A1.append(_Lo_refexp("(H-ECEPA)*Te"))
            LAG_A1.append(_Lo_refexp("(H-ECEPA)*Tr"))

            LAG_A1.append(_Lr_refexp("(H-ECEPA)*To"))

            LAG_A2.append(_Li_refexp("(H-ECEPA)*To"))

            LAG_A2.append(_Le_refexp("(H-ECEPA)*To"))

    elif singles == 7:
       LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
       LAG_E.append(_refexp("1/2*H*Te*Te"))
       LAG_E.append(_refexp("1/2*H*Ts*Ts"))

       LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ti"))
       LAG_A1.append(_Ls_refexp("H*Te"))
       LAG_A1.append(_Ls_refexp("H*Ts"))
       LAG_A1.append(_Ls_refexp("-Te*H"))
       LAG_A1.append(_Ls_refexp("-Ts*H"))
       LAG_A1.append(_Ls_refexp("1/2*H*Te*Te"))
       LAG_A1.append(_Ls_refexp("1/2*H*Ts*Ts"))
       LAG_A1.append(_Ls_refexp("-Te*H*Te"))
       LAG_A1.append(_Ls_refexp("-Ts*H*Ts"))
       LAG_A1.append(_Ls_refexp("1/2*Te*Te*H"))
       LAG_A1.append(_Ls_refexp("1/2*Ts*Ts*H"))

       LAG_A2.append(_Li_refexp("(H-ECEPA)*Ts"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))

       LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
       LAG_A2.append(_Le_refexp("H*Te"))
       LAG_A2.append(_Le_refexp("H*Ts"))
       LAG_A2.append(_Le_refexp("-Te*H"))
       LAG_A2.append(_Le_refexp("-Ts*H"))
       LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
       LAG_A2.append(_Le_refexp("1/2*H*Ts*Ts"))
       LAG_A2.append(_Le_refexp("-Te*H*Te"))
       LAG_A2.append(_Le_refexp("-Ts*H*Ts"))
       LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
       LAG_A2.append(_Le_refexp("1/2*Ts*Ts*H"))


LAG_E.set_rule()
LAG_A1.set_rule()
LAG_A2.set_rule()

#Replaces the Ti and Te strings, so that the equations can be parsed to the solver
if singles in [0,7]:

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

elif singles in [1,2,3,4,5,6]:

    if no_occ:

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',
             LABEL_IN:'FORM_MRCC_LAG_E',
             OP_LIST:['Ti','T2g','Te','T2g','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
             LABEL_IN:'FORM_MRCC_LAG_A1',
             OP_LIST:['Ti','T2g','Te','T2g','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
             LABEL_IN:'FORM_MRCC_LAG_A1',
             OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMr','LAM1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
             LABEL_IN:'FORM_MRCC_LAG_A2',
             OP_LIST:['Ti','T2g','Te','T2g','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
             LABEL_IN:'FORM_MRCC_LAG_A2',
             OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMr','LAM1']})
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})
    else:

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',
             LABEL_IN:'FORM_MRCC_LAG_E',
             OP_LIST:['Ti','T2g','Te','T2g','To','T1','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
             LABEL_IN:'FORM_MRCC_LAG_A1',
             OP_LIST:['Ti','T2g','Te','T2g','To','T1','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
             LABEL_IN:'FORM_MRCC_LAG_A1',
             OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMo','LAM1','LAMr','LAM1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
             LABEL_IN:'FORM_MRCC_LAG_A2',
             OP_LIST:['Ti','T2g','Te','T2g','To','T1','Tr','T1']})

        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
             LABEL_IN:'FORM_MRCC_LAG_A2',
             OP_LIST:['LAMi','LAM2g','LAMe','LAM2g','LAMo','LAM1','LAMr','LAM1']})
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
        PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})


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


if method in ['CEPT2','CCEPA','CEPA0']:
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
