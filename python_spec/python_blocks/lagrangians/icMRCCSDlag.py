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

#test-keyword
if method in ["CEPT2","CCEPA","TCPT2"]:
   singles = 0
if method in ["TCPT20","TCPT21","CEPT20","CEPT21","CEPT22","CEPT23","CEPT24","CEPT25","CEPT26","CEPT27","CCEPA0","CCEPA1","CCEPA2","CCEPA3","CCEPA4","CCEPA5","CCEPA6","CCEPA7"]:
   singles = int(method[5])
   print(singles)
   method_old = method
   method = method[0:5]  #"CCEPA" #THIS LINE NEEDS A REWORK!!!
#test-keyword

known_methods=["Yuri", "MRCCSD22", "CEPT2","CCEPA","TCPT2","BETACEPT2","CEPA0","CEPA0TeTi","PT2","BETACCEPA","SPLITCEPT2","SPLITCCEPA"]
if method not in known_methods :
    raise Exception(i_am+": unknown method:"+str(method))
print("Using the "+method+" method.")

    

# Select H0
if method in ['CEPT2','BETACEPT2','PT2','SPLITCEPT2','SPLITCCEPA','TCPT2']:
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

#---Test-for-new-keyword
if method in ["CEPT2","CCEPA","TCPT2"]:
#if method_old in ["CEPT20","CEPT21"]:
 #   known_singles=[0,1,2,3,4,5,6,7]
 #   singles = keywords.get('method.MP_P.singles')
    

#    if singles is None:
#       singles == 0

 #   if singles not in known_singles:
 #       raise Exception(i_am+": This single seperation is now known:"+str(singles))

#To -> Internal_Single_EXT
#Tr -> External_Single_EXT
    if singles == 1:
       print("i am defining To")
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
       print("me2")
       DEF_OP_FROM_OCC({LABEL:'To',
                        DESCR:to_shape})

       CLONE_OPERATOR({LABEL:'LAMo',
                       TEMPLATE:'To',
                       ADJOINT:True})
       print(to_shape)
       DEF_OP_FROM_OCC({LABEL:'Tr',
                        DESCR:tr_shape})

       CLONE_OPERATOR({LABEL:'LAMr',
                       TEMPLATE:'Tr',
                       ADJOINT:True})
#----------
#Define the Te and Ti operators
#with_singles_2
te_shape='PP,HV|PP,VV|PP,HH'
#--------------
#CIPT2 seperation ? "all excitations from the active space"
#te_shape='PP,VV|PV,VV|P,V'
#--------------

#with_singles_2   
ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'
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

#SPLIT_CEPT2_TEST

#if method in ["SPLITCEPT2","SPLITCCEPA"]:
#To -> Internal_Single_EXT
#Tr -> External_Single_EXT

#to_shape='V,H'      #1 
   #to_shape='P,V'      #2
   #to_shape='P,H'      #3
   #to_shape='V,H|P,V'  #4
   #to_shape='V,H|P,H'  #5
   #to_shape='P,V|P,H'  #6
#DEF_OP_FROM_OCC({LABEL:'To',
#                 DESCR:to_shape})
#CLONE_OPERATOR({LABEL:'LAMo',
#                TEMPLATE:'To',
#                ADJOINT:True})



#tr_shape='P,V|P,H'  #1
   #tr_shape='V,H|P,H'  #2
   #tr_shape='V,H|P,V'  #3
   #tr_shape='P,H'      #4
   #tr_shape='P,V'      #5
   #tr_shape='V,H'      #6 
#DEF_OP_FROM_OCC({LABEL:'Tr',
#                 DESCR:to_shape})
#CLONE_OPERATOR({LABEL:'LAMr',
#                TEMPLATE:'Tr',
#                ADJOINT:True})

#---------------  


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

if method in ['TCPT2','CEPT2','CCEPA','BETACEPT2','CEPA0TeTi','BETACCEPA','SPLITCEPT2','SPLITCCEPA']: #TeTistuff
   def _Le_refexp(x):
       return _refexp("LAMe(" + x + ")")
   
   def _Li_refexp(x):
       return _refexp("LAMi(" + x + ")")
   print("LAMe definiert")

if method in ['BETACEPT2','CEPA0TeTi','BETACCEPA']:  #THE Ts stuff
   def _Ls_refexp(x):
       return _refexp("LAMs(" + x + ")")
#keyword test
if singles in [0,7]:
   def _Ls_refexp(x):
       return _refexp("LAMs(" + x + ")")
#keyword test

if method in ["SPLITCEPT2","SPLITCCEPA"]:
   def _Lo_refexp(x):
       return _refexp("LAMo(" + x + ")")

   def _Lr_refexp(x):
       return _refexp("LAMr(" + x + ")")
#keyword test
if singles in [1,2,3,4,5,6]:
   def _Lo_refexp(x):
       return _refexp("LAMo(" + x + ")")

   def _Lr_refexp(x):
       return _refexp("LAMr(" + x + ")")
   print("LAMo definiert")
#keyword test

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("LAM2g(" + x + ")")


LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))

if method in ['BETACEPT2','CEPA0TeTi','BETACCEPA']: #Ts -> T1     Te,Ti -> T2g  --note:testidea1, check idea2

    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Ls_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
    LAG_A2.append(_Le_refexp("H"))

elif method in ['SPLITCEPT2','SPLITCCEPA']:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Lo_refexp("H"))
    LAG_A1.append(_Lr_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
    LAG_A2.append(_Le_refexp("H"))

#Keyword test
elif singles in [0,7]:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Ls_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
    LAG_A2.append(_Le_refexp("H"))

elif singles in [1,2,3,4,5,6]:
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _Lo_refexp("H"))
    LAG_A1.append(_Lr_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _Li_refexp("H"))
    LAG_A2.append(_Le_refexp("H"))


#keyword test


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
    

    if singles == 0: 
       LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
    
       if hamiltonian == "DYALL":
       
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))

          LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

          LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])"))
          LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

       else:
     
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  # exclude for  test
          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

    elif singles in [1,2,3,4,5,6]:

       LAG_E.append(_refexp("(H*(Te+Ti))+(H*(Tr+To))"))
       print("i am in singles1")
       if hamiltonian == "DYALL":

          print("i am in Dyall")
          LAG_A1.append(_Lo_refexp("(["+_h0_+",To])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Te])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Tr])"))

          LAG_A1.append(_Lr_refexp("(["+_h0_+",To])"))
          LAG_A1.append(_Lr_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Lr_refexp("(H-ECEPA)*Te"))
          LAG_A1.append(_Lr_refexp("(H-ECEPA)*Tr"))

          LAG_A2.append(_Li_refexp("(["+_h0_+",To])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Tr])"))

          LAG_A2.append(_Le_refexp("(["+_h0_+",To])"))
          LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Tr"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

       else:

          print("i am not in Dyall")
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*To"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

          LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*To"))
          LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Lr_refexp("(H-ECEPA)*Te"))
          LAG_A1.append(_Lr_refexp("(H-ECEPA)*Tr"))

          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*To"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test
          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Tr"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))


    elif singles == 7:
       LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
    
       if hamiltonian == "DYALL":
          LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))

          LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

          LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))
          LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

       else:
     
          LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))

          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

          LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))
          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

elif method == 'TCPT2':
    print("halli hallo")
    if singles == 0: 
       LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
       LAG_E.append(_refexp("1/2*H*Te*Te"))
    
       if hamiltonian == "DYALL":
       
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))

          LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

          LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])"))
          LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Le_refexp("H*Te"))
          LAG_A2.append(_Le_refexp("-Te*H"))
          LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
          LAG_A2.append(_Le_refexp("-Te*H*Te"))
          LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))

       else:
     
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  # exclude for  test
          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Le_refexp("H*Te"))
          LAG_A2.append(_Le_refexp("-Te*H"))
          LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
          LAG_A2.append(_Le_refexp("-Te*H*Te"))
          LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))

    if singles == 1: 
       LAG_E.append(_refexp("(H*(To+Tr))+(H*(Te+Ti))"))
       LAG_E.append(_refexp("1/2*H*Te*Te"))
       LAG_E.append(_refexp("1/2*H*Tr*Tr"))
    
       if hamiltonian == "DYALL":
       
          LAG_A1.append(_Lo_refexp("(["+_h0_+",To])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Ti])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Te])"))
          LAG_A1.append(_Lo_refexp("(["+_h0_+",Tr])"))

          LAG_A1.append(_Lr_refexp("(["+_h0_+",To])"))
          LAG_A1.append(_Lr_refexp("(["+_h0_+",Ti])"))
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

          LAG_A2.append(_Li_refexp("(["+_h0_+",To])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
          LAG_A2.append(_Li_refexp("(["+_h0_+",Tr])"))

          LAG_A2.append(_Le_refexp("(["+_h0_+",To])"))
          LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
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

       else:
     
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*To"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
          LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

          LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test
          LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
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

          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*To"))  
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
          LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test
          LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
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
elif method == 'BETACEPT2':
    #here goes cept2
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()
    
    LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))

    if hamiltonian == "DYALL":
       
       LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
       LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
       LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))

       LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
       LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
       LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))

       LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])"))
       LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))
       LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

    else:
     
       LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
       LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
       LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

       LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
       LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
       LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

       LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  # exclude for  test
       LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
       LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

#-------singles in Te------------
#    LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti")) #THIS IS WORKING RIGHT NOW, EQUATIONS SHOULD REFLECT THE CHANGE IN SINGLES!!!
#    LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))   # testfor with_singles_6
#    LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))   # 
#
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
#
#    LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
#    LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))
#    LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts"))   #test 


#-------------WORKING EQUATIONS--------
#    LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))

#    LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
#    LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
#    LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
#
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
#    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
#
#    LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))  # exclude for  test
#    LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
#    LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))

elif method == 'SPLITCEPT2':
    #here goes cept2
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

    LAG_E.append(_refexp("(H*(Tr+To))+(H*(Te+Ti))"))

    LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*To"))
    LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Ti")) 
    LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
    LAG_A1.append(_Lo_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

    LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*To"))
    LAG_A1.append(_Lr_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
    LAG_A1.append(_Lr_refexp("(H-ECEPA)*Te"))  
    LAG_A1.append(_Lr_refexp("(H-ECEPA)*Tr"))    

    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*To"))  
    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
    LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Tr"))

    LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*To"))  # exclude for  test
    LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*Te"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*Tr"))   #test 

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
       LAG_E.append(_refexp("(H*(To+Tr))+(H*(Te+Ti))"))
       LAG_E.append(_refexp("1/2*H*Te*Te"))
       LAG_E.append(_refexp("1/2*H*Tr*Tr"))

       LAG_A1.append(_Lo_refexp("(H-ECEPA)*To"))
       LAG_A1.append(_Lo_refexp("(H-ECEPA)*Ti"))
       LAG_A1.append(_Lo_refexp("(H-ECEPA)*Te"))
       LAG_A1.append(_Lo_refexp("(H-ECEPA)*Tr"))

       LAG_A1.append(_Lr_refexp("(H-ECEPA)*To"))
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

       LAG_A2.append(_Li_refexp("(H-ECEPA)*To"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))
       LAG_A2.append(_Li_refexp("(H-ECEPA)*Tr"))

       LAG_A2.append(_Le_refexp("(H-ECEPA)*To"))
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

elif method == 'BETACCEPA':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

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

#---Ts in ext

#    LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
#    LAG_E.append(_refexp("1/2*H*Te*Te"))
#    LAG_E.append(_refexp("1/2*H*Ts*Ts"))

#   LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ti"))
#   LAG_A1.append(_Ls_refexp("H*Te"))
#   LAG_A1.append(_Ls_refexp("H*Ts"))
#   LAG_A1.append(_Ls_refexp("-Te*H"))
#   LAG_A1.append(_Ls_refexp("-Ts*H"))
#   LAG_A1.append(_Ls_refexp("1/2*H*Te*Te"))
#   LAG_A1.append(_Ls_refexp("1/2*H*Ts*Ts"))
#   LAG_A1.append(_Ls_refexp("-Te*H*Te"))
#   LAG_A1.append(_Ls_refexp("-Ts*H*Ts"))
#   LAG_A1.append(_Ls_refexp("1/2*Te*Te*H"))
#   LAG_A1.append(_Ls_refexp("1/2*Ts*Ts*H"))

#   LAG_A2.append(_Li_refexp("(H-ECEPA)*Ts"))
#   LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
#   LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))

#   LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
#   LAG_A2.append(_Le_refexp("H*Te"))
#   LAG_A2.append(_Le_refexp("H*Ts"))
#   LAG_A2.append(_Le_refexp("-Te*H"))
#   LAG_A2.append(_Le_refexp("-Ts*H"))
#   LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
#   LAG_A2.append(_Le_refexp("1/2*H*Ts*Ts"))
#   LAG_A2.append(_Le_refexp("-Te*H*Te"))
#   LAG_A2.append(_Le_refexp("-Ts*H*Ts"))
#   LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
#   LAG_A2.append(_Le_refexp("1/2*Ts*Ts*H"))

elif method == 'SPLITCCEPA':
    DEF_SCALAR({LABEL:'ECEPA'})

    E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
    E_CEPA.set_rule()

    LAG_E.append(_refexp("(H*(To+Tr))+(H*(Te+Ti))"))
    LAG_E.append(_refexp("1/2*H*Te*Te"))
    LAG_E.append(_refexp("1/2*H*Tr*Tr"))

    LAG_A1.append(_Lo_refexp("(H-ECEPA)*To"))
    LAG_A1.append(_Lo_refexp("(H-ECEPA)*Ti"))
    LAG_A1.append(_Lo_refexp("(H-ECEPA)*Te"))
    LAG_A1.append(_Lo_refexp("(H-ECEPA)*Tr"))
    LAG_A1.append(_Lr_refexp("-Te*H"))
    LAG_A1.append(_Lr_refexp("-Tr*H"))
    LAG_A1.append(_Lr_refexp("1/2*H*Te*Te"))
    LAG_A1.append(_Lr_refexp("1/2*H*Tr*Tr"))
    LAG_A1.append(_Lr_refexp("-Te*H*Te"))
    LAG_A1.append(_Lr_refexp("-Tr*H*Tr"))
    LAG_A1.append(_Lr_refexp("1/2*Te*Te*H"))
    LAG_A1.append(_Lr_refexp("1/2*Tr*Tr*H"))

    LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
    LAG_A2.append(_Li_refexp("(H-ECEPA)*To"))
    LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))
    LAG_A2.append(_Li_refexp("(H-ECEPA)*Tr"))

    LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
    LAG_A2.append(_Le_refexp("(H-ECEPA)*To"))
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

LAG_E.set_rule()


if method in ['CEaPT2','CCaEPA']:
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

elif method in ['BETACEPT2','CEPA0TeTi','BETACCEPA']:

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


elif method in ['SPLITCEPT2','SPLITCCEPA']:
    LAG_A1.set_rule()
    LAG_A2.set_rule()
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_E', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

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

else:
    LAG_A1.set_rule()
    LAG_A2.set_rule()
    print("Rules set")
if singles in [0,7]:
    print("i am replacing Ts") 
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
    print("i am replacing To")
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
    PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})
    print("i am replacing To")
    print(LAG_A1)
    print("LAG_A2")
    print(LAG_A2)
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
    print("i am replacing To")
    print(LAG_A1)
    print("LAG_A2")
    print(LAG_A2)

#Replaces the Ti and Te strings, so that the equations can be parsed to the solver

#Make the Derivative with respect to LAM
print("JETZT LEITE ICH AB!! 1 1 11")
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


if method in ['CEPT2','CCEPA','CEPA0','CEPA0TeTi','BETACEPT2','BETACCEPA','SPLITCEPT2','SPLITCCEPA']:
   print("I AM DEFINING CEPA(0)")
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
   print("bi ba butzebar")
   OPTIMIZE({
           LABEL_OPT:'FOPT_MRCC_LAG',
           LABELS_IN:['FORM_MRCC_LAG_Amp2','FORM_MRCC_LAG_Amp1','FORM_MRCC_LAG_E']})


#-----
ref_relaxation.make_form_for_optref_minus3('FORM_MRCC_LAG_E', 'DEF_FORM_MRCC_LAG')
