"""

core part for setting up the different methods

"""
import python_interface.gecco_modules.string_to_form as stf
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *


#===================================================================================#
# helper routines
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


# Every term in the Lagrangian is enclosed by <C0^+ and C0>
def _refexp(x):
    return "<C0^+*(" + x + ")*C0>"

# The terms with the Lambda are always enclosed by <C0^+|LAM1 and C0> or <C0^+|LAM2g and C0>
def _L1_refexp(x):
    return _refexp("LAM1(" + x + ")")

def _L2_refexp(x):
    return _refexp("L2(" + x + ")")

def _L2g_refexp(x):
    return _refexp("LAM2g(" + x + ")")

def _L3g_refexp(x):
    return _refexp("LAM3g(" + x + ")")

def _L4g_refexp(x):
    return _refexp("LAM4g(" + x + ")")

def _Le_refexp(x):
    return _refexp("LAMe(" + x + ")")

def _Li_refexp(x):
    return _refexp("LAMi(" + x + ")")

def _Ls_refexp(x):
    return _refexp("LAMs(" + x + ")")

def _Lo_refexp(x):
    return _refexp("LAMo(" + x + ")")

def _Lr_refexp(x):
    return _refexp("LAMr(" + x + ")")



#===================================================================================#


###############################################################################################################################
def set_mrcc(nc_en,nc_rs,select,noVVV):
###############################################################################################################################

    if noVVV:
        T2_shape = 'VV,HH|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'  # skipped VVV amps here
    else:
        T2_shape = 'V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'

    # will later be replaced by "T2g" operators
    DEF_OP_FROM_OCC({LABEL:'T2',DESCR:T2_shape})
    CLONE_OPERATOR({LABEL:'L2',TEMPLATE:'T2',ADJOINT:True})

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
    #if nc_rs > 3:
    #    if not select:
    #        for nsingles in range(5):
    #            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
    #            for entryT in listT:
    #                print "Generating: "+_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]#+"]")
    #                LAG_A1.append(_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))###
    #
    ##        LAG_A1.append(_L1_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 4:
    #    if not select:
    #        for nsingles in range(6):
    #            listT = create_plist(nsingles,'T1',5-nsingles,'T2')
    #            for entryT in listT:
    #                print "Generating: "+_L1_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[#3]+"],"+entryT[4]+"]")
    #                LAG_A1.append(_L1_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"#+entryT[4]+"]"))
    ##        LAG_A1.append(_L1_refexp("1/120*[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 5:
    #    if not select:
    #        LAG_A1.append(_L1_refexp("1/720*[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    ## I think that singles can accomodate at most 6-fold
    ##if nc_rs > 6:
    ##    if not select:
    ##        LAG_A1.append(_L1_refexp("1/5040*[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    ##if nc_rs > 7:
    ##    if not select:
    ##        LAG_A1.append(_L1_refexp("1/40320*[[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))

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
    #    else:
    #        for nsingles in range(5):
    #            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
    #            for entryT in listT:
    #                print "Generating: "+_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]")
    #                LAG_A2.append(_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))

            #LAG_A2.append(_L2_refexp("1/24*[[[[H,T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 4:
    #    if not select:
    #        for nsingles in range(6):
    #            listT = create_plist(nsingles,'T1',5-nsingles,'T2')
    #            for entryT in listT:
    #                print "Generating: "+_L2_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]")
    #                LAG_A2.append(_L2_refexp("1/120*[[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"],"+entryT[4]+"]"))

    #       # LAG_A2.append(_L2_refexp("1/120*[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 5:
    #    if not select:
    #        LAG_A2.append(_L2_refexp("1/720*[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 6:
    #    if not select:
    #        LAG_A2.append(_L2_refexp("1/5040*[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))
    #if nc_rs > 7:
    #    if not select:
    #        LAG_A2.append(_L2_refexp("1/40320*[[[[[[[[H,T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2],T1+T2]"))

    PRINT({STRING:"Now expanding energy"})
    LAG_E.set_rule()
    PRINT({STRING:"Now expanding singles projection"})
    LAG_A1.set_rule()
    if nc_rs > 3 and not select:
        ngroups = 0
        for nsingles in range(5):
            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
            for entryT in listT:
                ngroups = ngroups+1
                print("Generating input for: "+_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))
                LAG_A1_C4 = stf.Formula("LAG_A1_C4_"+str(ngroups)+":MRCC_LAG_A1="+_L1_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))
                LAG_A1_C4.set_rule()
        groups = ["FORM_MRCC_LAG_A1"]
        for igrp in range(1,ngroups+1):
            groups.append("LAG_A1_C4_"+str(igrp))
        CONCAT({LABEL_RES:"FORM_MRCC_LAG_A1",LABEL_IN:groups})
        SUM_TERMS({LABEL_RES:"FORM_MRCC_LAG_A1",LABEL_IN:"FORM_MRCC_LAG_A1"})

    PRINT({STRING:"Now expanding doubles projection"})
    LAG_A2.set_rule()
    if nc_rs > 3 and not select:
        ngroups = 0
        for nsingles in range(5):
            listT = create_plist(nsingles,'T1',4-nsingles,'T2')
            for entryT in listT:
                ngroups = ngroups+1
                print("Generating input for: "+_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))
                LAG_A2_C4 = stf.Formula("LAG_A2_C4_"+str(ngroups)+":MRCC_LAG_A2="+_L2_refexp("1/24*[[[[H,"+entryT[0]+"],"+entryT[1]+"],"+entryT[2]+"],"+entryT[3]+"]"))
                LAG_A2_C4.set_rule()
        groups = ["FORM_MRCC_LAG_A2"]
        for igrp in range(1,ngroups+1):
            groups.append("LAG_A2_C4_"+str(igrp))
        CONCAT({LABEL_RES:"FORM_MRCC_LAG_A2",LABEL_IN:groups})
        SUM_TERMS({LABEL_RES:"FORM_MRCC_LAG_A2",LABEL_IN:"FORM_MRCC_LAG_A2"})


    REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',LABEL_IN:'FORM_MRCC_LAG_E',OP_LIST:['T2','T2g']})
    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',LABEL_IN:'FORM_MRCC_LAG_A1',OP_LIST:['T2','T2g']})
    REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',LABEL_IN:'FORM_MRCC_LAG_A2',OP_LIST:['T2','T2g','L2','LAM2g']})

###############################################################################################################################
#end set_mrcc
###############################################################################################################################

###############################################################################################################################
def set_mrcc_pt():
###############################################################################################################################
#   set the (T) correction
###############################################################################################################################
    depend('EVAL_HAM_D')  ### CHECK: check carefully that the Dyall was set up using C0 before any relaxation (ideally from C00)

    DEF_SCALAR({LABEL:'MRCC_EPT4'})
    DEF_SCALAR({LABEL:'MRCC_EPT5'})
    
    DEF_SCALAR({LABEL:'MRCC_LPT'})

    PT_E4 = stf.Formula('FORM_MRCC_PT_E4:MRCC_EPT4=<C0^+*T2g^+*([H,T3g])*C0>')
    PT_E4.set_rule()
    PT_E5 = stf.Formula('FORM_MRCC_PT_E5:MRCC_EPT5=<C0^+*T1^+*([H,T3g])*C0>')
    PT_E5.set_rule()

    PT_LAG = stf.Formula('FORM_MRCC_PT_LAG:MRCC_LPT=<C0^+*LAM3g*([HAM_D,T3g]+[H,T2g])*C0>')
    PT_LAG.set_rule()

    debug_FORM('FORM_MRCC_PT_E4',only_this=True)
    debug_FORM('FORM_MRCC_PT_E5',only_this=True)
    debug_FORM('FORM_MRCC_PT_LAG',only_this=True)



###############################################################################################################################
def set_mrcc_higher(maxexc,nc_en,nc_rs,select,noVVV):
###############################################################################################################################

    if not noVVV:
        raise Exception("Sorry, higher order MRCC only for CAS(2,2)")
        

    LAG_E = stf.Formula("FORM_MRCC_LAG_E:MRCC_LAG=" + _refexp("H"))
    LAG_A1 = stf.Formula("FORM_MRCC_LAG_A1:MRCC_LAG_A1=" + _L1_refexp("H"))
    LAG_A2 = stf.Formula("FORM_MRCC_LAG_A2:MRCC_LAG_A2=" + _L2g_refexp("H"))
    opT = 'T1+T2g'
    if (maxexc > 2):
        LAG_A3 = stf.Formula("FORM_MRCC_LAG_A3:MRCC_LAG_A3=" + _L3g_refexp("H"))
        opT += '+T3g'
    if (maxexc > 3):
        LAG_A4 = stf.Formula("FORM_MRCC_LAG_A4:MRCC_LAG_A4=" + _L3g_refexp("H"))
        opT += '+T4g'

    LAG_E.append(_refexp("[H,T1+T2g]"))
    if (maxexc>2):
        LAG_E.append(_refexp("[H,T3g]"))  # PPV,HHV can contribute
    if nc_en > 1:
       LAG_E.append(_refexp("1/2*[[H,T1+T2g],T1+T2g]"))
    if nc_en > 2:
       LAG_E.append(_refexp("1/6*[[[H,T1+T2g],T1+T2g],T1+T2g]"))
    if nc_en > 3:
       LAG_E.append(_refexp("1/24*[[[[H,T1+T2g],T1+T2g],T1+T2g],T1+T2g]"))


    LAG_A1.append(_L1_refexp("[H,"+opT+"]"))
    LAG_A2.append(_L2g_refexp("[H,"+opT+"]"))
    if (maxexc > 2):
        LAG_A3.append(_L3g_refexp("[H,"+opT+"]"))
    if (maxexc > 3):
        LAG_A4.append(_L3g_refexp("[H,"+opT+"]"))

    if nc_rs > 1:
        LAG_A1.append(_L1_refexp("0.5*[[H,"+opT+"]+"+opT+"]"))
        LAG_A2.append(_L2g_refexp("0.5*[[H,"+opT+"]+"+opT+"]"))
        if (maxexc > 2):
            LAG_A3.append(_L3g_refexp("0.5*[[H,"+opT+"]+"+opT+"]"))
        if (maxexc > 3):
            LAG_A4.append(_L4g_refexp("0.5*[[H,"+opT+"]+"+opT+"]"))

    if (nc_rs > 2 and not select):
        raise Exception("Higher-order CC with >2 commutators: Only for 'select' option!")

    if nc_rs > 2:
        LAG_A1.append(_L1_refexp("1/6*[[[H,T1],T1],T1]"))
        LAG_A2.append(_L2g_refexp("1/6*[[[H,T1],T1],T1]"))
        LAG_A2.append(_L2g_refexp("1/6*[[[H,T2g],T1],T1]"))
        LAG_A2.append(_L2g_refexp("1/6*[[[H,T1],T2g],T1]"))
        LAG_A2.append(_L2g_refexp("1/6*[[[H,T1],T1],T2g]"))   
        if maxexc > 2:
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T1],T1],T1]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T2g+T3g],T1],T1]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T1],T2g+T3g],T1]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T1],T1],T2g+T3g]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T2g],T2g],T1]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T2g],T1],T2g]"))
            LAG_A3.append(_L3g_refexp("1/6*[[[H,T1],T2g],T2g]"))
        if maxexc > 3:
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T2g+T3g+T4g],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T2g+T3g+T4g],T1]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T1],T2g+T3g+T4g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T2g],T2g],T1]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T2g],T1],T2g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T2g],T2g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T2g],T3g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T2g],T3g],T1]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T3g],T1],T2g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T2g],T1],T3g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T1],T3g],T2g]"))
            LAG_A4.append(_L4g_refexp("1/6*[[[H,T3g],T2g],T1]"))

    if nc_rs > 3:
        LAG_A2.append(_L2g_refexp("1/24*[[[[H,T1],T1],T1],T1]"))
        if maxexc > 2:
            LAG_A3.append(_L3g_refexp("1/24*[[[[H,T1],T1],T1],T1]"))
            LAG_A3.append(_L3g_refexp("1/24*[[[[H,T2g],T1],T1],T1]"))
            LAG_A3.append(_L3g_refexp("1/24*[[[[H,T1],T2g],T1],T1]"))
            LAG_A3.append(_L3g_refexp("1/24*[[[[H,T1],T1],T2g],T1]"))
            LAG_A3.append(_L3g_refexp("1/24*[[[[H,T1],T1],T1],T2g]"))
        if maxexc > 3: 
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T2g],T1],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T2g],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T2g],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T1],T2g]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T3g],T1],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T3g],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T3g],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T1],T3g]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T2g],T2g],T1],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T2g],T1],T2g],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T2g],T2g],T1]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T2g],T1],T1],T2g]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T2g],T1],T2g]"))
            LAG_A4.append(_L4g_refexp("1/24*[[[[H,T1],T1],T2g],T2g]"))
        
###############################################################################################################################
#end set_mrcc
###############################################################################################################################

    

###############################################################################################################################
def set_hybrids(method,separation,hamiltonian,singles,no_occ,noVVV):
###############################################################################################################################

    if hamiltonian != "":
        depend('H0') # not required, I think

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
            elif (no_occ and singles==6):
                DEF_OP_FROM_OCC({LABEL:'To',
                                 DESCR:tr_shape})

                CLONE_OPERATOR({LABEL:'LAMr',
                                TEMPLATE:'To',
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
    if separation == "SY":
        if (noVVV):
            te_shape='PP,HV|PP,VV|PP,HH'
            ti_shape='PV,HV|VV,HH|PV,HH|P,H'
            t2_shape=te_shape+'|'+ti_shape
        else:
            te_shape='PP,HV|PP,VV|PP,HH'
            ti_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'
            t2_shape=te_shape+'|'+ti_shape
    elif separation == "SY-INV":
        if (noVVV):
            ti_shape='PP,HV|PP,VV|PP,HH'
            te_shape='PV,HV|VV,HH|PV,HH|P,H'
            t2_shape=te_shape+'|'+ti_shape
        else:
            ti_shape='PP,HV|PP,VV|PP,HH'
            te_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'
            t2_shape=te_shape+'|'+ti_shape
    else:
        raise Exception("set_hybrid: unknown separation: "+str(separation))
        
    #Ti_shape='PP,HV|PP,VV|PP,HH' #needed for INVERSE
    #te_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H' #needed for INVERSE
    #inverse
    #ti_shape='PP,HV|PP,VV|PP,HH'
    #te_shape='VV,VH|PV,VV|PV,HV|VV,HH|PV,HH|V,H|P,V|P,H'
    
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

    #
    DEF_OP_FROM_OCC({LABEL:'T2',
                     DESCR:t2_shape})
    CLONE_OPERATOR({LABEL:'L2',
                    TEMPLATE:'T2',
                    ADJOINT:True})

    
    #-----------------TEST___TEST_for seperating singles (delete it and use T1 instead?)
    ts_shape='V,H|P,V|P,H'
    
    DEF_OP_FROM_OCC({LABEL:'Ts',
                     DESCR:ts_shape})
    CLONE_OPERATOR({LABEL:'LAMs',
                    TEMPLATE:'Ts',
                    ADJOINT:True})
    #----------------END_TEST
    

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

        LAG_E.append(_refexp("(H*T1)+(H*T2)"))
    
        if hamiltonian == 'DYALL':
    
            LAG_A1.append(_L1_refexp("(["+_h0_+",T1])"))
            LAG_A1.append(_L1_refexp("(["+_h0_+",T2])"))
    
            LAG_A2.append(_L2_refexp("(["+_h0_+",T1])"))
            LAG_A2.append(_L2_refexp("(["+_h0_+",T2])"))

        else:

            LAG_A1.append(_L1_refexp("("+_h0_+"-"+_h0exp_+")*T1"))
            LAG_A1.append(_L1_refexp("("+_h0_+"-"+_h0exp_+")*T2"))

            LAG_A2.append(_L2_refexp("("+_h0_+"-"+_h0exp_+")*T1"))
            LAG_A2.append(_L2_refexp("("+_h0_+"-"+_h0exp_+")*T2"))

    elif method=='Yuri':

        # a term-by-term setup for tests only
        
        #LAG_E.append(_refexp("[H,T1]"))
        #LAG_E.append(_refexp("[H,T2]"))
        #
        #LAG_A1.append(_L1_refexp("[H,T1]"))
        #LAG_A1.append(_L1_refexp("[H,T2]"))
        #LAG_A1.append(_L1_refexp("(1/2)*[[H,T1],T2]"))
        #
        #LAG_A2.append(_L2_refexp("[H,T1]"))
        #LAG_A2.append(_L2_refexp("[H,T2]"))
        #LAG_A2.append(_L2_refexp("(1/2)*[[H,T2],T2]"))
    
    
        LAG_E.append(_refexp("(H*T1)+(H*T2)"))
        LAG_E.append(_refexp("-(T1*H)-(T2*H)"))

        LAG_E.append(_refexp("(1/2)*(H*T1 *T1 )"))
        LAG_E.append(_refexp("(1/2)*(H*T1 *T2)"))
        LAG_E.append(_refexp("(1/2)*(H*T2*T1 )"))
        LAG_E.append(_refexp("(1/2)*(H*T2*T2)"))

        #LAG_E.append(_refexp("(1/2)*(H*T2'*T2'')"), avoid=["T2'","T2''"])
    
        LAG_E.append(_refexp("-(T1 *H*T1 )"))
        LAG_E.append(_refexp("-(T1 *H*T2)"))
        LAG_E.append(_refexp("-(T2*H*T1 )"))
        LAG_E.append(_refexp("-(T2*H*T2)"))
    
        LAG_E.append(_refexp("(1/2)*(T1 *T1 *H)"))
        LAG_E.append(_refexp("(1/2)*(T1 *T2*H)"))
        LAG_E.append(_refexp("(1/2)*(T2*T1 *H)"))
        LAG_E.append(_refexp("(1/2)*(T2*T2*H)"))
    
    

        LAG_A1.append(_L1_refexp("(H*T1)+(H*T2)"))
        LAG_A1.append(_L1_refexp("-(T1*H)-(T2*H)"))
    
        LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T1 )"))
        LAG_A1.append(_L1_refexp("(1/2)*(H*T1 *T2)"))
        LAG_A1.append(_L1_refexp("(1/2)*(H*T2*T1 )"))
        LAG_A1.append(_L1_refexp("(1/2)*(H*T2*T2)"))
    
        LAG_A1.append(_L1_refexp("-(T1 *H*T1 )"))
        LAG_A1.append(_L1_refexp("-(T1 *H*T2)"))
        LAG_A1.append(_L1_refexp("-(T2*H*T1 )"))
        LAG_A1.append(_L1_refexp("-(T2*H*T2)"))
    
        LAG_A1.append(_L1_refexp("(1/2)*(T1 *T1 *H)"))
        LAG_A1.append(_L1_refexp("(1/2)*(T1 *T2*H)"))
        LAG_A1.append(_L1_refexp("(1/2)*(T2*T1 *H)"))
        LAG_A1.append(_L1_refexp("(1/2)*(T2*T2*H)"))

    
    
        LAG_A2.append(_L2_refexp("(H*T1)+(H*T2)"))
        LAG_A2.append(_L2_refexp("-(T1*H)-(T2*H)"))

        LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T1 )"))
        LAG_A2.append(_L2_refexp("(1/2)*(H*T1 *T2)"))
        LAG_A2.append(_L2_refexp("(1/2)*(H*T2*T1 )"))
        LAG_A2.append(_L2_refexp("(1/2)*(H*T2*T2)"))

        LAG_A2.append(_L2_refexp("-(T1 *H*T1 )"))
        LAG_A2.append(_L2_refexp("-(T1 *H*T2)"))
        LAG_A2.append(_L2_refexp("-(T2*H*T1 )"))
        LAG_A2.append(_L2_refexp("-(T2*H*T2)"))

        LAG_A2.append(_L2_refexp("(1/2)*(T1 *T1 *H)"))
        LAG_A2.append(_L2_refexp("(1/2)*(T1 *T2*H)"))
        LAG_A2.append(_L2_refexp("(1/2)*(T2*T1 *H)"))
        LAG_A2.append(_L2_refexp("(1/2)*(T2*T2*H)"))

    elif method == 'CEPT2':
        #here goes cept2
        DEF_SCALAR({LABEL:'ECEPA'})
    
        E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
        E_CEPA.set_rule()
    
    
        if singles == 0:
            LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
            LAG_A2.append(_Le_refexp("(H-ECEPA)*Te")) #this term stays the same
    #        LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts")) #this term is only valid for INVERSE computations
    
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
    
                LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])")) #this term needs to be a comment for INVERSE
                LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

            else:
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
    
                LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts")) #this term needs to be a comment for INVERSE
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
    #        LAG_E.append(_refexp("1/2*H*Ts*Ts")) #INVERSE
    
            LAG_A2.append(_Le_refexp("H*Te"))
            LAG_A2.append(_Le_refexp("-Te*H"))
            LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
            LAG_A2.append(_Le_refexp("-Te*H*Te"))
            LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
    
    #        LAG_A2.append(_Le_refexp("H*Ts")) #INVERSE
    #        LAG_A2.append(_Le_refexp("-Ts*H")) #INVERSE
    #        LAG_A2.append(_Le_refexp("1/2*H*Ts*Ts")) #INVERSE
    #        LAG_A2.append(_Le_refexp("-Ts*H*Ts")) #INVERSE
    #        LAG_A2.append(_Le_refexp("1/2*Ts*Ts*H")) #INVERSE

            if hamiltonian == "DYALL": #exchanges the "PT2" terms for the DYALL Ham

                LAG_A1.append(_Ls_refexp("(["+_h0_+",Ts])"))
                LAG_A1.append(_Ls_refexp("(["+_h0_+",Ti])"))
                LAG_A1.append(_Ls_refexp("(["+_h0_+",Te])"))
    
                LAG_A2.append(_Li_refexp("(["+_h0_+",Ts])"))
                LAG_A2.append(_Li_refexp("(["+_h0_+",Ti])"))
                LAG_A2.append(_Li_refexp("(["+_h0_+",Te])"))
    
                LAG_A2.append(_Le_refexp("(["+_h0_+",Ts])")) #comment if INVERSE
                LAG_A2.append(_Le_refexp("(["+_h0_+",Ti])"))

            else:
    
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
    
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))

                LAG_A2.append(_Le_refexp("("+_h0_+"-"+_h0exp_+")*Ts")) #comment if INVERSE
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
    
                LAG_A1.append(_Ls_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
    
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ti"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Te"))
                LAG_A2.append(_Li_refexp("("+_h0_+"-"+_h0exp_+")*Ts"))
    
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
            LAG_A2.append(_Le_refexp("1/2*Ts*Ts*H"))

    elif method == 'CEPA0':
        #here goes CEPA0
        DEF_SCALAR({LABEL:'ECEPA'})
    
        E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
        E_CEPA.set_rule()
    
        LAG_E.append(_refexp("(H*T1)+(H*T2)"))
    
        LAG_A1.append(_L1_refexp("(H-ECEPA)*T1"))
        LAG_A1.append(_L1_refexp("(H-ECEPA)*T2"))
    
        LAG_A2.append(_L2_refexp("(H-ECEPA)*T1"))
        LAG_A2.append(_L2_refexp("(H-ECEPA)*T2"))
    elif method == 'CCEPA':
        #here goes the ccepa
        DEF_SCALAR({LABEL:'ECEPA'})

        E_CEPA=stf.Formula("FORM_ECEPA:ECEPA=<C0^+*H*C0>")
        E_CEPA.set_rule()
    
        if singles == 0:
           LAG_E.append(_refexp("(H*Ts)+(H*(Te+Ti))"))
           LAG_E.append(_refexp("1/2*H*Te*Te"))
    #       LAG_E.append(_refexp("1/2*H*Ts*Ts")) #INVERSE

           LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ts"))
           LAG_A1.append(_Ls_refexp("(H-ECEPA)*Ti"))
           LAG_A1.append(_Ls_refexp("(H-ECEPA)*Te"))
    
           LAG_A2.append(_Li_refexp("(H-ECEPA)*Ts"))
           LAG_A2.append(_Li_refexp("(H-ECEPA)*Ti"))
           LAG_A2.append(_Li_refexp("(H-ECEPA)*Te"))

           LAG_A2.append(_Le_refexp("(H-ECEPA)*Ts")) #comment if INVERSE
           LAG_A2.append(_Le_refexp("(H-ECEPA)*Ti"))
           LAG_A2.append(_Le_refexp("H*Te"))
           LAG_A2.append(_Le_refexp("-Te*H"))
           LAG_A2.append(_Le_refexp("1/2*H*Te*Te"))
           LAG_A2.append(_Le_refexp("-Te*H*Te"))
           LAG_A2.append(_Le_refexp("1/2*Te*Te*H"))
    #       LAG_A2.append(_Le_refexp("H*Ts")) #INVERSE
    #       LAG_A2.append(_Le_refexp("-Ts*H")) #INVERSE
    #       LAG_A2.append(_Le_refexp("1/2*H*Ts*Ts")) #INVERSE
    #       LAG_A2.append(_Le_refexp("-Ts*H*Ts")) #INVERSE
    #       LAG_A2.append(_Le_refexp("1/2*Ts*Ts*H")) #INVERSE
    
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
             OP_LIST:['L2','LAM2g','LAMi','LAM2g','LAMe','LAM2g','LAMs','LAM1']})

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
                 OP_LIST:['L2','LAM2g','LAMi','LAM2g','LAMe','LAM2g','LAMr','LAM1']})
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
                 OP_LIST:['L2','LAM2g','LAMi','LAM2g','LAMe','LAM2g','LAMo','LAM1','LAMr','LAM1']})
            PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A1', MODE:"SHORT"})
            PRINT_FORMULA({LABEL:'FORM_MRCC_LAG_A2', MODE:"SHORT"})

    else:
        REPLACE({LABEL_RES:'FORM_MRCC_LAG_E',
                LABEL_IN:'FORM_MRCC_LAG_E',
                OP_LIST:['T2','T2g','L2','LAM2g']})
        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A1',
                LABEL_IN:'FORM_MRCC_LAG_A1',
                OP_LIST:['T2','T2g','L2','LAM2g']})
        REPLACE({LABEL_RES:'FORM_MRCC_LAG_A2',
                LABEL_IN:'FORM_MRCC_LAG_A2',
                OP_LIST:['T2','T2g','L2','LAM2g']})


#end 
#########################################################################################    

    
