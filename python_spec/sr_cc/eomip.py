#This is a python interface for GeCCo to calculate ionization energies using EOM theory

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

_orb = Orb_Info()

_env = keywords.env

mkDyson = keywords.is_keyword_set('calculate.ionization.Dyson')

# _multd2h will be required later to find out the spatial symmetry of what will be called 'R_q' ####

_multd2h = [[1,2,3,4,5,6,7,8],
[2,1,4,3,6,5,8,7],
[3,4,1,2,7,8,5,6],
[4,3,2,1,8,7,6,5],
[5,6,7,8,1,2,3,4],
[6,5,8,7,2,1,4,3],
[7,8,5,6,3,4,1,2],
[8,7,6,5,4,3,2,1]]

########### constructed _multd2h ############

### Now the necessary information will be obtained from the input and orb_info ###


_s2_0 = _orb.get('imult')
_isym_0 = _orb.get('lsym')
_ms_0 = _orb.get('ims')
_nsym = _orb.get('nsym')


if ((_ms_0 == 0) and ((_s2_0-1 % 4) == 0)):
    _msc_0 = 1
elif ((_ms_0 == 0) and ((_s2_0+1 % 4) == 0)):
    _msc_0 = -1
else:
    _msc_0 = 0

if _ms_0 != 0:
    quit_error('Must have closed shell reference in EOM-IP module!')

_sym = keywords.get('calculate.ionization.sym')

print("In EOM-IP module")

if (_sym == None):
    quit_error('Input missing: spatial symmetry of the excited state has not been set in the input')
else:
    _ncnt = len(_sym)

################################

#################################
# Defining the operators required
#################################

new_target('EOMIP_OPS')
depend('CC_OPS')

# read main info about method
minexc = keywords.get('method.CC.minexc')
minexc = int(minexc) if minexc is not None else 1
maxexc = int(keywords.get('method.CC.maxexc'))
maxexc = int(maxexc) if maxexc is not None else 2
truncate = keywords.get('method.CC.truncate')

DEF_EXCITATION({LABEL:'RIP',MIN_RANK:minexc,MAX_RANK:maxexc,CHARGE:-1})  # -1 = 1 electron less
CLONE_OPERATOR({LABEL:'LIP',TEMPLATE:'RIP',ADJOINT:True})
CLONE_OPERATOR({LABEL:'DIAG_IP',TEMPLATE:'RIP'})
CLONE_OPERATOR({LABEL:'DIAG_LIP',TEMPLATE:'LIP'})  # not sure that ordering of LIP and RIP is the same
CLONE_OPERATOR({LABEL:'A_RIP',TEMPLATE:'RIP'})
CLONE_OPERATOR({LABEL:'LIP_A',TEMPLATE:'RIP',ADJOINT:True})

DEF_SCALAR({LABEL:'L_IP'})
DEF_SCALAR({LABEL:'E_IP'})

DEF_SCALAR({LABEL:'NORM_IP'})

if mkDyson:
    new_target("EOMIP_DYSON")
    depend('SOLVE_LAMBDA')  # trigger LAMBDA equations
    # DYSON ORBITALS
    DEF_EXCITATION({LABEL:'DY_R',MIN_RANK:1,MAX_RANK:1,CHARGE:-1})  # -1 = 1 electron less
    CLONE_OPERATOR({LABEL:'DY_L',TEMPLATE:'DY_R',ADJOINT:True})

    # equation to evaluate

    # place holder for annihilation operator (adjoint is creation, obviously)
    DEF_EXCITATION({LABEL:'ANNI',MIN_RANK:1,MAX_RANK:1,CHARGE:-1})

    # generate equations with place holder
    form_dyR = stf.Formula("F_DY_R_0:LCC=<ANNI^+*RIP + LAM*ANNI^+*RIP>") # LCC is dummy
    form_dyR.set_rule()
    form_dyL = stf.Formula("F_DY_L_0:LCC=<LIP*ANNI>")
    form_dyL.set_rule()

    # remove the place holder
    DERIVATIVE({LABEL_RES:'F_DY_R',LABEL_IN:'F_DY_R_0',OP_RES:'DY_R',OP_DERIV:'ANNI^+'})
    DERIVATIVE({LABEL_RES:'F_DY_L',LABEL_IN:'F_DY_L_0',OP_RES:'DY_L',OP_DERIV:'ANNI'})



new_target("DIPOLE_OPS")

DEF_HAMILTONIAN({LABEL:'DM',MIN_RANK:1,MAX_RANK:1})
DEF_ME_LIST({LIST:'ME_DMX',OPERATOR:'DM',IRREP:2,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_DMY',OPERATOR:'DM',IRREP:3,'2MS':0,AB_SYM:+1})
DEF_ME_LIST({LIST:'ME_DMZ',OPERATOR:'DM',IRREP:1,'2MS':0,AB_SYM:+1})
IMPORT({LIST:'ME_DMX',TYPE:'XDIPLEN',ENV:_env})
IMPORT({LIST:'ME_DMY',TYPE:'YDIPLEN',ENV:_env})
IMPORT({LIST:'ME_DMZ',TYPE:'ZDIPLEN',ENV:_env})



new_target("EOMIP_EQS")
depend('EOMIP_OPS','CC_LAGRANGIAN','CC_EQS')

DERIVATIVE({LABEL_RES:'F_A_RIP',LABEL_IN:'F_CC_OMG',OP_RES:'A_RIP',OP_DERIV:'T',OP_MULT:'RIP'})

form_e_ip = stf.Formula("F_E_IP:E_IP=<RIP^+*A_RIP>")
form_e_ip.set_rule()



new_target("EOMIP_L_EQS")
depend('EOMIP_OPS','CC_LAGRANGIAN','CC_LAMBDA_EQS')

DERIVATIVE({LABEL_RES:'F_LIP_A',LABEL_IN:'F_CC_LAM_A',OP_RES:'LIP_A',OP_DERIV:'LAM',OP_MULT:'LIP'})

form_norm = stf.Formula("F_NORM_IP:NORM_IP=<LIP*RIP>")
form_norm.set_rule()


if mkDyson:
    new_target("EOMIP_DYTM")
    depend('DIPOLE_OPS','EOMIP_OPS')

    DEF_EXCITATION({LABEL:'DYTM_R',MIN_RANK:1,MAX_RANK:1,CHARGE:-1})  # -1 = 1 electron less
    CLONE_OPERATOR({LABEL:'DYTM_L',TEMPLATE:'DYTM_R',ADJOINT:True})

    form_dytmr = stf.Formula("F_DYTM_R_0:LCC=<ANNI^+*DM*RIP>")
    form_dytmr.append("<LAM*ANNI^+*DM*RIP>")
    form_dytmr.append("<ANNI^+*[DM,T]*RIP>")
    form_dytmr.append("<LAM*ANNI^+*[DM,T]*RIP>")
    # does that contribute?
    form_dytmr.append("<(1/2)*ANNI^+*[[DM,T],T]*RIP>")
    form_dytmr.append("<(1/2)*LAM*ANNI^+*[[DM,T],T]*RIP>")
    
    form_dytmr.set_rule()

    DERIVATIVE({LABEL_RES:'F_DYTM_R',LABEL_IN:'F_DYTM_R_0',OP_RES:'DYTM_R',OP_DERIV:'ANNI^+'})

    form_dytml = stf.Formula("F_DYTM_L_0:LCC=<LIP*DM*ANNI>")
    form_dytml.append("<LIP*[DM,T]*ANNI>")
    form_dytml.append("<(1/2)*LIP*[[DM,T],T]*ANNI>")

    form_dytml.set_rule()

    DERIVATIVE({LABEL_RES:'F_DYTM_L',LABEL_IN:'F_DYTM_L_0',OP_RES:'DYTM_L',OP_DERIV:'ANNI'})


#Here the iterations start for all symmetries

# ionize sym may appear several times (does it make sense?)
for _icnt in range (0,_ncnt):

    _sym_arr = _sym[_icnt]

    for _isym in range (0,_nsym):

        _no_root = int(_sym_arr[_isym])

        if( _no_root == 0):
            continue

        _msc = 0
        _ms_0 = 1

        _extension = '_' + str(_isym+1) + '_' + str(_msc)

        _list_rspns = 'LIST_RSPNS_' + _extension


### Defining the lists for all the operators involved

        new_target('EOMIP_LISTS'+_extension)
        depend('EOMIP_OPS','EOMIP_EQS')

        _op_list={'RIP':[_isym+1,_ms_0,_msc],
                  'A_RIP':[_isym+1,_ms_0,_msc],
                  'DIAG_IP':[_isym+1,_ms_0,_msc],
                  'E_IP':[1,0,0]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extension,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],MIN_REC:1,MAX_REC:_no_root})

        PRECONDITIONER({LIST_PRC:'ME_DIAG_IP'+_extension,
                        LIST_INP:'H0'})
        

### Lists needed for LHS equations:
        new_target('EOMIP_L_LISTS'+_extension)
        depend('EOMIP_OPS','EOMIP_EQS')

        # it seems that the internal logics requires LIP to have MS=-1
        _op_list={'LIP':[_isym+1,-_ms_0,_msc],
                  'LIP_A':[_isym+1,-_ms_0,_msc],
                  'DIAG_LIP':[_isym+1,-_ms_0,_msc],
                  'NORM_IP':[1,0,0]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extension,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],MIN_REC:1,MAX_REC:_no_root})

        PRECONDITIONER({LIST_PRC:'ME_DIAG_LIP'+_extension,
                        LIST_INP:'H0'})
        

### Optimising

        new_target('EOMIP_OPT_F'+_extension)

        depend('EOMIP_EQS','EOMIP_LISTS'+_extension)

        OPTIMIZE({LABEL_OPT:'FOPT_A_RIP'+_extension,
                  LABELS_IN:'F_A_RIP'})
        OPTIMIZE({LABEL_OPT:'FOPT_E_IP'+_extension,LABELS_IN:'F_E_IP'})
      
        if mkDyson:
            # add code to evaluate the dyson orbitals
            # stay on the same target, but add new dependencies:
            depend('EOMIP_DYSON')
            DEF_ME_LIST({LIST:'ME_DY_R'+_extension,OPERATOR:'DY_R',
                IRREP:_isym+1,'2MS':_ms_0,AB_SYM:0,MIN_REC:1,MAX_REC:_no_root})
            OPTIMIZE({LABEL_OPT:'FOPT_DY_R'+_extension,LABELS_IN:'F_DY_R'})

            if _isym != 3: ### to be adapted
                depend('EOMIP_DYTM')
                ASSIGN_ME2OP({LIST:'ME_DMZ',OPERATOR:'DM'})
                #if _isym == 0:
                #    ASSIGN_ME2OP({LIST:'ME_DMZ',OPERATOR:'DM'})
                #if _isym == 1:
                #    ASSIGN_ME2OP({LIST:'ME_DMX',OPERATOR:'DM'})
                #if _isym == 2:
                #    ASSIGN_ME2OP({LIST:'ME_DMY',OPERATOR:'DM'})

                DEF_ME_LIST({LIST:'ME_DYTM_R'+_extension,OPERATOR:'DYTM_R',
                    IRREP:_isym+1,'2MS':_ms_0,AB_SYM:0,MIN_REC:1,MAX_REC:_no_root})
                OPTIMIZE({LABEL_OPT:'FOPT_DYTM_R'+_extension,LABELS_IN:'F_DYTM_R'})


        new_target('EOMIP_L_OPT_F'+_extension)

        depend('EOMIP_L_EQS','EOMIP_L_LISTS'+_extension)

        OPTIMIZE({LABEL_OPT:'FOPT_LIP_A'+_extension,
                 LABELS_IN:'F_LIP_A'})
        OPTIMIZE({LABEL_OPT:'FOPT_NORM_IP'+_extension,LABELS_IN:'F_NORM_IP'})

        if mkDyson:
           depend('EOMIP_DYSON')
           DEF_ME_LIST({LIST:'ME_DY_L'+_extension,OPERATOR:'DY_L',
                IRREP:_isym+1,'2MS':_ms_0,AB_SYM:0,MIN_REC:1,MAX_REC:_no_root})
           OPTIMIZE({LABEL_OPT:'FOPT_DY_L'+_extension,LABELS_IN:'F_DY_L'})


#### Finally solving the eigenvalue equation to get the excitation energies

        new_target('EOMIP_SOLVE'+_extension)
        required()
        depend('EOMIP_OPT_F'+_extension,'EOMIP_LISTS'+_extension)

        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_RIP'+_extension]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_IP'+_extension]
        _solve_evp_basis[OP_MVP]=['A_RIP']
        _solve_evp_basis[OP_SVP]=['RIP']
        _solve_evp_basis[FORM]='FOPT_A_RIP'+_extension
        _solve_evp_basis[MODE]='DIA'
        _solve_evp_basis[N_ROOTS]=_no_root

        PRINT({STRING: ' ' })
        PRINT({STRING: '======================================= ' })
        PRINT({STRING: '= Calculating ionization in irrep:  ' + str(_isym+1)+' ='})
        PRINT({STRING: '======================================= ' })
        PRINT({STRING: ' ' })

        PRINT_MEL_INFO({LIST:'ME_RIP'+_extension})

        SOLVE_EVP(_solve_evp_basis)

        PRINT({STRING:''})
        PRINT({STRING:'Results:'})
        PRINT({STRING:'--------'})

        # make sure that this assignment holds
        ASSIGN_ME2OP({LIST:'ME_A_RIP'+_extension,OPERATOR:'A_RIP'})

        for i in range(1,_no_root +1):
            SET_STATE({LISTS:['ME_A_RIP'+_extension,
                              'ME_RIP'+_extension,
                              'ME_E_IP'+_extension] ,
                       ISTATE:i})
            EVALUATE({FORM:'FOPT_A_RIP'+_extension})

            EVALUATE({FORM:'FOPT_E_IP'+_extension})

            PRINT_MEL({LIST:'ME_E_IP'+_extension,COMMENT:'EOM IP STATE '+str(_isym+1)+'.'+str(i)+' excitation energy: ',FORMAT:'SCAL F24.14'})
            PUSH_RESULT({LIST:'ME_E_IP'+_extension,COMMENT:'EOM-IP_'+str(_isym+1)+'.'+str(i), FORMAT:"SCAL F20.14"})

        PRINT({STRING:''})
        PRINT({STRING:'Analysis of excitation vectors: '})
        ANALYZE_MEL({LISTS:'ME_RIP'+_extension,
                     LISTS_CV:'ME_RIP'+_extension})

        PRINT({STRING:''})
        if mkDyson:
            for i in range(1,_no_root+1):
                SET_STATE({LISTS:['ME_RIP'+_extension,'ME_DY_R'+_extension],ISTATE:i})
                EVALUATE({FORM:'FOPT_DY_R'+_extension})
                PRINT_MEL({LIST:'ME_DY_R'+_extension,
                    COMMENT:'Dyson orbital for state '+str(_isym+1)+'.'+str(i), 
                    FORMAT:'LIST'})
                PRINT({STRING:''})

            for i in range(1,_no_root+1):
                SET_STATE({LISTS:['ME_RIP'+_extension,'ME_DYTM_R'+_extension],ISTATE:i})
                EVALUATE({FORM:'FOPT_DYTM_R'+_extension})
                PRINT_MEL({LIST:'ME_DYTM_R'+_extension,
                    COMMENT:'Dyson dipole TMs for state '+str(_isym+1)+'.'+str(i),
                    FORMAT:'LIST'})
                PRINT({STRING:''})

        
### Left equations to be solved (if req.)
        new_target('EOMIP_LSOLVE'+_extension)
        if mkDyson:
            required()
        depend('EOMIP_SOLVE'+_extension,'EOMIP_L_OPT_F'+_extension,'EOMIP_L_LISTS'+_extension)

        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_LIP'+_extension]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_LIP'+_extension]
        _solve_evp_basis[OP_MVP]=['LIP_A']
        _solve_evp_basis[OP_SVP]=['LIP']
        _solve_evp_basis[FORM]='FOPT_LIP_A'+_extension
        _solve_evp_basis[MODE]='DIA'
        _solve_evp_basis[N_ROOTS]=_no_root

        PRINT({STRING: ' ' })
        PRINT({STRING: 'Calculating ionization (left eqs.) in irrep:  ' + str(_isym+1)})

        PRINT_MEL_INFO({LIST:'ME_LIP'+_extension})

        SOLVE_EVP(_solve_evp_basis)

        # have to renormalize such that <L|R> = 1
        for i in range(1,_no_root+1):
            SET_STATE({LISTS:['ME_LIP'+_extension,'ME_RIP'+_extension,
                'ME_NORM_IP'+_extension],ISTATE:i})
            # compute overlap
            EVALUATE({FORM:'FOPT_NORM_IP'+_extension})
            # and scale
            SCALE({LIST_RES:'ME_LIP'+_extension,LIST_INP:'ME_LIP'+_extension,
                   LIST_SCAL:'ME_NORM_IP'+_extension,FAC:1.0,INV:True})

        PRINT({STRING:''})
        PRINT({STRING:'Analysis of excitation vectors (biorth. norm): '})
        ANALYZE_MEL({LISTS:'ME_RIP'+_extension,
                     LISTS_CV:'ME_LIP'+_extension})
        PRINT({STRING:''})

        if mkDyson:
            for i in range(1,_no_root+1):
                SET_STATE({LISTS:['ME_LIP'+_extension,'ME_DY_L'+_extension],ISTATE:i})
                EVALUATE({FORM:'FOPT_DY_L'+_extension})
                PRINT_MEL({LIST:'ME_DY_L'+_extension,
                    COMMENT:'Dyson orbital (left) for state '+str(_isym+1)+'.'+str(i),
                    FORMAT:'LIST'})
                PRINT({STRING:''})



export_targets();
