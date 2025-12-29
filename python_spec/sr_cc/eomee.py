#This is a python interface for GeCCo to calculate excitation energies using EOM theory
# (and LR is planned as well :) ... )

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *
import python_interface.gecco_modules.string_to_form as stf

_orb = Orb_Info()

_env = keywords.env

# multd2h will be required later to find out the spatial symmetry of what will be called 'R_q' ####

multd2h = [[1,2,3,4,5,6,7,8],
[2,1,4,3,6,5,8,7],
[3,4,1,2,7,8,5,6],
[4,3,2,1,8,7,6,5],
[5,6,7,8,1,2,3,4],
[6,5,8,7,2,1,4,3],
[7,8,5,6,3,4,1,2],
[8,7,6,5,4,3,2,1]]


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
    quit_error('Must have closed shell reference in EOM-EE / LR module!')

_sym = keywords.get('calculate.excitation.sym')
_mult = keywords.get('calculate.excitation.mult')

moments = keywords.is_keyword_set('calculate.excitation.moments')

# special switch for transition moments from ground state: TODO - more detailed control
moments0 = moments
export = True

use_z = keywords.is_keyword_set('calculate.excitation.z-vector')

use_z0 = use_z


print("In EOM-EE module")

if (_sym == None):
    quit_error('Input missing: spatial symmetry of the excited state has not been set in the input')
else:
    _ncnt = len(_sym)

if (_mult == None):
    print('Info: Multiplicity not specified - assuming singlet')
    _mult = []
    for ii in range(_ncnt):
        _mult.append(1)


################################

#################################
# Defining the operators required
#################################

new_target('EOMEE_OPS')
depend('CC_OPS')

# read main info about method
minexc = keywords.get('method.CC.minexc')
minexc = int(minexc) if minexc is not None else 1
maxexc = int(keywords.get('method.CC.maxexc'))
maxexc = int(maxexc) if maxexc is not None else 2
truncate = keywords.get('method.CC.truncate')

formalism=keywords.get('method.CC.xs_formalism')
formalism=formalism.upper() if formalism is not None else "EOM"

#DEF_EXCITATION({LABEL:'REX',MIN_RANK:minexc,MAX_RANK:maxexc}) 
CLONE_OPERATOR({LABEL:'REX',TEMPLATE:'T'})
CLONE_OPERATOR({LABEL:'LEX',TEMPLATE:'REX',ADJOINT:True})
CLONE_OPERATOR({LABEL:'LEX_T',TEMPLATE:'REX'})  # transpose of LEX, needed for analysis routine
CLONE_OPERATOR({LABEL:'DIAG_EX',TEMPLATE:'REX'})
CLONE_OPERATOR({LABEL:'DIAG_LEX',TEMPLATE:'LEX'})  # not sure that ordering of LEX and REX is the same
CLONE_OPERATOR({LABEL:'A_REX',TEMPLATE:'REX'})
CLONE_OPERATOR({LABEL:'LEX_A',TEMPLATE:'REX',ADJOINT:True})

DEF_SCALAR({LABEL:'L_EX'})
DEF_SCALAR({LABEL:'E_EX'})

# EOM formalism: reference state contribution to excited state
DEF_SCALAR({LABEL:'REX0'})

DEF_SCALAR({LABEL:'NORM_EX'})

DEF_HAMILTONIAN({LABEL:'OP1',MIN_RANK:1,MAX_RANK:1})
DEF_HAMILTONIAN({LABEL:'OP01',MIN_RANK:0,MAX_RANK:1})


new_target("EOMEE_EQS")
depend('EOMEE_OPS','CC_LAGRANGIAN','CC_EQS')

DERIVATIVE({LABEL_RES:'F_A_REX',LABEL_IN:'F_CC_OMG',OP_RES:'A_REX',OP_DERIV:'T',OP_MULT:'REX'})

# expression for energy
form_e_ex = stf.Formula("F_E_EX:E_EX=<REX^+*A_REX>")
form_e_ex.set_rule()

# expression for reference state contribution:
form_e_rex0 = stf.Formula("F_REX0:REX0=<H*REX>")
form_e_rex0.append("<[H,T]*REX>")
form_e_rex0.set_rule()
# Note: to get the actual REX0, it must be divided by the excitation energy
# For all non-totally symmetric excitations, a zero should be provided or the contribution is skipped

new_target("EOMEE_L_EQS")
depend('EOMEE_OPS','CC_LAGRANGIAN','CC_LAMBDA_EQS')

DERIVATIVE({LABEL_RES:'F_LEX_A',LABEL_IN:'F_CC_LAM_A',OP_RES:'LEX_A',OP_DERIV:'LAM',OP_MULT:'LEX'})

form_norm = stf.Formula("F_NORM_EX:NORM_EX=<LEX*REX>")
form_norm.set_rule()


new_target("EOMEE_DEN")
depend('EOMEE_OPS','CC_LAGRANGIAN', 'EOMEE_Z')

DEF_OP_FROM_OCC({LABEL:'CC_EOM_D',JOIN:2,DESCR:'H,;,H|,;P,H|H,P;,|,P;P,'})

# this is response-like
#DERIVATIVE({LABEL_RES:'F_EOM_D_0A',LABEL_IN:'F_CC_LAG',OP_RES:'LCC',OP_DERIV:'LAM',OP_MULT:'LEX'})
#DERIVATIVE({LABEL_RES:'F_EOM_D_0B',LABEL_IN:'F_EOM_D_0A',OP_RES:'LCC',OP_DERIV:'T',OP_MULT:'REX'})
#DERIVATIVE({LABEL_RES:'F_EOM_D',LABEL_IN:'F_EOM_D_0B',OP_RES:'CC_EOM_D',OP_DERIV:'H'})

# more EOM like:
form_exprop = stf.Formula("F_EOM_D_0B:LCC=<LEX*OP1*REX>")
form_exprop.append("<LEX*[OP1,T]*REX>")
form_exprop.append("<LEX*(1/2)*[[OP1,T],T]*REX>")  # is at most quadratic in T
form_exprop.set_rule()

DERIVATIVE({LABEL_RES:'F_EOM_D',LABEL_IN:'F_EOM_D_0B',OP_RES:'CC_EOM_D',OP_DERIV:'OP1'})

if use_z:
    CLONE_OPERATOR({LABEL:'CC_EOM_DZ',TEMPLATE:'CC_EOM_D'})
    form_expropZ = stf.Formula("F_EOM_DZ_0B:LCC=<LEX*OP1*REX>")
    form_expropZ.append("<LEX*[OP1,T]*REX>")
    form_expropZ.append("<LEX*(1/2)*[[OP1,T],T]*REX>")  
    form_expropZ.append("<Z*OP1>")  
    form_expropZ.append("<Z*[OP1,T]>")  
    form_expropZ.append("<Z*(1/2)*[[OP1,T],T]>")  
    form_expropZ.set_rule()

    DERIVATIVE({LABEL_RES:'F_EOM_DZ',LABEL_IN:'F_EOM_DZ_0B',OP_RES:'CC_EOM_DZ',OP_DERIV:'OP1'})


# transition moments from ground state
new_target("EOMEE_TM0F")
depend('EOMEE_OPS','CC_LAGRANGIAN')

DEF_OP_FROM_OCC({LABEL:'CC_EOM_TMF0',JOIN:2,DESCR:'H,;,H|,;P,H|H,P;,|,P;P,'})
DEF_OP_FROM_OCC({LABEL:'CC_EOM_TM0F',JOIN:2,DESCR:'H,;,H|,;P,H|H,P;,|,P;P,'})
# needed in some cases:
DEF_OP_FROM_OCC({LABEL:'CC_EOM_TM0F_SUM',JOIN:2,DESCR:'H,;,H|,;P,H|H,P;,|,P;P,'})

# f0 transition moment
# is the same for LR and EOM
form_tmf0 = stf.Formula("F_EOM_TMF0_0:LCC=<LEX*OP1>")
form_tmf0.append("<LEX*[OP1,T]>")
form_tmf0.append("<LEX*(1/2)*[[OP1,T],T]>")
form_tmf0.set_rule()

DERIVATIVE({LABEL_RES:'F_EOM_TMF0',LABEL_IN:'F_EOM_TMF0_0',OP_RES:'CC_EOM_TMF0',OP_DERIV:'OP1'})

# 0f transition moment
if formalism == "EOM":
   form_tm0f = stf.Formula("F_EOM_TM0F_0:LCC=<OP1*REX>")
   form_tm0f.append("<[OP1,T]*REX>") # no contribution
   form_tm0f.append("<(1/2)*[[OP1,T],T]*REX>") # no contribution
   form_tm0f.append("<LAM*OP1*REX>")
   form_tm0f.append("<LAM*[OP1,T]*REX>")
   form_tm0f.append("<LAM*(1/2)*[[OP1,T],T]*REX>")
   form_tm0f.set_rule()
else:
   form_tm0f = stf.Formula("F_EOM_TM0F_0:LCC=<OP1*REX>")
   form_tm0f.append("<LAM*[OP1,REX]>")
   form_tm0f.append("<LAM*[[OP1,T],REX]>")
   form_tm0f.append("<LAM*(1/2)*[[[OP1,T],T],REX]>")
   form_tm0f.set_rule()

DERIVATIVE({LABEL_RES:'F_EOM_TM0F',LABEL_IN:'F_EOM_TM0F_0',OP_RES:'CC_EOM_TM0F',OP_DERIV:'OP1'})

if formalism=="EOM":
   # R0 contribution to 0f transition moment:
   DEF_OP_FROM_OCC({LABEL:'CC_EOM_TM0F_R0',JOIN:2,DESCR:',;,|H,;,H|,;P,H|H,P;,|,P;P,'})
   form_tm0f_r0 = stf.Formula("F_EOM_TM0F_R0_0:LCC=<OP01*REX0>")
   form_tm0f_r0.append("<[OP01,T]*REX0>")
   form_tm0f_r0.append("<(1/2)*[[OP01,T],T]*REX0>")
   form_tm0f_r0.append("<LAM*OP01*REX0>")
   form_tm0f_r0.append("<LAM*[OP01,T]*REX0>")
   form_tm0f_r0.append("<LAM*(1/2)*[[OP01,T],T]*REX0>")
   form_tm0f_r0.set_rule()

   DERIVATIVE({LABEL_RES:'F_EOM_TM0F_R0',LABEL_IN:'F_EOM_TM0F_R0_0',OP_RES:'CC_EOM_TM0F_R0',OP_DERIV:'OP01'})

if formalism=="LR":
    depend("LREE_M")
    DEF_OP_FROM_OCC({LABEL:'CC_EOM_TM0F_M',JOIN:2,DESCR:'H,;,H|,;P,H|H,P;,|,P;P,'})
    form_tm0f_m = stf.Formula("F_EOM_TM0F_M_0:LCC=<M*OP1>")
    form_tm0f_m.append("<M*[OP1,T]>")
    form_tm0f_m.append("<M*(1/2)*[[OP1,T],T]>")
    form_tm0f_m.set_rule()

    DERIVATIVE({LABEL_RES:'F_EOM_TM0F_M',LABEL_IN:'F_EOM_TM0F_M_0',OP_RES:'CC_EOM_TM0F_M',OP_DERIV:'OP1'})


new_target("LREE_M")
depend('EOMEE_OPS','EOMEE_L_EQS','CC_LAMBDA_EQS')

CLONE_OPERATOR({LABEL:'MRHS',TEMPLATE:'LEX'})
CLONE_OPERATOR({LABEL:'M',TEMPLATE:'LEX'})
CLONE_OPERATOR({LABEL:'M_A',TEMPLATE:'LEX'})

# transformation
DERIVATIVE({LABEL_RES:'F_LR_M_A',LABEL_IN:'F_CC_LAM_A',OP_RES:'M_A',OP_DERIV:'LAM',OP_MULT:'M'})

# rhs
DERIVATIVE({LABEL_RES:'F_LR_MRHS0',LABEL_IN:'F_CC_LAG',OP_RES:'NORM_EX',OP_DERIV:'T',OP_MULT:'REX'})  # NORM_EX is dummy
DERIVATIVE({LABEL_RES:'F_LR_MRHS',LABEL_IN:'F_LR_MRHS0',OP_RES:'MRHS',OP_DERIV:'T'})


new_target("EOMEE_Z")
depend('EOMEE_OPS','EOMEE_L_EQS','CC_LAMBDA_EQS')

CLONE_OPERATOR({LABEL:'ZRHS',TEMPLATE:'LEX'})
CLONE_OPERATOR({LABEL:'Z',TEMPLATE:'LEX'})
CLONE_OPERATOR({LABEL:'Z_A',TEMPLATE:'LEX'})

# for the transformation:
DERIVATIVE({LABEL_RES:'F_Z_A',LABEL_IN:'F_CC_LAM_A',OP_RES:'Z_A',OP_DERIV:'LAM',OP_MULT:'Z'})

# for the rhs:
# this is in principle the modified Lagrangian, but we also need to include the GS derivative:
form_zrhs = stf.Formula("F_EOM_ZRHS0:LCC=<LEX*H*REX>")
form_zrhs.append("<LEX*[H,T]*REX>")
form_zrhs.append("<LEX*(1/2)*[[H,T],T]*REX>")
form_zrhs.append("<LEX*(1/6)*[[[H,T],T],T]*REX>")
form_zrhs.append("<LEX*(1/24)*[[[[H,T],T],T],T]*REX>")
form_zrhs.set_rule()

DERIVATIVE({LABEL_RES:'F_EOM_ZRHS',LABEL_IN:'F_EOM_ZRHS0',OP_RES:'ZRHS',OP_DERIV:'T'})



#Here the iterations start for all symmetries

#  sym may appear several times (does it make sense?)
for _icnt in range (0,_ncnt):

    _sym_arr = _sym[_icnt]
    _s2 =int(_mult[_icnt])

    _ms = 0
    if ((_ms == 0) and (_s2  == 1)):
        _msc = 1
    elif ((_ms == 0) and (_s2 == 3)):
        _msc = -1
    else:
        quit_error('EOMEE knows only singlets and triplets, found: (2S+1) = '+str(_s2))

    for _isym in range (0,_nsym):

        n_root = int(_sym_arr[_isym])

        if( n_root == 0):
            continue

        _ms_0 = 0

        _extension = '_' + str(_isym+1) + '_' + str(_msc)

        _list_rspns = 'LIST_RSPNS_' + _extension


### Defining the lists for all the operators involved

        new_target('EOMEE_LISTS'+_extension)
        depend('EOMEE_OPS','EOMEE_EQS')

        _op_list={'REX':[_isym+1,_ms_0,_msc],
                  'A_REX':[_isym+1,_ms_0,_msc],
                  'DIAG_EX':[_isym+1,_ms_0,_msc],
                  'E_EX':[1,0,0]}

        if _isym==0 and _msc==1:
            _op_list['REX0']=[1,0,0]

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extension,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],MIN_REC:1,MAX_REC:n_root})

        PRECONDITIONER({LIST_PRC:'ME_DIAG_EX'+_extension,
                        LIST_INP:'H0'})
        

### Lists needed for LHS equations:
        new_target('EOMEE_L_LISTS'+_extension)
        depend('EOMEE_OPS','EOMEE_EQS')

        _op_list={'LEX':[_isym+1,_ms_0,_msc],
                  'LEX_T':[_isym+1,_ms_0,_msc],
                  'LEX_A':[_isym+1,_ms_0,_msc],
                  'DIAG_LEX':[_isym+1,-_ms_0,_msc],
                  'NORM_EX':[1,0,0]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extension,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],MIN_REC:1,MAX_REC:n_root})

        PRECONDITIONER({LIST_PRC:'ME_DIAG_LEX'+_extension,
                        LIST_INP:'H0'})
        

### Optimising

        new_target('EOMEE_OPT_F'+_extension)

        depend('EOMEE_EQS','EOMEE_LISTS'+_extension)

        OPTIMIZE({LABEL_OPT:'FOPT_A_REX'+_extension,
                  LABELS_IN:'F_A_REX'})
        OPTIMIZE({LABEL_OPT:'FOPT_E_EX'+_extension,LABELS_IN:'F_E_EX'})

        if _isym==0 and _msc==1:
            OPTIMIZE({LABEL_OPT:'FOPT_REX0'+_extension,LABELS_IN:'F_REX0'})

      
        new_target('EOMEE_L_OPT_F'+_extension)

        depend('EOMEE_L_EQS','EOMEE_L_LISTS'+_extension)

        ASSIGN_ME2OP({LIST:'ME_REX'+_extension,OPERATOR:'REX'})
        ASSIGN_ME2OP({LIST:'ME_LEX'+_extension,OPERATOR:'LEX'})
        ASSIGN_ME2OP({LIST:'ME_LEX_A'+_extension,OPERATOR:'LEX_A'})

        OPTIMIZE({LABEL_OPT:'FOPT_LEX_A'+_extension,
                 LABELS_IN:'F_LEX_A'})
        OPTIMIZE({LABEL_OPT:'FOPT_NORM_EX'+_extension,LABELS_IN:'F_NORM_EX'})


#### Finally solving the eigenvalue equation to get the excitation energies

        new_target('EOMEE_SOLVE'+_extension)
        required()
        depend('EOMEE_OPT_F'+_extension,'EOMEE_LISTS'+_extension)

        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_REX'+_extension]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_EX'+_extension]
        _solve_evp_basis[OP_MVP]=['A_REX']
        _solve_evp_basis[OP_SVP]=['REX']
        _solve_evp_basis[FORM]='FOPT_A_REX'+_extension
        _solve_evp_basis[MODE]='DIA'
        _solve_evp_basis[N_ROOTS]=n_root

        PRINT({STRING: ' ' })
        PRINT({STRING: '======================================= ' })
        PRINT({STRING: '= Calculating excitation in irrep:  ' + str(_isym+1)+' ='})
        PRINT({STRING: '======================================= ' })
        PRINT({STRING: ' ' })

        PRINT_MEL_INFO({LIST:'ME_REX'+_extension})

        SOLVE_EVP(_solve_evp_basis)

        PRINT({STRING:''})
        PRINT({STRING:'Results:'})
        PRINT({STRING:'--------'})

        # make sure that this assignment holds
        ASSIGN_ME2OP({LIST:'ME_A_REX'+_extension,OPERATOR:'A_REX'})

        for i in range(1,n_root +1):
            SET_STATE({LISTS:['ME_A_REX'+_extension,
                              'ME_REX'+_extension,
                              'ME_E_EX'+_extension] ,
                       ISTATE:i})
            EVALUATE({FORM:'FOPT_A_REX'+_extension})

            EVALUATE({FORM:'FOPT_E_EX'+_extension})

            if _isym==0 and _msc==1:
                SET_STATE({LISTS:['ME_REX0'+_extension],
                           ISTATE:i})
                EVALUATE({FORM:'FOPT_REX0'+_extension})
                SCALE({LIST_RES:'ME_REX0'+_extension,LIST_INP:'ME_REX0'+_extension,
                       LIST_SCAL:'ME_E_EX'+_extension,FAC:1.0,INV:True})


            PRINT_MEL({LIST:'ME_E_EX'+_extension,COMMENT:'EOM EE STATE '+str(_isym+1)+'.'+str(i)+' excitation energy: ',FORMAT:'SCAL F24.14'})
            PUSH_RESULT({LIST:'ME_E_EX'+_extension,COMMENT:'EOM-EE_'+str(_isym+1)+'.'+str(i), FORMAT:"SCAL F20.14"})
            if _isym==0 and _msc==1 and formalism=="EOM":
                PRINT_MEL({LIST:'ME_REX0'+_extension,COMMENT:'                             <0|R> = ',FORMAT:'SCAL F24.14'})

        PRINT({STRING:''})
        PRINT({STRING:'Analysis of excitation vectors: '})
        ANALYZE_MEL({LISTS:'ME_REX'+_extension,
                     LISTS_CV:'ME_REX'+_extension})


        
### Left equations to be solved (if req.)
        new_target('EOMEE_LSOLVE'+_extension)
        depend('EOMEE_SOLVE'+_extension,'EOMEE_L_OPT_F'+_extension,'EOMEE_L_LISTS'+_extension)

        # copy REX to LEX as start vectors:
        for i in range(1,n_root+1):
            SET_STATE({LISTS:['ME_LEX'+_extension,'ME_REX'+_extension],ISTATE:i})
            COPY_LIST({LIST_RES:'ME_LEX'+_extension,LIST_INP:'ME_REX'+_extension,FAC:1.0,ADJOINT:True})

        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_LEX'+_extension]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_LEX'+_extension]
        _solve_evp_basis[OP_MVP]=['LEX_A']
        _solve_evp_basis[OP_SVP]=['LEX']
        _solve_evp_basis[FORM]='FOPT_LEX_A'+_extension
        _solve_evp_basis[MODE]='DIA'
        _solve_evp_basis[N_ROOTS]=n_root
        _solve_evp_basis[INIT]=True

        PRINT({STRING: ' ' })
        PRINT({STRING: 'Calculating excitations (left eqs.) in irrep:  ' + str(_isym+1)})

        PRINT_MEL_INFO({LIST:'ME_LEX'+_extension})

        SOLVE_EVP(_solve_evp_basis)

        # have to renormalize such that <L|R> = 1
        for i in range(1,n_root+1):
            SET_STATE({LISTS:['ME_LEX'+_extension,'ME_LEX_T'+_extension,'ME_REX'+_extension,
                'ME_NORM_EX'+_extension],ISTATE:i})
            # compute overlap
            EVALUATE({FORM:'FOPT_NORM_EX'+_extension})
            # and scale
            SCALE({LIST_RES:'ME_LEX'+_extension,LIST_INP:'ME_LEX'+_extension,
                    LIST_SCAL:'ME_NORM_EX'+_extension,FAC:1.0,INV:True})
      
            # for the analysis we have to write a transposed list:
            COPY_LIST({LIST_RES:'ME_LEX_T'+_extension,LIST_INP:'ME_LEX'+_extension,FAC:1.0,ADJOINT:True})


        PRINT({STRING:''})
        PRINT({STRING:'Analysis of excitation vectors (biorth. norm): '})
        ANALYZE_MEL({LISTS:'ME_REX'+_extension,
                     LISTS_CV:'ME_LEX_T'+_extension})
        PRINT({STRING:''})


# do here a loop to compute all moments from ground to excited states:
if moments0:
    for _icnt in range (0,_ncnt):

        _sym_arr = _sym[_icnt]
        _s2 =int(_mult[_icnt])

        _ms = 0
        if ((_ms == 0) and (_s2  == 1)):
            _msc = 1
        elif ((_ms == 0) and (_s2 == 3)):
            _msc = -1
        else:
            quit_error('EOMEE knows only singlets and triplets, found: (2S+1) = '+str(_s2))

        for _isym in range (0,_nsym):

            n_root = int(_sym_arr[_isym])

            if( n_root == 0):
                continue

            _ms_0 = 0

            _extension = '_' + str(_isym+1) + '_' + str(_msc)    

            new_target('EOMEE_TM0'+_extension)
            required()

            depend('EOMEE_TM0F','EOMEE_LSOLVE'+_extension,'EOMEE_SOLVE'+_extension)

            n_dens = n_root

            PRINT({STRING:f'Transition moments in {formalism} formalism:'})

            DEF_ME_LIST({LIST:'ME_CC_EOM_TM0F'+_extension,OPERATOR:'CC_EOM_TM0F',
                IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
            DEF_ME_LIST({LIST:'ME_CC_EOM_TMF0'+_extension,OPERATOR:'CC_EOM_TMF0',
                IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})

            if _isym==0 and _msc==1 and formalism=="EOM":
                DEF_ME_LIST({LIST:'ME_CC_EOM_TM0F_R0'+_extension,OPERATOR:'CC_EOM_TM0F_R0',
                    IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
            
            ASSIGN_ME2OP({LIST:'ME_LEX'+_extension,OPERATOR:'LEX'})
            ASSIGN_ME2OP({LIST:'ME_REX'+_extension,OPERATOR:'REX'})


            if formalism=="LR":
                depend("LREE_M")

                DEF_ME_LIST({LIST:'ME_CC_EOM_TM0F_M'+_extension,OPERATOR:'CC_EOM_TM0F_M',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                DEF_ME_LIST({LIST:'ME_MRHS'+_extension,OPERATOR:'MRHS',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                DEF_ME_LIST({LIST:'ME_M'+_extension,OPERATOR:'M',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                for istate in range(n_root):
                    SET_STATE({LISTS:['ME_M'+_extension,'ME_E_EX'+_extension],ISTATE:istate+1})
                    SET_FREQ({LIST:'ME_M'+_extension,LIST_FREQ:'ME_E_EX'+_extension})
                    PRINT_MEL({LIST:'ME_E_EX'+_extension,COMMENT:'Shift for root '+str(_isym+1)+'.'+str(i)+': ',FORMAT:'SCAL F24.14'})
                DEF_ME_LIST({LIST:'ME_M_A'+_extension,OPERATOR:'M_A',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                DEF_ME_LIST({LIST:'ME_DIAG_M'+_extension,OPERATOR:'DIAG_LEX',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1})

                OPTIMIZE({LABEL_OPT:'FOPT_LR_M'+_extension,LABELS_IN:['F_LR_MRHS','F_LR_M_A']})
                #OPTIMIZE({LABEL_OPT:'FOPT_LR_TM0F_M'+_extension,LABELS_IN:'F_EOM_TM0F_M'})

                PRECONDITIONER({LIST_PRC:'ME_DIAG_M'+_extension,
                        LIST_INP:'H0'})

                PRINT({STRING:'Solving M equations for symmetry '+str(_isym+1)})
                SOLVE_LEQ({LIST_OPT:'ME_M'+_extension,
                        LIST_PRC:'ME_DIAG_M'+_extension,
                        OP_MVP:'M_A',
                        OP_SVP:'M',
                        OP_RHS:'MRHS',
                        FORM:'FOPT_LR_M'+_extension,
                        MODE:'DIA',
                        N_ROOTS:n_dens})

            if formalism=="LR" or (formalism=="EOM" and _isym==0 and _msc==1):
                DEF_ME_LIST({LIST:'ME_CC_EOM_TM0F_SUM'+_extension,OPERATOR:'CC_EOM_TM0F_SUM',
                        IRREP:_isym+1,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})

            label_list = ['F_EOM_TMF0','F_EOM_TM0F']
            if _isym==0 and _msc==1 and formalism=="EOM":
                label_list.append('F_EOM_TM0F_R0')
            if formalism=="LR":
                label_list.append('F_EOM_TM0F_M')
            OPTIMIZE({LABEL_OPT:'FOPT_EOM_TM0F'+_extension,LABELS_IN:label_list})


            for istate in range(n_root):
                SET_STATE({LISTS:'ME_LEX'+_extension,ISTATE:istate+1})
                SET_STATE({LISTS:'ME_REX'+_extension,ISTATE:istate+1})

                SET_STATE({LISTS:'ME_CC_EOM_TM0F'+_extension,ISTATE:istate+1})
                SET_STATE({LISTS:'ME_CC_EOM_TMF0'+_extension,ISTATE:istate+1})

                if  _isym==0 and _msc==1 and formalism=="EOM":
                    SET_STATE({LISTS:'ME_REX0'+_extension,ISTATE:istate+1})
                    SET_STATE({LISTS:'ME_CC_EOM_TM0F_R0'+_extension,ISTATE:istate+1})
                    SET_STATE({LISTS:'ME_CC_EOM_TM0F_SUM'+_extension,ISTATE:istate+1})
                if formalism=="LR":
                    SET_STATE({LISTS:'ME_M'+_extension,ISTATE:istate+1})
                    SET_STATE({LISTS:'ME_CC_EOM_TM0F_M'+_extension,ISTATE:istate+1})
                    SET_STATE({LISTS:'ME_CC_EOM_TM0F_SUM'+_extension,ISTATE:istate+1})

                EVALUATE({FORM:'FOPT_EOM_TM0F'+_extension})

                PRINT({STRING:''})

                PRINT({STRING:'0F transition moment (ground state to state '+str(_isym+1)+'.'+str(istate+1)+')'})
                EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TM0F'+_extension,ADD_REF:False})

                if  _isym==0 and _msc==1 and formalism=="EOM":
                    PRINT({STRING:'0F transition moment (ground state to state '+str(_isym+1)+'.'+str(istate+1)+') - R0 contr.'})
                    EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TM0F_R0'+_extension,ADD_REF:False})
                    ADD({LIST_SUM:'ME_CC_EOM_TM0F_SUM'+_extension,
                             LISTS:['ME_CC_EOM_TM0F'+_extension,'ME_CC_EOM_TM0F_R0'+_extension],FAC:[1.0,1.0]})
                    PRINT({STRING:'0F transition moment (ground state to state '+str(_isym+1)+'.'+str(istate+1)+') - sum'})
                    EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TM0F_SUM'+_extension,ADD_REF:False})

                if formalism=="LR":
                    PRINT({STRING:'0F transition moment (ground state to state '+str(_isym+1)+'.'+str(istate+1)+') - M contr.'})
                    EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TM0F_M'+_extension,ADD_REF:False})
                    ADD({LIST_SUM:'ME_CC_EOM_TM0F_SUM'+_extension,
                         LISTS:['ME_CC_EOM_TM0F'+_extension,'ME_CC_EOM_TM0F_M'+_extension],FAC:[1.0,1.0]})
                    PRINT({STRING:'0F transition moment (ground state to state '+str(_isym+1)+'.'+str(istate+1)+') - sum'})
                    EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TM0F_SUM'+_extension,ADD_REF:False})

                PRINT({STRING:'F0 transition moment (state '+str(_isym+1)+'.'+str(istate+1)+' to ground state)'})
                EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_TMF0'+_extension,ADD_REF:False})

                if export:
                    if formalism=="EOM" and _isym==0 and _msc==1:
                        EXPORT_DAO({RANK:1,DENS:'ME_CC_EOM_TM0F_SUM'+_extension,
                            OUTPUT:'export_tm0f',ADD_REF:False})
                    elif formalism=="EOM":
                        EXPORT_DAO({RANK:1,DENS:'ME_CC_EOM_TM0F'+_extension,
                            OUTPUT:'export_tm0f',ADD_REF:False})
                    else:
                        EXPORT_DAO({RANK:1,DENS:'ME_CC_EOM_TM0F_SUM'+_extension,
                            OUTPUT:'export_tm0f',ADD_REF:False})
                    EXPORT_DAO({RANK:1,DENS:'ME_CC_EOM_TMF0'+_extension,
                        OUTPUT:'export_tmf0',ADD_REF:False})

            PRINT({STRING:''})



# do another loop to compute all possible moments
if moments:
    for _icnt in range (0,_ncnt):

        _sym_arr = _sym[_icnt]
        _s2 =int(_mult[_icnt])

        _ms = 0
        if ((_ms == 0) and (_s2  == 1)):
            _msc = 1
        elif ((_ms == 0) and (_s2 == 3)):
            _msc = -1
        else:
            quit_error('EOMEE knows only singlets and triplets, found: (2S+1) = '+str(_s2))

        for _isym in range (0,_nsym):

            n_root = int(_sym_arr[_isym])

            if( n_root == 0):
                continue

            _ms_0 = 0

            _extension = '_' + str(_isym+1) + '_' + str(_msc)

            for _jsym in range (0,_nsym):

                n_root2 = int(_sym_arr[_jsym])

                if ( n_root2 == 0):
                    continue

                _extension2 = '_' + str(_jsym+1) + '_' + str(_msc)

                new_target('EOMEE_TM'+_extension+_extension2)
                required()

                depend('EOMEE_DEN','EOMEE_LSOLVE'+_extension,'EOMEE_SOLVE'+_extension2)

                n_dens = n_root * n_root2
                ijsym = multd2h[_isym][_jsym]

                _extension12 = '_' + str(ijsym) + '_' + str(_msc)

                DEF_ME_LIST({LIST:'ME_CC_EOM_D'+_extension+_extension2,OPERATOR:'CC_EOM_D',
                    IRREP:ijsym,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})

                ASSIGN_ME2OP({LIST:'ME_LEX'+_extension,OPERATOR:'LEX'})
                ASSIGN_ME2OP({LIST:'ME_REX'+_extension2,OPERATOR:'REX'})

                OPTIMIZE({LABEL_OPT:'FOPT_EOM_D'+_extension+_extension2,LABELS_IN:'F_EOM_D'})

                if use_z:
                    depend('EOMEE_Z')
                    DEF_ME_LIST({LIST:'ME_CC_EOM_DZ'+_extension+_extension2,OPERATOR:'CC_EOM_DZ',
                        IRREP:ijsym,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                    DEF_ME_LIST({LIST:'ME_ZRHS'+_extension+_extension2,OPERATOR:'ZRHS',
                        IRREP:ijsym,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                    DEF_ME_LIST({LIST:'ME_Z'+_extension+_extension2,OPERATOR:'Z',
                        IRREP:ijsym,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                    DEF_ME_LIST({LIST:'ME_Z_A'+_extension+_extension2,OPERATOR:'Z_A',
                        IRREP:ijsym,'2MS':0,AB_SYM:+1,MIN_REC:1,MAX_REC:n_dens})
                    DEF_ME_LIST({LIST:'ME_DIAG_Z'+_extension+_extension2,OPERATOR:'DIAG',
                        IRREP:ijsym,'2MS':0,AB_SYM:+1})

                    OPTIMIZE({LABEL_OPT:'FOPT_EOM_Z'+_extension+_extension2,LABELS_IN:['F_EOM_ZRHS','F_Z_A']})
                    OPTIMIZE({LABEL_OPT:'FOPT_EOM_DZ'+_extension+_extension2,LABELS_IN:'F_EOM_DZ'})

                    PRECONDITIONER({LIST_PRC:'ME_DIAG_Z'+_extension+_extension2,
                        LIST_INP:'H0'})

                    PRINT({STRING:'Solving Z equations for symmetry pair '+str(_isym+1)+'/'+str(_jsym+1)})
                    SOLVE_LEQ({LIST_OPT:'ME_Z'+_extension+_extension2,
                        LIST_PRC:'ME_DIAG_Z'+_extension+_extension2,
                        OP_MVP:'Z_A',
                        OP_SVP:'Z',
                        OP_RHS:'ZRHS',
                        FORM:'FOPT_EOM_Z'+_extension+_extension2,
                        MODE:'DIA',
                        N_ROOTS:n_dens})


                for istate in range(n_root):
                    SET_STATE({LISTS:'ME_LEX'+_extension,ISTATE:istate+1})

                    for jstate in range(n_root2):
                        SET_STATE({LISTS:'ME_REX'+_extension2,ISTATE:jstate+1})

                        idens = istate*n_root2+jstate
                        SET_STATE({LISTS:'ME_CC_EOM_D'+_extension+_extension2,ISTATE:idens+1})
                        if use_z:
                            SET_STATE({LISTS:'ME_CC_EOM_DZ'+_extension+_extension2,ISTATE:idens+1})
                            SET_STATE({LISTS:'ME_Z'+_extension+_extension2,ISTATE:idens+1})
                        
                        EVALUATE({FORM:'FOPT_EOM_D'+_extension+_extension2})
                        if use_z:
                            EVALUATE({FORM:'FOPT_EOM_DZ'+_extension+_extension2})

                        PRINT({STRING:''})
                        
                        if _isym==_jsym and istate==jstate:
                            PRINT({STRING:'Excited state properties for state '+str(_isym+1)+'.'+str(istate+1)})
                        else:
                            PRINT({STRING:'Transition moments for states '+str(_isym+1)+'.'+str(istate+1)+'/'+str(_jsym+1)+'.'+str(jstate+1)})
                        EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_D'+_extension+_extension2})
                        if (use_z):
                            PRINT({STRING:'Including Z contribution:'})
                            EVAL_PROP({RANK:1,DENS:'ME_CC_EOM_DZ'+_extension+_extension2})

                PRINT({STRING:''})


export_targets();
