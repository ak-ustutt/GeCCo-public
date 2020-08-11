from python_interface.gecco_interface import *
import python_interface.gecco_modules.string_to_form as stf


from python_interface.gecco_modules.NoticeUtil import *
_op_list={'R1_prime_q':'T1',
          'R1_q':'T1',
          'R2g_prime_q':'T2g',
          'R2g_q':'T2g',
          'R_mu':'C0',
          'AR1_rspns_q':'O1',
          'AR2g_rspns_q':'O2g',
          'AR_rspns_mu':'H_C0',
          'SR1_rspns_q':'O1',
          'SR2g_rspns_q':'O2g',
          'SR_rspns_mu':'C0'}

_inp = GeCCo_Input()
_orb = Orb_Info()

_s2_0 = _orb.get('imult')
_isym_0 = _orb.get('lsym')
_ms_0 = _orb.get('ims')
_nsym = _orb.get('nsym')

_multd2h = [[1,2,3,4,5,6,7,8],
            [2,1,4,3,6,5,8,7],
            [3,4,1,2,7,8,5,6],
            [4,3,2,1,8,7,6,5],
            [5,6,7,8,1,2,3,4],
            [6,5,8,7,2,1,4,3],
            [7,8,5,6,3,4,1,2],
            [8,7,6,5,4,3,2,1]]



if ((_ms_0 == 0) and ((_s2_0-1 % 4) == 0)):
    _msc_0 = 1
elif ((_ms_0 == 0) and ((_s2_0+1 % 4) == 0)):
    _msc_0 = -1
else:
    _msc_0 = 0


relaxref = False
if keywords.is_keyword_set('method.MRCC2.relaxref'):
    relaxref = True

_sym = _inp.get('calculate.excitation.sym')
_mult = _inp.get('calculate.excitation.mult')

if (_sym == None):
    quit_error('Input missing: spatial symmetry of the excited state has not been set in the input')
else:
    _ncnt = len(_sym)

if (_mult == None):
    quit_error('Input missing: spin multiplicity of the excited state has not been set in the input')


_choice= _inp.get('calculate.solve.eigen.guess')
if (_choice == None):
    _choice = 0

_method="MRCC2"

new_target('INPUT_INFO')

PRINT({STRING: 'Doing icMRCC response in the ' + str(_method) + ' framework' })

PRINT({STRING: 'Irrep, S2, Ms of the reference state: ' + str(_isym_0) + ', ' +  str(_s2_0) + ', ' + str(_msc_0)})

PRINT({STRING: 'Factor for spin-combination: ' + str(_msc_0)})

new_target("FORM_EXCITED_ENERGY")
depend('DEF_RESPONSE_OPs')




DEF_SCALAR({LABEL:"Exc_En"})
#EE = stf.GenForm("FORM_EXCITED_ENERGY","Exc_En")
#EE += "<(AR1_rspns_q')*R1_q^+*(AR1_rspns_q')>"
#EE += "<(AR2g_rspns_q')*R2g_q^+*(AR2g_rspns_q')>"
#EE += "<R_mu^+*(AR_rspns_mu')>"
#EE.set_rule()

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_ENERGY',NEW:True,OP_RES:'Exc_En',FAC:1.0,FIX_VTX:True,
                    OPERATORS:["AR1_rspns_q", "R1_q^+", "AR1_rspns_q"],
                    IDX_SV   :[1, 2, 1],
                    AVOID    :[1,3]})

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_ENERGY',NEW:False,OP_RES:'Exc_En',FAC:1.0,FIX_VTX:True,
                   OPERATORS:["AR2g_rspns_q", "R2g_q^+", "AR2g_rspns_q"],
                   IDX_SV   :[1, 2, 1],
                   AVOID    :[1,3]})

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_ENERGY',NEW:False,OP_RES:'Exc_En',FAC:1.0,FIX_VTX:True,
                   OPERATORS:["R_mu^+", "AR_rspns_mu"],
                   IDX_SV   :[1, 2],
                   CONNECT  :[1,2]})


DEF_SCALAR({LABEL:"Exc_Sr"})

#SR = stf.GenForm("FORM_EXCITED_OVERLAPP","Exc_Sr")
#SR += "<(SR1_rspns_q')*R1_q^+*(SR1_rspns_q')>"
#SR += "<(SR2g_rspns_q')*R2g_q^+*(SR2g_rspns_q')>"
#SR += "<R_mu^+*(SR_rspns_mu')>"
#SR.set_rule()

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_OVERLAPP',NEW:True,OP_RES:'Exc_Sr',FAC:1.0,FIX_VTX:True,
                   OPERATORS:["SR1_rspns_q", "R1_q^+", "SR1_rspns_q"],
                   IDX_SV   :[1, 2, 1],
                   AVOID    :[1,3]})

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_OVERLAPP',NEW:False,OP_RES:'Exc_Sr',FAC:1.0,FIX_VTX:True,
                   OPERATORS:["SR2g_rspns_q", "R2g_q^+", "SR2g_rspns_q"],
                   IDX_SV   :[1, 2, 1],
                   AVOID    :[1,3]})

EXPAND_OP_PRODUCT({LABEL:'FORM_EXCITED_OVERLAPP',NEW:False,OP_RES:'Exc_Sr',FAC:1.0,FIX_VTX:True,
                   OPERATORS:["R_mu^+", "SR_rspns_mu"],
                   IDX_SV   :[1, 2],
                   CONNECT  :[1,2]})


new_target("FORM_LR_TRANSFORM")
depend("DEF_RESPONSE_OPs")
depend("MakeRefState")
INVARIANT({LABEL_RES:'FORM_R2g_q',
           LABEL_IN:'FORM_T2_orth',
           OP_RES:'R2g_q',
           OPERATORS:'DUMMY'})

REPLACE({LABEL_RES:'FORM_R2g_q',
         LABEL_IN:'FORM_R2g_q',
         OP_LIST:['T2_orth','R2g_prime_q']})


INVARIANT({LABEL_RES:'FORM_R1_q',
           LABEL_IN:'FORM_T1_orth',
           OP_RES:'R1_q',
           OPERATORS:'DUMMY'})

REPLACE({LABEL_RES:'FORM_R1_q',
         LABEL_IN:'FORM_R1_q',
         OP_LIST:['T1_orth','R1_prime_q']})
_solve_eqn_arr=[]

_first_iter = True

for _icnt in range (0,_ncnt):

    _sym_arr = _sym[_icnt]
    _s2 =int(_mult[_icnt])

    _ms = 0
    if ((_s2 % 2) == 0):
        _ms = 1
    if ((_ms == 0) and ((_s2 % 4) == 1)):
        _msc = 1
    elif ((_ms == 0) and ((_s2 % 4) == 3)):
        _msc = -1
    else:
        _msc = 0
    for _isym in range (0,_nsym):

        _no_root = int(_sym_arr[_isym])

        if( _no_root == 0):
            continue

        _isym_r = _multd2h[_isym][_isym_0-1]

        if (_s2 == _s2_0):
            s2_r = 1
        elif (abs(_s2-_s2_0) == 2):
            s2_r = 3
        else:
            quit_error('cannot handle this S2 difference')

        _ms_r = _ms - _ms_0
        _msc_r = 0

        if ((_ms_r == 0) and (s2_r == 1)):
                _msc_r = 1
        elif ((_ms_r == 0) and (s2_r == 3)):
                _msc_r = -1

        _extnsn = str(_isym+1) + '_' + str(_msc+1)

        _list_rspns = 'LIST_RSPNS_' + _extnsn


### Defining the lists for all the operators involved

        new_target(_list_rspns)
        depend('DEF_RESPONSE_OPs','BUILD_PRECON',
               'MakeOrthBasis','H0','DEF_PRECON_C0','DEF_FORM_PT_LAG2')
#               ,'DEF_ME_T','DIA_C0','DEF_ME_E(MR)',depend('DEF_FORM_PT_LAG2')'DEF_ME_X_TRM')

        _op_list={'R1_q':[_isym_r,_ms_r,_msc_r],
                  'R2g_q':[_isym_r,_ms_r,_msc_r],
                  'R_mu':[_isym+1,_ms_0,_msc],
                  'R2g_prime_q':[_isym_r,_ms_r,0],
                  'R1_prime_q':[_isym_r,_ms_r,0]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
                        '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],MIN_REC:1,MAX_REC:_no_root})

#       DEF_ME_LIST({LIST:'ME_DIAG_t'+_extnsn,OPERATOR:'DAI_T',IRREP:_op_list[_op][0],
#                    '2MS':0,AB_SYM:_op_list[_op][1],MIN_REC:1,MAX_REC:_no_root})
#
#       DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
#                    '2MS':0,AB_SYM:_op_list[_op][1],MIN_REC:1,MAX_REC:_no_root})

        _op_list={'AR1_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'AR2g_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'AR_rspns_mu':[_isym+1,_ms_0,_msc],
                  'SR1_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'SR2g_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'SR_rspns_mu':[_isym+1,_ms_0,_msc],}
#                  'INT_PPr':[_isym_r,_ms_r,_msc_r]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
                         '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2],
                         MIN_REC:1,
                         MAX_REC:_no_root})

        _op_list={'PRECON1':['ME_DIAG_t1',_isym_r,_ms_r],
                  'PRECON':['ME_DIAG_t2g',_isym_r,_ms_r],
                  'PRECONC0':['ME_DIAG_c',_isym+1,_ms_0]}

        for _op in _op_list:
            DEF_ME_LIST({LIST:_op_list[_op][0]+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][1],
                        '2MS':_op_list[_op][2]})

        DEF_ME_LIST({LIST:'ME_MINEN'+_extnsn,
                     OPERATOR:'PT_LAG',
                     IRREP:1,
                     '2MS':0})

        ASSIGN_ME2OP({LIST:'PT_LAG_LST',OPERATOR:'PT_LAG'})

        _diag_cal_q = 'DIAG_CAL_q_' + _extnsn

### Setting up the preconditioner corresponding to the equations to solve 'R_q'

        new_target(_diag_cal_q)

        depend(_list_rspns)
        depend('MAKE_FOCK_REF','MAKE_A_TRF',"EVAL_E0")

        debug_FORM('FORM_FOCK_REF')

        debug_MEL("ME_FOCK_REF")

        PRECONDITIONER({LIST_PRC:'ME_DIAG_t1'+_extnsn,
                        LIST_INP:'ME_FOCK_REF'})


        ASSIGN_ME2OP({LIST:'ME_X_TRM',
                     OPERATOR:'X_TRM'})

        debug_MEL("ME_X_TRM")

        EVALUATE({
                FORM:'FOPT_A_TRF_FINAL'
                })

        debug_MEL('A_TRF_LST')

        EXTRACT_DIAG({
                LIST_RES:'ME_DIAG_t1'+_extnsn,
                LIST_IN:'A_TRF_LST',
                MODE:'extend'})

        debug_MEL("ME_DIAG_t1"+_extnsn)
        # preconditioner for T2g

        PRECONDITIONER({LIST_PRC:'ME_DIAG_t2g'+_extnsn,
                        LIST_INP:'ME_FOCK_REF'})

        debug_MEL("ME_DIAG_t2g"+_extnsn)

        EXTRACT_DIAG({LIST_RES:'ME_DIAG_t2g'+_extnsn,
                      LIST_IN:'A_TRF_LST',
                      MODE:'extend'})


        debug_MEL("ME_DIAG_t2g"+_extnsn)

        _diag_cal_mu = 'DIAG_CAL_mu_' + _extnsn

### Setting up the preconditioner corresponding to the equations to solve 'R_q'

        new_target(_diag_cal_mu)

        depend(_list_rspns)

        PRECONDITIONER({LIST_PRC:'ME_DIAG_c'+_extnsn,
                        LIST_INP:'H0',
                        MODE:'dia-H'})

        SCALE_COPY({LIST_RES:'ME_MINEN'+_extnsn,
                    LIST_INP:'ME_E0',
                    FAC:-1.0})

        EXTRACT_DIAG({LIST_RES:'ME_DIAG_c'+_extnsn,
                     LIST_IN:'ME_MINEN'+_extnsn,
                     MODE:'ext_act'})

        debug_MEL("ME_DIAG_c"+_extnsn)
        _rspns_opt = 'RSPNS_OPT_' + _extnsn


### Optimising all the Formula.

        new_target(_rspns_opt)

        depend('DEF_FORMS_METRIC',_list_rspns,'DEF_FORM_AR_RSPNS_q','DEF_FORM_AR_RSPNS_mu',
               'MakeRefState',"FORM_LR_TRANSFORM")

        OPTIMIZE({LABEL_OPT:'RSPNS_OPT'+_extnsn,
                  LABELS_IN:['FORM_AR1_RSPNS_q','FORM_AR2g_RSPNS_q','FORM_AR_RSPNS_mu',
                             'FORM_S1','FORM_S2g','FORM_SR_RSPNS_mu',
                             'FORM_R1_q','FORM_R2g_q']})
#                  INTERM:'F_PPrint'})


#### projecting out the elements of the ground state during the solution for the
#### states with same symmetry as the ground state

        if (_first_iter):
            _prj_form = 'PRJ_FORM_'
            new_target(_prj_form)
            depend("DEF_RESPONSE_OPs",'RefState-Operators')
            EXPAND_OP_PRODUCT({LABEL:'F_prj',
                               OP_RES:'R_mu',
                               OPERATORS:['R_mu','C0','C0^+','R_mu','R_mu'],
                               IDX_SV:[1,2,3,4,1],
                               AVOID:[1,4,3,5],
                               FAC:-1.0})

            debug_FORM("F_prj" )

            OPTIMIZE({LABEL_OPT:'FOPT_prj',
                      LABELS_IN:'F_prj'})


            _first_iter = False

#### Finally solving the eigen value equation to get the excitation energies

        _solve_eqn = 'SOLVE_EQN_' + _extnsn
        _solve_eqn_arr.append(_solve_eqn)

        new_target(_solve_eqn)
        if (relaxref):depend("MAKE_GAM0_HMRCC2")
        depend('INPUT_INFO')
        depend(_rspns_opt)
        depend(_diag_cal_q)
        depend(_diag_cal_mu)
        depend(_prj_form)
        depend('MAKE_MRCC2')
        _solve_evp_basis={}
        _solve_evp_basis[LIST_OPT]=['ME_R1_q'+_extnsn,'ME_R2g_q'+_extnsn,'ME_R_mu'+_extnsn]
        _solve_evp_basis[LIST_PRC]=['ME_DIAG_t1'+_extnsn,'ME_DIAG_t2g'+_extnsn,'ME_DIAG_c'+_extnsn]
        _solve_evp_basis[OP_MVP]=['AR1_rspns_q','AR2g_rspns_q','AR_rspns_mu']
        _solve_evp_basis[OP_SVP]=['SR1_rspns_q','SR2g_rspns_q','SR_rspns_mu']
        _solve_evp_basis[FORM]='RSPNS_OPT'+_extnsn
        _solve_evp_basis[LIST_SPC]=['ME_R1_prime_q'+_extnsn,'ME_X_TRM','ME_X_TRM_DAG',"ME_R2g_prime_q"+_extnsn]
        _solve_evp_basis[MODE]='TRF TR0 PRJ'
        _solve_evp_basis[FORM_SPC]='FOPT_prj'
        _solve_evp_basis[N_ROOTS]=_no_root
        _solve_evp_basis[CHOICE_OPT]=_choice
        _solve_evp_basis[SOLVER] = "NEW"
        PRINT({STRING: ''})
        PRINT({STRING: 'Calculating excitation to irrep:   ' + str(_isym+1) +
                       '  with spin multiplicity:   ' + str(_s2)})

        #PRINT({STRING: 'isym_r, msc_r:' + str(_isym_r) + ',  ' + str(_msc_r)})

        SOLVE_EVP(_solve_evp_basis)


        _get_mel_inf = 'GET_MEL_INF' + _extnsn

        new_target(_get_mel_inf,True)
#
        depend(_solve_eqn)
        depend("FORM_EXCITED_ENERGY")
        DEF_ME_LIST({LIST:'ME_Exc_En'+_extnsn,
                     OPERATOR:'Exc_En',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:0,
                     MIN_REC:1,MAX_REC:_no_root})
        DEF_ME_LIST({LIST:'ME_Exc_Sr'+_extnsn,
                     OPERATOR:'Exc_Sr',
                     IRREP:1,
                     '2MS':0,
                     AB_SYM:0,
                     MIN_REC:1,MAX_REC:_no_root})
##

        _op_list={'R1_q':[_isym_r,_ms_r,_msc_r],
                  'R2g_q':[_isym_r,_ms_r,_msc_r],
                  'R_mu':[_isym+1,_ms_0,_msc],
                  'AR1_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'AR2g_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'AR_rspns_mu':[_isym+1,_ms_0,_msc],
                  'SR1_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'SR2g_rspns_q':[_isym_r,_ms_r,_msc_r],
                  'SR_rspns_mu':[_isym+1,_ms_0,_msc],
                  }
        #                  'INT_PPr':[_isym_r,_ms_r,_msc_r]}

#        _op='AR1_rspns_q'
#        ASSIGN_ME2OP({LIST:'ME_'+_op+_extnsn,
#                      OPERATOR:_op })

#        DEF_ME_LIST({LIST:'ME_'+_op+_extnsn,OPERATOR:_op,IRREP:_op_list[_op][0],
#                     '2MS':_op_list[_op][1],AB_SYM:_op_list[_op][2]})

        for _op in _op_list:
            ASSIGN_ME2OP({LIST:'ME_'+_op+_extnsn,
                          OPERATOR:_op })


        OPTIMIZE({LABEL_OPT:'FOPT_SR'+_extnsn,
                  LABELS_IN:['FORM_AR1_RSPNS_q','FORM_AR2g_RSPNS_q','FORM_AR_RSPNS_mu',
                             'FORM_S1','FORM_S2g','FORM_SR_RSPNS_mu']})
        debug_FORM("FORM_EXCITED_OVERLAPP")
        debug_FORM("FORM_EXCITED_ENERGY")

        OPTIMIZE({LABEL_OPT:"FOPT_Exc_En"+_extnsn,
                  LABELS_IN:["FORM_EXCITED_OVERLAPP","FORM_EXCITED_ENERGY"]})



        for i in xrange(1,_no_root +1):
            SET_STATE({LISTS:['ME_AR1_rspns_q'+_extnsn,
                              'ME_AR2g_rspns_q'+_extnsn,
                              'ME_AR_rspns_mu'+_extnsn,
                              'ME_SR1_rspns_q'+_extnsn,
                              'ME_SR2g_rspns_q'+_extnsn,
                              'ME_SR_rspns_mu'+_extnsn,
                              'ME_R1_q'+_extnsn,
                              'ME_R2g_q'+_extnsn,
                              'ME_R_mu'+_extnsn] ,
                       ISTATE:i})
            EVALUATE({FORM:'FOPT_SR'+_extnsn})
            debug_MEL('ME_R1_q'+_extnsn )
            debug_MEL('ME_R2g_q'+_extnsn )
            debug_MEL('ME_R_mu'+_extnsn )
            debug_MEL('ME_AR_rspns_mu'+_extnsn )

            EVALUATE({FORM:'FOPT_Exc_En'+_extnsn})

            debug_MEL('ME_Exc_En'+_extnsn)
            debug_MEL('ME_Exc_Sr'+_extnsn)

            SCALE({LIST_RES:'ME_Exc_En'+_extnsn,LIST_INP:'ME_Exc_En'+_extnsn,
                   LIST_SCAL:'ME_Exc_Sr'+_extnsn,FAC:1.0,INV:True})

            PUSH_RESULT({LIST:'ME_Exc_En'+_extnsn,COMMENT:'MRCC_STATE_'+str(_isym+1)+'.'+str(i), FORMAT:"SCAL F20.14"})


        for i in xrange(1,_no_root +1):
            SET_STATE({LISTS:['ME_R2g_q'+_extnsn,
                              ],
                       ISTATE:i})
            SCALE_COPY({LIST_INP:'ME_R2g_q'+_extnsn,
                        LIST_RES:'ME_R2g_q'+_extnsn,
                        LIST_SHAPE:'ME_SR2g_rspns_q'+_extnsn,
                        FAC:1,
                        MODE:"scale"
            })


        ANALYZE_MEL({LISTS:['ME_R_mu'+_extnsn,'ME_R1_q'+_extnsn,'ME_R2g_q'+_extnsn],
                     LISTS_CV:['ME_SR_rspns_mu'+_extnsn,'ME_SR1_rspns_q'+_extnsn,'ME_SR2g_rspns_q'+_extnsn]})
