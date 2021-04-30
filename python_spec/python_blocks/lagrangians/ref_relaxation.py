""" Formulas for reference relaxation

History:
Yuri Aoto, August 2017 - Creation, based on Arne's code

"""
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

def _me_list_label(root, sym, spin, ms, msc, spinadapt):
    if spinadapt :
        return "{root}G{sym:1}S{spin:0>2:}M{ms:0>2}".format(root = root,
                                                            sym=sym,
                                                            spin=spin,
                                                            ms=ms,
        )
    else:
        return {0:"{root}G{sym:1}SxxM{ms:0>2}",
                1:"{root}G{sym:1}C+1M{ms:0>2}",
                -1:"{root}G{sym:1}C-1MM{ms:0>2}",}[msc].format(root = root,
                                                               sym=sym,
                                                               ms=ms,)
def _ab_sym(ms, mult):
    if (ms == 0) and ((mult-1)%4 == 0) :
        return 1
    elif (ms == 0) and ((mult+1)%4 == 0) :
        return -1
    else :
        return 0


def make_form_for_optref_minus3(E_form, E_target):
    """Creates the formula necessary for optref=-3

    Creates also some necessary ME-list, as required in the 'quick & dirty'
    implementation
    
    E_form (str)     Formula for the Energy
    E_target (str)   Target where the Energy is defined
    """
    new_target("FOPT_OMG_C0")
    depend(E_target)
    depend("DEF_A_C0")
    depend("MakeRefState")
    
    DERIVATIVE({
            LABEL_IN:E_form,
            LABEL_RES:'FORM_A_C0',
            OP_RES:'A_C0',
            OP_DERIV:'C0^+'})
    
    debug_FORM('FORM_A_C0')
    
    OPTIMIZE({
            LABEL_OPT:'FOPT_OMG_C0',
            LABELS_IN:'FORM_A_C0'})

    DEF_ME_LIST({LIST:_me_list_label("DIA",orbitals.get('lsym'),0,0,0,False)+"C0",
                 OPERATOR:"D0",
                 IRREP:int(orbitals.get('lsym')),
                 "2MS":int(orbitals.get('ims')),
                 AB_SYM:_ab_sym(int(orbitals.get('ims')),int(orbitals.get('imult')))
                 })

    SCALE_COPY({LIST_RES:_me_list_label("DIA",orbitals.get('lsym'),0,0,0,False)+"C0",
                LIST_INP:"ME_D0",
                FAC:1,})
