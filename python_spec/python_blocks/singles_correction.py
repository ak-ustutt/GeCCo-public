#
# Evaluate singles correction according to
# Kong and Valeev, JCP 133, 174126 (2010)
#
# Created by Andreas
# Yuri, Nov 2016 - commented and translated to python
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import*

gecco_input = GeCCo_Input()


print_level = 0
alt_res_en = gecco_input.get('method.R12.SC.en_from_res')

# Check a sequence of things needed for singles correction
# Maybe turn on the necessary stuffs automaticaly?
if ( not( gecco_input.is_keyword_set('calculate.skip_E'))):
    quit_error('Use "calculate skip_E" for singles correction.')

if ( gecco_input.get('method.MR.maxexc') != '1'):
    quit_error('Use "method MR maxexc=1" for singles correction.')

if ( gecco_input.get('calculate.solve.non_linear.optref') != '0'):
    quit_error('Use "calculate solve non_linear optref=0" for singles correction.')




# -----
# Necessary operators
#
new_target( 'CSC_OPS')
DEF_SCALAR({LABEL:'LCSC'}) # energy

DEF_OP_FROM_OCC({LABEL:'C1', # single excitations to the complete virtual space
                 DESCR:'P,H|X,H|P,V|X,V'})

DEF_OP_FROM_OCC({LABEL:'O1', # omega, the residual
                 JOIN: 2,
                 DESCR:',,[PX],H,,V,[PX],'})

CLONE_OPERATOR({LABEL:   'C1tr',
                TEMPLATE:'C1'})


# -----
# ME lists of above operators
#
new_target( 'CSC_LISTS')
depend( 'CSC_OPS')
for op in ['LCSC','C1','C1tr','O1']:
    DEF_ME_LIST({LIST:'ME_'+op,
                 OPERATOR:op,
                 IRREP:1,
                 '2MS':0,
                 AB_SYM:+1
                 })


# -----
# Residual, eq. 15 from Kong and Valeev
#
new_target( 'CSC_RES')
depend( 'Favg', 'CSC_OPS')

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:True,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','C1^+','Favg','C1','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4   ,5]
                   })

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:False,
                   FAC:-1.0,
                   OP_RES:'LCSC',
                   FIX_VTX:True,
                   OPERATORS:['C0^+','C1^+','C0^+','Favg','C0','C1','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4     ,5   ,6   ,7],
                   AVOID:[1,5,  1,4,  3,7,  4,7,  2,4,  2,5,  3,6,  4,6]
                   })

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:False,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','C1^+','Favg','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4],
                   LABEL_DESCR:['3,,X,[HV]'], # only f_i^{\alpha'} and f_u^{\alpha'} contribute
                   CONNECT:[2,3]
                   })

DERIVATIVE({LABEL_RES:'CSC_RES',
            LABEL_IN:'CSC_FORM',
            OP_RES:'O1',
            OP_DERIV:'C1^+'
            })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'CSC_RES'})



# -----
# Energy (the singles correction), eq. 16 from Kong and Valeev
#
new_target( 'CSC_EN')
depend( 'Favg', 'CSC_OPS')

EXPAND_OP_PRODUCT({LABEL:'CSC_EN',
                   NEW:True,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','Favg','C1','C0'],
                   IDX_SV:   [1     ,2     ,3   ,4],
                   LABEL_DESCR:['2,,[HV],X'], # only f_{\alpha'}^i and f_{\alpha'}^u contribute
                   CONNECT:[2,3]
                   })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'CSC_EN'})




# -----
# Alternative definition of energy and residuals:
# energy and residuals from CSC_FORM
# -----
new_target( 'CSC_FORM')
depend( 'Favg', 'CSC_OPS')

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:True,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','C1^+','Favg','C1','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4   ,5]
                   })

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:False,
                   FAC:-1.0,
                   OP_RES:'LCSC',
                   FIX_VTX:True,
                   OPERATORS:['C0^+','C1^+','C0^+','Favg','C0','C1','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4     ,5   ,6   ,7],
                   AVOID:[1,5,  1,4,  3,7,  4,7,  2,4,  2,5,  3,6,  4,6]
                   })

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:False,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','C1^+','Favg','C0'],
                   IDX_SV:   [1     ,2     ,3     ,4],
                   LABEL_DESCR:['3,,X,[HV]','3,,[HV],X'],
                   CONNECT:[2,3]
                   })

EXPAND_OP_PRODUCT({LABEL:'CSC_FORM',
                   NEW:False,
                   OP_RES:'LCSC',
                   OPERATORS:['C0^+','Favg','C1','C0'],
                   IDX_SV:   [1     ,2     ,3   ,4],
                   LABEL_DESCR:['2,,X,[HV]','2,,[HV],X'],
                   CONNECT:[2,3]
                   })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'CSC_FORM'})
    
# -----
# Alternative definition of energy and residuals:
# energy and residuals from CSC_FORM
# -----
new_target( 'CSC_RES_alt')
depend( 'CSC_FORM')

DERIVATIVE({LABEL_RES:'CSC_RES',
            LABEL_IN:'CSC_FORM',
            OP_RES:'O1',
            OP_DERIV:'C1^+'
            })

INVARIANT({LABEL_RES:'CSC_RES',
           LABEL_IN:'CSC_RES',
           OP_RES:'O1',
           OPERATORS:'O1',
           REORDER:True
           })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'CSC_RES'})

# -----
# Alternative definition of energy and residuals:
# energy and residuals from CSC_FORM
# -----
new_target( 'CSC_EN_alt')
depend( 'CSC_RES_alt')

FACTOR_OUT({LABEL_RES:'CSC_EN',
            LABEL_IN:'CSC_FORM',
            INTERM:'CSC_RES'
            })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'CSC_EN'})





# -----
# Opt. formula
#
new_target( 'CSC_OPT')
if (alt_res_en):
    alt_R_E = '_alt'
else:
    alt_R_E = ''

depend( 'CSC_EN' + alt_R_E, 'CSC_RES' + alt_R_E, 'CSC_LISTS', 'Favg-INT', 'DEF_ME_C0')



OPTIMIZE({LABEL_OPT:'CSC_OPT',
          LABELS_IN:['CSC_RES', 'CSC_EN']
          })



# -----
# Preconditioner
#
new_target( 'DIAG_SC')
depend( 'CSC_OPS','Favg-INT','DEF_ME_C0','FOPT_Atr')

CLONE_OPERATOR({LABEL:'DIA1',
                TEMPLATE:'C1'
                })

DEF_ME_LIST({LIST:'DIAG_SC',
             OPERATOR:'DIA1',
             IRREP:1,
             '2MS':0,
             AB_SYM:+1
             })

PRECONDITIONER({LIST_PRC:'DIAG_SC',
                LIST_INP:'Favg-INT'
                })

EVALUATE({FORM:'FOPT_Atr'})

EXTRACT_DIAG({LIST_RES:'DIAG_SC',
              LIST_IN:'ME_A',
              MODE:'extend'
              })

if (print_level > 0):
    PRINT_MEL({LIST:'DIAG_SC'})

#PRINT_MEL({LIST:'Favg-INT'})
#PRINT_MEL({LIST:'ME_Fvvtr'})
#PRECONDITIONER({LIST_PRC:'DIAG_SC',
#                LIST_INP:'ME_Fvvtr'})





# -----
# Lists and formulas for mode TRF
#
new_target( 'IC_TRAF')
depend( 'CSC_OPS', 'CSC_LISTS', 'C0', 'EVAL_E_REF')

DEF_OP_FROM_OCC({LABEL:'D1',
                 JOIN: 2,
                 DESCR:',V;V,'})

CLONE_OPERATOR({LABEL:'Xtr',
                TEMPLATE:'D1'})

CLONE_OPERATOR({LABEL:'D1tr',
                TEMPLATE:'D1'})

CLONE_OPERATOR({LABEL:'D1trdag',
                TEMPLATE:'D1tr'})

EXPAND_OP_PRODUCT({LABEL: 'F_D1',
                   OP_RES:'D1',
                   OPERATORS:['D1','C0^+','D1','D1','C0','D1'],
                   IDX_SV:   [1   ,2     ,1   ,1   ,3   ,1]
                   })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'F_D1'})

EXPAND_OP_PRODUCT({LABEL:'F_C1tr',
                   OP_RES:'C1',
                   OPERATORS:['C1','C1tr','C1'],
                   IDX_SV:   [1   ,2     ,1],
                   LABEL_DESCR:'2,,[PX],H'
                   })

EXPAND_OP_PRODUCT({LABEL:'F_C1tr',
                   OP_RES:'C1',
                   OPERATORS:['C1','Xtr','C1tr','Xtr','C1'],
                   IDX_SV:   [1   ,2    ,3     ,2    ,1],
                   AVOID: [2,4],
                   NEW: False
                   })

#EXPAND_OP_PRODUCT({LABEL:'F_C1tr',
#                   OP_RES:'C1',
#                   OPERATORS:['C1','D1','C1tr','D1','C1'],
#                   IDX_SV:   [1   ,2   ,3     ,2   ,1],
#                   AVOID: [2,4],
#                   NEW: False
#                   })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'F_C1tr'})

DEF_ME_LIST({LIST:'ME_D1',
             OPERATOR:'D1',
             IRREP:1,
             '2MS':0,
             AB_SYM:+1
             })

DEF_ME_LIST({LIST:'ME_Xtr',
             OPERATOR:'Xtr',
             IRREP:1,
             '2MS':0,
             AB_SYM:+1
             })

OPTIMIZE({LABEL_OPT:'FOPT_D1',
          LABELS_IN:'F_D1'
          })

EVALUATE({FORM:'FOPT_D1'})

INVERT({LIST_INV:'ME_D1',
        LIST:'ME_Xtr',
        MODE:'invsqrt'})

OPTIMIZE({LABEL_OPT:'FOPT_C1tr',
          LABELS_IN:'F_C1tr'
          })

# -----
# Lists and formulas for mode TRF (?) Seem not to be needed
#
new_target( 'FTRAF')
depend('CSC_OPS','Favg-INT','IC_TRAF')

# formula for transformed Fock in VV space
DEF_OP_FROM_OCC({LABEL:'Fvvtr',
                 DESCR:'H,H|P,P|V,V|X,X',
                 })

DEF_ME_LIST({LIST:'ME_Fvvtr',
             OPERATOR:'Fvvtr',
             IRREP:1,
             '2MS':0,
             AB_SYM:+1
             })

EXPAND_OP_PRODUCT({LABEL:'F_FTRAF',
                   NEW:True,
                   OP_RES:'Fvvtr',
                   OPERATORS:['Fvvtr','Favg','Fvvtr'],
                   IDX_SV:   [1      ,2     ,1],
                   LABEL_DESCR:'2,,[HPX],[HPX]'
                   })

EXPAND_OP_PRODUCT({LABEL:'F_FTRAF',
                   NEW:False,
                   OP_RES:'Fvvtr',
                   FAC:1.0,
                   OPERATORS:['Fvvtr','Xtr','Xtr','Favg','Xtr','Xtr','Fvvtr'],
                   IDX_SV:   [1      ,2    ,3    ,4     ,2    ,3    ,1],
                   AVOID:[2,5,2,6,3,5,3,6]
                   })

if (print_level > 0):
    PRINT_FORMULA({LABEL:'F_FTRAF'})

OPTIMIZE({LABEL_OPT:'FOPT_FTRAF',
          LABELS_IN:'F_FTRAF'})




# -----
# Solve the equations
#
new_target( 'Singles_Correction', True)
depend( 'DIAG_SC','CSC_OPT','EVAL_E_REF','IC_TRAF' )

SOLVE_NLEQ({LIST_OPT:'ME_C1',
            LIST_RESID:'ME_O1',
            LIST_PRC:'DIAG_SC',
            MODE:'TRF',
            LIST_E:'ME_LCSC',
            FORM:'CSC_OPT',
            LIST_SPC:['ME_C1tr','ME_Xtr','ME_Xtr'],
            FORM_SPC:'FOPT_C1tr'
            })

