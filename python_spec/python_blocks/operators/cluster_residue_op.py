
from python_interface.gecco_interface import *
from python_interface.gecco_modules.NoticeUtil import *

from python_interface.gecco_modules.omg_generator import *

i_am='cluster_residue'

orbitals=Orb_Info()
keywords=GeCCo_Input()

spinadapt=0
if keywords.is_keyword_set('calculate.routes.spinadapt'):
    spinadapt=int(keywords.get('calculate.routes.spinadapt'))

shift_terms = False
if keywords.is_keyword_set('method.MRCC2.shift_terms'):
    shift_terms=True if keywords.get('method.MRCC2.shift_terms') == 'T' else False


#-----------------------------------------------------------------#
# operators associated with T
#-----------------------------------------------------------------#
minexc=1
maxexc=2
if keywords.is_keyword_set('method.MR.minexc'):
  minexc=int(keywords.get('method.MR.minexc'))
if keywords.is_keyword_set('method.MR.maxexc'):
  maxexc=int(keywords.get('method.MR.maxexc'))
if (minexc==1 and maxexc==2):
  t1_shape='V,H|P,V|P,H'
  t2g_shape='V,H|VV,VH|VV,HH|P,V|PV,VV|P,H|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH' #compatible with Matthias
  t2ps_shape='V,H|VV,VH|P,V|PV,VV|P,H|PV,HV' #pseudo-doubles
  to0_shape=t1_shape
  to1_shape=t2g_shape
  to0_frm='To0=T1'
  to1_frm='To1=T2g'
  lamo0_frm='LAMo0=LAM1'
  lamo1_frm='LAMo1=LAM2g'
  useT1=True
  if (shift_terms):
      to0_shape='V,H|VV,VH|P,V|PV,VV|P,H|PV,HV'
      to1_shape='VV,HH|PV,HH|PP,VV|PP,HV|PP,HH'
      to0_frm='To0=T1+T2g'
      to1_frm='To1=T2g'
      lamo0_frm='LAMo0=LAM1+LAM2g'
      lamo1_frm='LAMo1=LAM2g'
elif (minexc==2 and maxexc==2):
  t2g_shape='VV,VH|VV,HH|PV,VV|PV,HV|PV,HH|PP,VV|PP,HV|PP,HH'
  t1_shape=',' # just to not leave it undefined
  t2ps_shape=',' # just to not leave it undefined
  useT1=False
  to0_shape=t1_shape
  to1_shape=t2g_shape
  # ... not clear that this works ... so ...
  raise Exception(i_am+": minexc==maxexc==2 -> check this case!")
  to0_frm='To0=T1'
  to1_frm='To1=T2g'
  lamo0_frm='LAMo0=LAM1'
  lamo1_frm='LAMo1=LAM2g'
else:
  raise Exception(i_am+": covered only minexc=1,2 and maxexc=2 so far ...")



#-----------------------------------------------------------------#
# targets to be included from other places
#-----------------------------------------------------------------#

new_target('DEF_T')
print("passing DEF_T!!")
depend('DEF_T2g')
if (useT1):
   depend('DEF_T1')


new_target('DEF_O')
depend('DEF_O2g')
if (useT1):
   depend('DEF_O1')


new_target('DEF_LAM')
depend('DEF_LAM2g')
if (useT1):
   depend('DEF_LAM1')



#-----------------------------------------------------------------#
# specific targets (to be included with care)
#-----------------------------------------------------------------#

# formal operators for perturbation order exp.
# + formulae for replacement
new_target('DEF_ToX')
depend('DEF_T')
depend('DEF_LAM')
depend('DEF_O')
PRINT({STRING:''})
PRINT({STRING:'Info on pertubation definition:'})
PRINT({STRING:'Zeroth order: '+to0_shape})
PRINT({STRING:'First order : '+to1_shape})

DEF_OP_FROM_OCC({
        LABEL:'To0',
        DESCR:to0_shape})

# pseudo-doubles
DEF_OP_FROM_OCC({
        LABEL:'T2ps',
        DESCR:t2ps_shape})

DEF_OP_FROM_OCC({
        LABEL:'To1',
        DESCR:to1_shape})

DEF_FORMULA({LABEL:'FORM_To0',FORMULA:to0_frm})
DEF_FORMULA({LABEL:'FORM_To1',FORMULA:to1_frm})

CLONE_OPERATOR({
        LABEL:'LAMo0',
        TEMPLATE:'To0',
        ADJOINT:True})

CLONE_OPERATOR({
        LABEL:'LAM2ps',
        TEMPLATE:'T2ps',
        ADJOINT:True})

CLONE_OPERATOR({
        LABEL:'LAMo1',
        TEMPLATE:'To1',
        ADJOINT:True})

DEF_FORMULA({LABEL:'FORM_LAMo0',FORMULA:lamo0_frm})
DEF_FORMULA({LABEL:'FORM_LAMo1',FORMULA:lamo1_frm})



new_target('DEF_T2g')
comment('Cluster Operators')


DEF_OP_FROM_OCC({
        LABEL:'T2g',
        DESCR:t2g_shape})

ME_param={
        LIST:'ME_T2g',
        OPERATOR:'T2g  ',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_T2g',info_only=True)





new_target('DEF_T1 ')

DEF_OP_FROM_OCC({
        LABEL:'T1',
        DESCR:t1_shape})

ME_param={
        LIST:'ME_T1',
        OPERATOR:'T1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_T1',info_only=True)




new_target('DEF_O2g')
comment('Residue Operators')

DEF_OP_FROM_OCC({
        LABEL:'O2g',
        JOIN:2,
        DESCR:omg_generator(t2g_shape)})

ME_param={
        LIST:'ME_O2g',
        OPERATOR:'O2g',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)



debug_MEL('ME_O2g',info_only=True)





new_target('DEF_O1')

DEF_OP_FROM_OCC({
        LABEL:'O1',
        JOIN:2,
        DESCR:omg_generator(t1_shape)})

ME_param={
        LIST:'ME_O1',
        OPERATOR:'O1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('ME_O1',info_only=True)




new_target('DEF_LAM2g')
comment('Lambda Operators')
depend('DEF_T2g')


CLONE_OPERATOR({
        LABEL:'LAM2g',
        TEMPLATE:'T2g  ',
        ADJOINT:True})

ME_param={
        LIST:'ME_LAM2g',
        OPERATOR:'LAM2g',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0


DEF_ME_LIST(ME_param)



debug_MEL('ME_LAM2g',info_only=True)





new_target('DEF_LAM1')
comment('Lambda Operators')
depend('DEF_T1')

CLONE_OPERATOR({
        LABEL:'LAM1',
        TEMPLATE:'T1',
        ADJOINT:True})

ME_param={
        LIST:'ME_LAM1',
        OPERATOR:'LAM1',
        IRREP:1,
        '2MS':0,
        AB_SYM:+1}
if spinadapt>=2:
    ME_param['S2']= 0

DEF_ME_LIST(ME_param)

debug_MEL('ME_LAM1',info_only=True)
