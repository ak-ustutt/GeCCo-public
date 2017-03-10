#single entry file for all methods implemented by Arne Bargholz

import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *


print( "python start"+"-"*50)

new_target("do all",True)

print keywords.data

if ( keywords.is_keyword_set("method.MR_P") ) :
    print "setting MR python code"
    import python_blocks.MR_P

if ( keywords.is_keyword_set("method.MRCC2") ):
    print "begin setting MRCC2 targets"
    import python_blocks.MRCC2

if ( keywords.is_keyword_set("method.MRCCPT2") ):
    print "begin setting MRCCPT2 targets"
    import python_blocks.MRCCPT2

if ( keywords.is_keyword_set("method.unit_test") ):
    print "begin setting unit-test targets"
    import python_blocks.unit_test
    
if ( keywords.is_keyword_set("method.MRCC2.excite")) : 
    print "begin setting MRCC2 response targets"
    import python_blocks.response

if ( keywords.is_keyword_set("method.R12.SC")) :
    print "begin setting SC targets"
    import python_blocks.singles_correction


export_targets();

print ("python target_setting ends"+"-"*50)
