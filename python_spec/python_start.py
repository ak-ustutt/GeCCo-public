
import sys,os
sys.path=sys.path+[os.getenv("GECCO_DIR")+"/python_interface"]

from gecco_interface import *


print("python start"+"-"*50)

new_target("do all",True)

if ( keywords.is_keyword_set("method.MR_P") ) :
    print "setting MR python code"
    import python_blocks.MR_P

if ( keywords.is_keyword_set("method.MRCC2") ):
    print "begin setting MRCC2 targets"
    import python_blocks.MRCC2

if ( keywords.is_keyword_set("method.MRCCPT2") ):
    print "begin setting MRCCPT2 targets"
    import python_blocks.MRCCPT2



export_targets();

print "python target_setting ends"+"-"*50 
