"""single entry file for all methods implemented in Python

Hystory:
Arne Bargholz, January 2016 - Creation
"""

import sys,os
import traceback

sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *

print ("-"*50+"\npython_start")

new_target("do_all",True)

try:

    if ( keywords.is_keyword_set("method.MR_P") ) :
        print "setting MR python code"
        import python_blocks.MR_P

    if ( keywords.is_keyword_set("method.MRCC2") ):
        print "begin setting MRCC2 targets"
        import python_blocks.MRCC2

    if ( keywords.is_keyword_set("method.MRCC_new") ):
        print "begin setting icMRCC targets"
        import python_blocks.icMRCC

    if ( keywords.is_keyword_set("method.MRCCPT2") ):
        print "begin setting MRCCPT2 targets"
        import python_blocks.MRCCPT2

    if ( keywords.is_keyword_set("method.unit_test") ):
        print "begin setting unit-test targets"
        import python_blocks.unit_test

    if ( keywords.is_keyword_set("method.MRCC2.excite")):
        print "begin setting MRCC2 response targets"
        import python_blocks.response

    if ( keywords.is_keyword_set("method.R12.SC")) :
        print "begin setting SC targets"
        import python_blocks.singles_correction

    if ( keywords.is_keyword_set("method.MRCC_new") ):
        print "begin setting icMRCC targets"
        import python_blocks.icMRCC

    if ( keywords.is_keyword_set("method.MRCCPT2") ):
        print "begin setting MRCCPT2 targets"
        import python_blocks.MRCCPT2

    if ( keywords.is_keyword_set("method.ITF") ):
        print "begin setting icMRCC targets"
        import python_blocks.ITF

    print ("end of python_start\n"+"-"*50)

except Exception as ex:

    exc_type, exc_value, exc_traceback = sys.exc_info()
    print
    print 'Exception caught in python_start:'
    print "  Value: " + str(exc_value)
    print "  Type:  " + str(exc_type)
    print "  Traceback:"
    for tb in traceback.extract_tb(exc_traceback):
        print "    in {0}, line {1:d}: {2}, {3}".format(tb[0], tb[1], str(tb[2]), tb[3])
    print
    flog.close()

else:
    export_targets()

