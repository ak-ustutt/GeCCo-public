from python_interface.gecco_interface import *
import sys
import unittest as ut
import python_interface.gecco_modules.stf_pack.stf_tests as test


suite=ut.TestLoader().loadTestsFromTestCase( test.Test_OPProduct )
suite.addTests(
    ut.TestLoader().loadTestsFromTestCase(
        test.Test_Bracket
    )
)
ut.TextTestRunner(stream=sys.stdout, verbosity=1).run(suite)
