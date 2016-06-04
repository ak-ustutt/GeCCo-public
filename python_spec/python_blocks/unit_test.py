import re
from gecco_interface import *


testing=keywords.get('method.unit_test.tests')

if re.search("MODIFY_BLOCK", testing) is not None:
    import tests.mel

if re.search("string_to_form", testing) is not None:
    import tests.stf






