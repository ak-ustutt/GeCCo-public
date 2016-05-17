import re
from gecco_interface import *


tests=keywords.get('method.unit_test.tests')
if re.search("MODIFY_MEL", tests) is not None:
    import tests.mel

if re.search("string_to_form", tests) is not None:
    import tests.stf






