import re
from python_interface.gecco_interface import *


testing=keywords.get('method.unit_test.tests')

if re.search("SCALE_COPY", testing) is not None:
    import tests.scale_copy

if re.search("MODIFY_BLOCK", testing) is not None:
    import tests.mel

if re.search("ADD_UNITY", testing) is not None:
    import tests.add_unity

if re.search("string_to_form", testing) is not None:
    import tests.stf






