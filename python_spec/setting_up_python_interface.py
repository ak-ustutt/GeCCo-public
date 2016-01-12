# A dummy target file, to make the preliminary checks and print warnings
#
# yuri, nov 2014
try:
    from gecco_interface import *
except ImportError as ex:
    print("warning gecco_interface could not be imported from $PYTHONPATH" )
    print("The following error message was recieved") 
    print(ex.message)
import sys,os
sys.path=sys.path+[os.getenv("GECCO_DIR")+"/python_interface"]

from gecco_interface import *
inp = GeCCo_Input()
orb = Orb_Info()

export_targets();
