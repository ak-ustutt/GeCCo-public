# A dummy target file, to make the preliminary checks and print warnings
#
# yuri, nov 2014
import sys,os
interface_path=os.path.join(os.getenv("GECCO_DIR"),"python_interface")
sys.path=[interface_path]+sys.path
print "using python version:",sys.version
from gecco_interface import *
print 
inp = GeCCo_Input()
orb = Orb_Info()

export_targets();
