# A dummy target file, to make the preliminary checks and print warnings
#
# yuri, nov 2014
import sys,os
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *

print "Using python version: ", sys.version
inp = GeCCo_Input()
orb = Orb_Info()

export_targets();
