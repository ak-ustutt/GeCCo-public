# A dummy target file, to make the preliminary checks and print warnings
#
# yuri, nov 2014
import sys,os
import subprocess
sys.path=[os.getenv("GECCO_DIR")]+sys.path
from python_interface.gecco_interface import *
try:
    gitproc = subprocess.Popen(['git','rev-parse','HEAD'],
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=os.getenv("GECCO_DIR") )
except OSError as ex:
    print "Could not get current version information"
else:
    err = gitproc.stderr.read()
    if (err==""):
        print "Current git revision: "+gitproc.stdout.read(),
    else:
        print "Could not get current version information"
        
print "Using python version: ", sys.version
inp = GeCCo_Input()
orb = Orb_Info()

export_targets();
