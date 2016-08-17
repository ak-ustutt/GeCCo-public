#!/usr/bin/env python
# library imports:
import sys, math, re, os
from optparse import OptionParser

# import of own functions:
from rt_parse_chk import parse_chk
from rt_runtest_kernel import runtest_kernel
from rt_check import check

OK = 0
FAILED = 1
USAGE_ERROR = 2
SETUP_ERROR = 3

DEBUG = False

# defaults from the environment
user = os.environ.get('USER','user_not_set')
scrdir = os.environ.get('GECCO_TMP','')
if scrdir == '':
    # worth a try
    scrdir = os.path.join('/work',user)

# define the options
parse = OptionParser()
parse.add_option('-r','--run',dest='run',action='store_true',
                 default=False,
                 help='run the test')
parse.add_option('-c','--check',dest='check',action='store_true',
                 default=False,
                 help='check vs. reference file')
parse.add_option('-p','--program',dest='program',
                 default='not defined',
                 help='program to be tested (must be given if --run is active)')
parse.add_option('-s','--scratch',dest='scratch',
                 default=scrdir,
                 help='scratch directory')
parse.add_option('-m','--mol_dir',dest='mol_dir',default='molecules',
                 help='directory with molecule input data')
parse.add_option('-i','--inp_dir',dest='inp_dir',default='inputs',
                 help='directory with GeCCo input files')
# read command line
options, args = parse.parse_args(sys.argv[1:])

# some checks
try:
    basename = args[0]
except:
    print "ERROR: presumably you forgot to specify the test basename"
    parse.print_help()
    sys.exit(1)
if options.run and options.program == 'not defined':
    print "ERROR: the program to be tested was not given and must be given"
    parse.print_help()
    sys.exit(1)

# assign names
testname = "%s.out" % (basename)
chk_name = "%s.chk" % (basename)
err_name = "%s.err" % (basename)
chk_commands = {}
chk_commands = parse_chk(chk_name)

if DEBUG:
    print "resulting chk_commands structure:"
    print chk_commands

runtoken_dict={OK:"run_OK",
               FAILED:"run_failed",
               USAGE_ERROR:"setup_buggy",
               SETUP_ERROR:"setup_buggy",
}
if options.run:
    code = runtest_kernel(chk_commands,options,testname,basename)
    with open(err_name,"w") as f:
        f.write(runtoken_dict.get(code, "Internal ERROR: undefined return code"+str(code)))


checktoken_dict={OK:"check_OK",
                FAILED:"check_failed",
                USAGE_ERROR:"check_buggy",
                SETUP_ERROR:"check_buggy",
}
if options.check:
    code = check(chk_commands,options,basename)
    with open(err_name,"w") as f:
        f.write(checktoken_dict.get(code, "Internal ERROR: undefined return code"+str(code)))
if code != OK:
	sys.exit(1)
sys.exit(0)
