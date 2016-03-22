#-----------------------------------------------------------------------------#
#  perform the test
#-----------------------------------------------------------------------------#
def runtest_kernel(commands,options,testname):

    import sys,os,shutil

    RUN_OK = 0
    RUN_FAILED = 1
    USAGE_ERROR = 2
    SETUP_ERROR = 3

    print 'Next test:'
    for line in commands['title'][-1]:
        print line
    print 'Output in: ' + testname

    if options.program.lstrip() == ' ':
        print 'the program to test was obviously not specified ...'
        return USAGE_ERROR
        
    curdir = os.getcwd()
    scrdir = options.scratch
    scrdir_save = scrdir   # make VERY sure that we remove this and
                           # only this directory in the end

    # in particular, let's make sure that this is a new directory
    if os.path.isdir(scrdir):
        print "Error: directory does already exist!"
        print str(scrdir)
        return SETUP_ERROR

    out_name = os.path.join(curdir,testname)

    
    try:
        os.mkdir(scrdir)
    except:
        print "Error: could not create scratch directory!"
        print str(scrdir)
        return SETUP_ERROR

    moldir = os.path.join(options.mol_dir,commands['molecule'][-1][0].lstrip())
    mol_files = os.listdir(moldir)
    os.chdir(moldir)
    for mol_file in mol_files:
        try:
            shutil.copy(mol_file,scrdir)
        except:
            print "Error: could not copy file " + str(mol_file)
            return SETUP_ERROR
    os.chdir(curdir)

    input_name = os.path.join(options.inp_dir,commands['input'][-1][0].lstrip())
    try:
        shutil.copy(input_name,scrdir)
    except:
        print "Error: could not copy file " + str(input_name)
        return SETUP_ERROR

    os.chdir(scrdir)

    cmd = options.program + ' ' + commands['input'][-1][0].split()[0] + ' > ' + out_name
   # QUICK FIX for current MAC OS 10.11.3, should not affect other OSs
    ret_code = os.system('export DYLD_LIBRARY_PATH=$LIBRARY_PATH;' + cmd)

    os.chdir(curdir)

    print "now removing: " + scrdir
    if scrdir != scrdir_save:
        print 'CAVEAT: scrdir has changed during execution of rt_runtest_kernel.py!!'
        print '        scrdir      = ' + scrdir
        print '        scrdir_save = ' + scrdir_save
        print '        I will NOT remove any of the above directories!'
        return SETUP_ERROR
        
    shutil.rmtree(scrdir)

    if (ret_code == 0):
        return RUN_OK
    else:
        return RUN_FAILED
