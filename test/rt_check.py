def check(chk_commands,options,basename):

    import re,sys,math

    CHECK_OK = 0
    CHECK_FAILED = 1
    SETUP_ERROR = 2
    USAGE_ERROR = 3

    DEBUG=False

    # file names
    testname = "%s.out" % (basename)
    ref_name = "%s.ref" % (basename)

    res_name = "%s.res" % (basename)

    # try to open the three files
    try:
        testfile = open(testname, 'r')
    except:
        print("Error: File " + testname + " was not found") ; return USAGE_ERROR
        
    try:
        ref_file = open(ref_name, 'r')
    except:
        print("Error: File " + ref_name + " was not found") ; return SETUP_ERROR

    try:
        res_file = open(res_name, 'w')
    except:
        print("Error: File " + res_name + " could not be opened") ; return SETUP_ERROR

    res_file.write(chk_commands.get('title',[['*** no title ***']])[-1][-1] + '\n')

    # overall result flag
    all_OK = True

    # loop through check commands
    for check in chk_commands['check']:

        if DEBUG:
            print('top of outer loop')
            print('check = ' + str(check))

        check_OK = True
        line_OK  = True

        args = []
        # loop over sub commands of check block
        for subkey, args in check:

            if DEBUG:
                print('current subkey: ' + subkey)
                print('current args:   ' + str(args))
            
            if subkey == 'rewind' :
                # rewind the two files
                testfile.seek(0)
                ref_file.seek(0)
            elif subkey == 'find' :
                # find the next occurrence of the given REGEX
                pattern = args[0]
                if DEBUG:
                    print('Looking for regex: ' + pattern)

                success1 = False
                success2 = False
                while 1:
                    ref_line = ref_file.readline()
                    if ref_line == '': break
                    if re.search(pattern,ref_line):
                        success1 = True; break
                while 1:
                    testline = testfile.readline()
                    if testline == '': break
                    if re.search(pattern,testline):
                        success2 = True; break

                if DEBUG:
                    print('success: ' + str(success1) + str(success2))
                    print('line in ref:  ' + str(ref_line))
                    print('line in test: ' + str(testline))

                # no success for reference file? pathologic:
                if not success1:
                    print('ERROR: pattern not found in reference file')
                    print(' the pattern was: ' + pattern)
                    if success2:
                        print('  NOTE: the pattern was found in the test file')
                        return SETUP_ERROR
                elif not success2:
                    line_OK = False

            elif subkey == 'skip':
                # skip the given amount of (non-blank) lines
                items = args[0].split()
                skip = int(items[0])

                if DEBUG:
                    print('skipping ' + str(skip) + ' lines')

                count = 0
                while count < skip:
                    testline = testfile.readline()
                    count += 1
                count = 0
                while count < skip:
                    ref_line = ref_file.readline()
                    count += 1 

                if DEBUG:
                    print('line in ref:  ' + str(ref_line))
                    print('line in test: ' + str(testline))

            elif subkey == 'label' :

                label = args[0]

            elif subkey == 'compare' :

                if DEBUG:
                    print('COMPARE:')
                    print('line in ref:  ' + str(ref_line))
                    print('line in test: ' + str(testline))

                items = args[0].split()

                # did we find the line?
                if line_OK: 

                    try:
                        type = items[0]
                        col  = int(items[1])-1
                        if type == 'real' :
                            tol  = float(items[2])
                    except:
                        print('SYNTAX: compare <type> <column> [<tolerance>]')
                        return SETUP_ERROR

                    try:
                        if col < 0:
                            ref_col = ref_line
                        else :
                            ref_items = ref_line.split()
                            ref_col   = ref_items[col]
                    except:
                        print('ERROR: column '+ str(col+1) + ' not found in reference')
                        print('  line was: ' + str(ref_line))         
                        return SETUP_ERROR

                    try:
                        if col < 0:
                            testcol = ref_line
                        else :
                            testitems = testline.split()
                            testcol   = testitems[col]
                    except:
                        line_OK = False

                    if line_OK:
                        if type == 'real' :
                            try :
                                check_OK = (check_OK and
                                            math.fabs(float(ref_col)-
                                                      float(testcol)) <= tol)
                            except:
                                check_OK = False
                        elif type == 'int' :
                            try :
                                check_OK = (check_OK and
                                            int(ref_col) == int(testcol))
                            except:
                                check_OK = False
                        elif type == 'str' :
                            check_OK = check_OK and ref_col == testcol
                        
            else :
                print("Unknown keyword: " + subkey); return SETUP_ERROR

        if not line_OK :
            result = 'FAILED - LINE NOT FOUND'
        elif not check_OK :
            result = 'FAILED'
        else :
            result = 'OK'

        message = " *** %-40s -- %-30s" % (label,result)
        print(message)
        res_file.write(message + "\n")

        all_OK = all_OK and check_OK and line_OK
    

    ref_file.close()
    testfile.close()
    res_file.close()

    if all_OK:
        return CHECK_OK
    else:
        return CHECK_FAILED
    
