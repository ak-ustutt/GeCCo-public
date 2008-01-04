#-----------------------------------------------------------------------------#
#  read through the file chk_name and store the commands ....
#
#-----------------------------------------------------------------------------#
def parse_chk(chk_name):

    import sys

    key_list = ('title','molecule','input','check','end')
    subkeys = {}
    subkeys['check'] = ('label','find','rewind','skip','compare')

    # open the chk file
    try:
        chk_file = open(chk_name, 'r')
    except:
        print "Error: File " + chk_name + " was not found" ; sys.exit(2)

    commands = {}  # create empty dict

    check_for_subkey = False
    while 1:
        # read a new line from file
        line  = chk_file.readline()
        if line == '': break
        items = line.split(None,1)
        if len(items) == 0: continue
        
        # skip comment lines
        if items[0] == '#' or items[0] == '!':
            continue
        
        # a keyword is always the first item
        key = items[0].lower()

        # not present in the list of main keywords?
        if not key in key_list:
            # if we were not checking for a subkeyword this must be an error
            if not check_for_subkey:
                print "Error: Unknown keyword -- " + items[0]
                print "  list of known keywords: "
                print "    " + str(key_list)
                sys.exit(2)
            if not key in subkeys[current_key]:
                print "Error: Unknown keyword -- " + items[0]
                print "  list of known keywords: "
                print "    " + str(key_list)
                print "  allowed subkeywords for " + current_key
                print "    " + str(subkeys[current_key])
                sys.exit(2)
        else :
            check_for_subkey = False

        # code to be executed for new subkeyword entry
        if check_for_subkey:

            # append to list referenced by latest occurrence of the
            # main keyword
            #commands[current_key][-1].append(key)
            # append a list consisting of the parameters
            try:
                commands[current_key][-1].append([key,[items[1].rstrip()]])
            except:
                # empty list, if no parameters present
                commands[current_key][-1].append([key,[]])

        # code to be executed for new main keyword entry
        if not check_for_subkey:

            if not key in commands:
                # create new list: one entry for each time the command
                # has appeared in the input file
                commands[key] = []

            # is this a keyword with or without subkey?
            check_for_subkey = key in subkeys
            if check_for_subkey:
                # remember main keyword
                current_key = key
                ## add dictionary
                #commands[key].append({})
                # add list
                commands[key].append([])
            
        if not check_for_subkey:
            # keywords without subkeys
            # store further items in this line or next line
            commands[key].append([])
            try:
                commands[key][-1].append(chk_file.readline().lstrip().rstrip())
            except:
                print "Error: confused while reading arguments to " + key

    chk_file.close()
    
    return commands
