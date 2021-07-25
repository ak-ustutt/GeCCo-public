class itf_line:
    # TODO: Finish class
    """
    A line of ITF algo code

    Member data
        line:   Complete line of ITF code
        parts:  Line split into seperate elements
        output: Name of output file
    """

    def __init__(self, line, output):
        """ Return new itf_line object"""
        self.line = line
        self.parts = line.split()
        self.output = output

    def rename_line(self):
        """Add generic index to tensors within a line"""
        # Need to replace residual and amplitude (R and T) tensor names
        # with their names defined in C++, ie. R:eecc
        names = ["R", "G", "T", "K", "f", "K4E", "INTkx"]
        for i in range(0, len(names)):
            self.rename_line_names(names[i])
        self.line = " ".join(self.parts)

    def rename_line_names(self, name):
        """Add generic index to tensor name"""
        if (name+"[" in self.line):
            for i in range(0, len(self.parts)):
                if (name+"[" in self.parts[i]):
                    self.parts[i] = self.parts[i].replace(name, name + ":" + "".join(generic_index(self.parts[i])))


def print_inter(prev_lines, init=False):
    # Load, contract, drop tensors involved with intermediates
    # Loop over previously stored lines
    global tab

    for i in range(0, len(prev_lines)):
        words = prev_lines[i].split()
        print_loop(prev_lines[i], words)
        print_result(prev_lines[i], tab, init)
    tab = False


def change_line_names(name, line, words):
    # Add generic index to tensor name

    if (name+"[" in line):
        for i in range(0, len(words)):
            if (name+"[" in words[i]):
                words[i] = words[i].replace(name, name + ":" + "".join(generic_index(words[i])))


def change_line(line_o):
    # Add generic index to tensors within a line, return a new string

    # Need to replace residual and amplitude (R and T) tensor names
    # with there names defined in C++, ie. R:eecc
    words=line_o.split()
    change_line_names("R", line_o, words)
    change_line_names("G", line_o, words)
    change_line_names("T", line_o, words)
    change_line_names("K", line_o, words)
    change_line_names("J", line_o, words)
    change_line_names("f", line_o, words)
    change_line_names("K4E", line_o, words)
    change_line_names("INTkx", line_o, words)
    line = " ".join(words)
    return line


def print_drop_tensors(load_tensors, indent=False):
    # Just reverse the load_ten line for drop_ten
    load_split = load_tensors.replace(',','').split()
    drop_tensors = "drop "
    for i in range(len(load_split)-1, 0, -1):
        drop_tensors = drop_tensors + load_split[i]
        if (i != 1):
            drop_tensors = drop_tensors + ", "

    if (indent):
        print("    "+drop_tensors.rstrip(), file=out)
    else:
        print(drop_tensors.rstrip(), file=out)


# TODO: merge this with above
def print_drop_init_tensors(load_tensors, indent=False):
    # Just reverse the load_ten line for drop_ten
    load_split = load_tensors.replace(',','').split()
    drop_tensors = "drop "
    for i in range(len(load_split)-1, 0, -1):
        drop_tensors = drop_tensors + load_split[i]
        if (i != 1):
            drop_tensors = drop_tensors + ", "

    if (indent):
        print("    "+drop_tensors.rstrip(), file=init_res_temp)
    else:
        print(drop_tensors.rstrip(), file=init_res_temp)



def print_result(line, indent=False, init=False):
    # Load, contract, drop tensors involved with result tensors

    words=line.split()

    if (indent):
        line = "    " + line

    # Load tensors, cases depend on how many tensors are on the right
    if len(words)==3:
        # Either a simple adding or copying case
        if "TIN" not in words[2]:
            load_ten="load " + words[2].split('*',1)[-1]

            print(load_ten.strip(), file=out)
            if (init): print(load_ten.strip(), file=init_res_temp)

        print(line.strip(), file=out)
        if (init): print(line.strip(), file=init_res_temp)

        # Drop tensors
        if "TIN" not in words[2]:
            print_drop_tensors(load_ten)
            if (init): print_drop_init_tensors(load_ten)


    else:
        # Get list of non-intermediate tesnsors in a line
        t = []
        for i in range (1,len(words)):
            if ('+' in words[i] or '-' in words[i] or 'TIN' in words[i]):
                continue
            else:
                t.append(words[i].split('*',1)[-1].replace('(','').replace(')',''))

        # The line contains a contraction between two intermedites, so no
        # need to print a load/drop line
        if (len(t)==0):
            print(line, file=out)
            return

        # Remove tensor which are the same (ie. don't need to load them twice)
        seen = []
        result = []
        for item in t:
            if item.split('[',1)[0] not in seen:
                seen.append(item.split('[',1)[0])
                result.append(item)

        # Check if the line is within a loop
        if (indent):
            load_ten="    load "
            # Check if (T - T^+) tensors need to be loaded inside the loop
            load_extra_loop_tensors(line, result)
        else:
            load_ten="load "

        # Construct load line
        for tensor in result[:-1]:
            load_ten = load_ten + tensor + ', '
        load_ten = load_ten + result[-1]

        # Print out load, line, and drop
        print(load_ten, file=out)
        print(line, file=out)
        print_drop_tensors(load_ten, indent)

        if (init):
            print(load_ten, file=init_res_temp)
            print(line, file=init_res_temp)
            print_drop_init_tensors(load_ten, indent)


def add_to_global(word,declare_ten,declare_ten_index,declare_ten_name):
    # Add tensor to global lists, so it can be declared at the start of the algo file

    if word.split('*',1)[-1] not in declare_ten and "TIN" not in word:
        index=list(word[word.find("[")+1:word.find("]")])

        # Construct generic index
        generic = generic_index(word)
        #generic = word.split(':',1)[1].split('[',1)[0]

        declared=False
        for i in range(0, len(declare_ten)):
            if word.split('[',1)[0].split('*',1)[-1] == declare_ten_name[i]:
                if generic == declare_ten_index[i]:
                    # Generic index must be at same position as name it belongs to - dangerous!
                    # Load previous tensor
                    declared=True
                    break
            else:
                continue

        if not declared:
            # Add result to global list
            declare_ten.append(word.split('*',1)[-1])
            declare_ten_index.append(generic)
            declare_ten_name.append(word.split('[',1)[0].split('*',1)[-1])


def print_loop(line, words):
    global tab
    global old_loop

    tab = False
    if ("J:eee" in line or "K:eee" in line):
        for i in range(0, len(words)):
            if ("J:eee" in words[i] or "K:eee" in words[i]):
                tmp = words[i].split('[',1)[1].split(']',1)[0]
                loop = "for ["+tmp[3:4]+"]:"
                if (loop != old_loop):
                    print(loop, file=out)
                old_loop = loop
                tab = True
                break


def generic_index(tensor):
    # Construct generic index representation of specific tensor index
    # Used to check if tensors index the same space, even if they have
    # different indices

    index=list(tensor[tensor.find("[")+1:tensor.find("]")])

    hole     = ['i','j','k','l','m','n','o']
    particle = ['a','b','c','d','e','f','g','h']
    valence  = ['p','q','r','s','t','u','v','w','x','y','z']

    gen=[]
    for i in range (0,len(index)):
        if index[i] in particle:
            gen.append('e')
        elif index[i] in valence:
            gen.append('a')
        elif index[i] in hole:
            gen.append('c')
    return gen


def declare_existing_tensors(declare_list, name, tensor, energy=False):
    c=":"
    if (energy):
        c=""

    print(file=f2)
    print("// " + name, file=f2)
    for i in range(0, len(declare_list)):
        if (tensor + c in declare_list[i]):
            tmp_ten = declare_list[i] + ", " + declare_list[i].split('[',1)[0]

            print("tensor:", tmp_ten, file=f2)


# Check if exchange of a tensor (T[apij] - T[apji]) needs to be loaded
# inside a loop. This is the case if the loop index is
# the index that is permuted so as to generate T[apji]
def load_extra_loop_tensors(line, result):

    words=line.split()

    # Find the loop index
    loop = ''
    for i in range(0, len(words)):
        if ("J:eee" in words[i] or "K:eee" in words[i]):
            tmp = words[i].split('[',1)[1].split(']',1)[0]
            loop = tmp[3:4]

    # See if the loop index is in the same position for the (T - T^+) tensors
    if '(' in words[2] and '(' in words[5]:
        if "eaac" not in words[2] and "eaca" not in words[2]:
            if (words[2].split('(')[1].find(loop) != words[4].find(loop)):
                result.append(words[4].replace(')',''))
        if "eaac" not in words[5] and "eaca" not in words[5]:
            if (words[5].split('(')[1].find(loop) != words[7].find(loop)):
                result.append(words[7].replace(')',''))
    elif '(' in words[3]:
        if "eaac" not in words[3] and "eaca" not in words[3]:
            if (words[3].split('(')[1].find(loop) != words[5].find(loop)):
                result.append(words[5].replace(')',''))
    elif '(' in words[2] :
        if "eaac" not in words[2] and "eaca" not in words[2]:
            if (words[2].split('(')[1].find(loop) != words[4].find(loop)):
                result.append(words[4].replace(')',''))


def print_code_block(code_block, gecco_dir, output):
    with open (gecco_dir+'/itf_python/code_blocks/'+code_block, 'r') as reader:
        print(reader.read(), file=output)


# =========================================================================================
# Main program starts here
# =========================================================================================
import argparse     # Parse arguments
import datetime     # Get tima and date
import time         # For timings
import tempfile     # For temporary files
import os           # For environment variables

# Parse arguments from gecco
parser = argparse.ArgumentParser(
                description="""Process ITF binary contraction file and output ITF algo file
                """,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF binary contraction file')
parser.add_argument('-o','--output',default=None,help='ITF algo file')
parser.add_argument('--multi',dest='multi',action='store_true')
parser.add_argument('--no-multi',dest='multi',action='store_false')
parser.add_argument('--kext',dest='kext',action='store_true')
parser.add_argument('--no-kext',dest='kext',action='store_false')
parser.add_argument('--init-res',dest='initalise',action='store_true')
parser.add_argument('--no-init-res',dest='initalise',action='store_false')
parser.set_defaults(multi=True)
parser.set_defaults(kext=False)
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
outp = args.output
multi = args.multi
kext = args.kext
initalise = args.initalise

# Open bcontr.tmp file to read from
f=open(inp,"r")
# Open first output file
#out=open(outp, "w+")
output=open(outp, "w+")

# Open tempfile to catch INTpp
kext_temp = tempfile.TemporaryFile(mode='w+t')

# Open tempfile for Initalist_Residual
init_res_temp = tempfile.TemporaryFile(mode='w+t')
init_res = False
init_alloc = False

# Declare lists needed in program
prev_lines=[]       # Previous intermediate lines which belong to next result block
prev_inter=[]       # List of previous intemediates, used to alloc/drop
prev_res='#####'    # Result of previous line
prev_generic=[]     # Previous generic result

declare_res=[]            # Global list of result tensors
declare_index=[]          # Global list of result indices
declare_name=[]           # Global of just result tensor names
# TODO: Now all residuals are declared for the mutli cases, don't need all of them, all of the time
declare_res_multi=[
        "R:I[I]",
        "R:ac[pi]","R:ec[ai]","R:ea[ap]",
        "R:aacc[pqij]","R:aaac[pqri]",
        "R:eacc[apij]","R:eaac[apqi]","R:eaca[apiq]","R:eaaa[apqr]",
        "R:eecc[abij]","R:eeac[abpi]","R:eeaa[abpq]"]      # List of all residuals in multireference case
declare_amp_multi=[
        "T:I[I]",
        "T:ac[pi]","T:ec[ai]","T:ea[ap]",
        "T:aacc[pqij]","T:aaac[pqri]",
        "T:eacc[apij]","T:eaac[apqi]","T:eaca[apiq]","T:eaaa[apqr]",
        "T:eecc[abij]","T:eeac[abpi]","T:eeaa[abpq]"]      # List of all residuals in multireference case

declare_inter=[]        # Global list of intermediates
declare_inter_index=[]  # Global list of intermediates
declare_inter_name=[]   # Global list of intermediates

declare_ten=[]          # Global list of tensors involved in binary contractions
declare_ten_index=[]    # Global list of tensor indicies
declare_ten_name=[]     # Global list of tensor names

prev_K4E_lines={0:"start"}
K4E_count=1

tab = False         # True if codes need to be indented with '    '
old_loop = ''       # Stores value of the previous for[x] loop, used so as not to repeat the for loop

# Spin summed family = One or more equations the arise from the spin summation
# proccedure on one result tensor contraction line
begin=False         # Marks the start of a spin summed family of contractions
end=False           # Marks the end of a spin summed family of contractions
old_spin_iter=[]    # Stores list of intermediates used throughout the spin summed family

# Use to redirect output into the main output file or a temperary file
# This can be expanded upon if different parts need to end up in different ---code blocks
special_begin=False
special_end=False
# instead of the above, use CODE_BLOCKS
code_blocks = []
code_block_tmp = {}
current_code_block = ''
# Don't store tensors from these code blocks; instead these are printed out at the end
dont_store=False

# Flag that code block changed
switched_block=False

# Flag if dealing with triples
triples = False

# Read each line of bcontr.tmp and process it
for line_o in f:

    # Try and organise code better
    #tensor_line = itf_line(line_o, out)
    #tensor_line.rename_line()

    #out = output

    # Change names of external tensors (add : + generic index to name)
    line = change_line(line_o)
    words=line.split()

    # Check if we are in a special block which will end up in its own ---code block
    # In this case we are checking for the INTpp interemediate that is passed to Kext
    if (words[0]=='BEGIN_INTPP'):
        special_begin = True
        special_end = False
        continue

    if (words[0]=='END_INTPP'):
        special_begin = False
        special_end = True
        dont_store=True
        continue

    # switch code block
    if (words[0]=='CODE_BLOCK:'):
        # quickly check that this name is no duplicate
        new_block = True
        for block_name in code_blocks:
            if (words[1]==block_name):
                new_block = False

        # keep info about previous temp file to finish last block on that file
        if (current_code_block!=''):
            out_prev = out
            switched_block = True
        else:
            switched_block = False

        # set current code block and switch output
        current_code_block = words[1]
        if (new_block):
            # add to list and open corresp. temp. file
            code_blocks.append(current_code_block)
            code_block_tmp[current_code_block] = tempfile.TemporaryFile(mode='w+t')

        out = code_block_tmp[current_code_block]

        continue

    # make sure that at this point a code block is assigned
    if (current_code_block==''):
        print("Error in python ITF processor: Must start with CODE_BLOCK statement")
        exit(1)
#    # Decide which file to write to
#    if (special_begin and not special_end):
#        out = kext_temp
#    else:
#        out = output


    # Check for spin summed block
    if (words[0]=='BEGIN'):
        # BEGIN marks the begining of a spin summed series of contractions

        begin=True
        end=False
        old_loop = ''
        continue

    if (words[0]=='END'):
        # END marks the end of a spin summed series of contractions

        if old_spin_iter and begin:
            print("drop ", end="", flush=True, file=out)
            print(*list(reversed(old_spin_iter)),sep=", ", file=out)

        old_spin_iter=[]
        end=True
        begin=False
        continue


    # Catch 4-external integrals and replace line with K4E intermediate
    # Note: this might not be general enough to handle all cases
    # ie. with brackets
    if ("K:eeee" in line):
        for j in prev_K4E_lines:
            if (line==prev_K4E_lines[j]):
                # Already printed this line, just change name
                for i in range(0, len(words)):
                    if ("K:eeee" in words[i]):
                        if ("*" in words[i]):
                            # There is a factor which we should retain
                            words[i] = words[i].split("*",1)[0] + "*" + "K4E" + str(j) + ":eecc[abij]"
                        else:
                            words[i] = "K4E" + str(j) + ":eecc[abij]"

                        words[2] = words[i]
                words = words[:3]


        if (line not in prev_K4E_lines.values()):
            for i in range(0, len(words)):
                if ("K:eeee" in words[i]):
                    if ("*" in words[i]):
                        # There is a factor which we should retain
                        words[i] = words[i].split("*",1)[0] + "*" + "K4E" + str(K4E_count) + ":eecc[abij]"
                    else:
                        words[i] = "K4E" + str(K4E_count) + ":eecc[abij]"

                    words[2] = words[i]
            words = words[:3]
            prev_K4E_lines.update({K4E_count:line})
            K4E_count=K4E_count+1

        print("// Replacing line with 4-external integrals", file=out)
        print("// ", line, file=out)
        line = " ".join(words)

    # Check if the line should be in the Initalise_Residual algo. This corresponds
    # to the first iteration where all amps = 0
    init_res = False
    init_inter_res = False
    if (("R:" in line or "G:" in line) and not "STIN" in line and initalise):
        if ("K:" in line or "f" in line):
            init_res = True
    if ("ITIN" in line and not "STIN" in line and initalise):
        if ("K:" in line or "f" in line):
            init_inter_res = True


    # Check if brackets in the binary contraction
    if (len(words)>=4):
        if len(words)==5:
            # I += (K - J)
            if '(' in words[2] and '(' in words[4]:
                add_to_global(words[2].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[4].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
        else:
            if '(' in words[2] and '(' in words[5]:
                # Check if first tensor needs to be declared
                add_to_global(words[2].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[4].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
                # Check if second tensor needs to be declared
                add_to_global(words[5].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[7].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
            elif '(' in words[3]:
                # Check if first tensor needs to be declared
                add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[3].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[5].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
            elif '(' in words[2] :
                # Check if first tensor needs to be declared
                add_to_global(words[2].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[4].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
                add_to_global(words[5],declare_ten,declare_ten_index,declare_ten_name)
            else:
                # No brackets
                # Check if first tensor needs to be declared
                add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)
                # Check if second tensor needs to be declared
                add_to_global(words[3],declare_ten,declare_ten_index,declare_ten_name)
    elif (len(words)==3):
        # Simple add or assign line
        if ("K:[]" not in words[2]):
            # Don't want to declare a tensor for the reference energy
            add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)


    # Process the line, print it out and decide what to alloc/load/drop/store
    if "TIN" in words[0]:  ## <<--- this fragment seems too generic to distinguish form other tensor names
        # The result of the line is an intermediate
        # Check if contraction forms and intermediate
        prev_lines.append(line)

        # Add intermediate name to prev_inter list, this is then used to alloc/drop intermediate
        # tensors. Hack(?) used to avoid adding the an intermediate to the list if it has already
        # been added. This is no longer safe if two intermediates need to be allocated, have the
        # same name, yet a different index. I assume this will not happen.
        if len(prev_inter)==0:
            prev_inter.append(words[0].replace('.',''))
        else:
            found = False
            for i in range(0, len(prev_inter)):
                if words[0].split('[',1)[0].replace('.','') in prev_inter[i].split('[',1)[0].replace('.',''):
                    generic1=generic_index(words[0])
                    generic2=generic_index(prev_inter[i])
                    if (generic1==generic2):
                        found = True

            if (not found):
                prev_inter.append(words[0].replace('.',''))


        if words[0].replace('.','') not in declare_inter:
            # Add intermedate to global list
            generic=generic_index(words[0])

            declared=False
            for i in range(0, len(declare_inter)):
                if words[0].split('[',1)[0].replace('.','') == declare_inter_name[i]:
                    if generic == declare_inter_index[i]:
                        # Generic index must be at same position as name it belongs to - dangerous!
                        # TODO: use dict instead of declare_index and declare_name
                        # Load previous tensor
                        declared=True
                        break
                else:
                    continue

            if not declared:
                # Add result to global list
                declare_inter.append(words[0].replace('.',''))
                declare_inter_index.append(generic)
                declare_inter_name.append(words[0].split('[',1)[0].replace('.',''))

    else:
        # Convert to generic index
        # Compare generic indicies to see if belongs in same block
        # + Load instead of alloc

        generic=generic_index(words[0])

        if "R:eaac" in prev_res and "R:eaca" in words[0] or "R:eaca" in prev_res and "R:eaac" in words[0]:

            # Still within the same result block (ie. still with the same residual)
            # This is a special case of the code below - R:eaac and R:eaca share the same intermediates and
            # so should be consider the 'same' residual. This means both of this will form one residual
            # block.
            # The code below is repeated from the end of else statment. In the old days, you would use a
            # goto instead...

            # Alloc intermediates if needed
            if prev_inter:
                print("alloc ", end="", flush=True, file=out)
                print(*prev_inter,sep=", ", file=out)

            # Print intermediates and load/drop relavant tensors
            print_inter(prev_lines, init_inter_res)

            # Print loop for 3-external integrals
            if (begin and not end):
                print_loop(line, words)

            # Print result line
            #print_result(line, tab)
            print_result(line, tab, init_res)
            tab = False

            # Drop intermediates if needed, don't drop if needed again in
            # the next contraction which is part of the same spin summed family
            if prev_inter and not begin:
                print("drop ", end="", flush=True, file=out)
                print(*list(reversed(prev_inter)),sep=", ", file=out)

            # Store intermediate so as to drop at the end of spin summed family
            if begin and not end:
                if not old_spin_iter:
                    old_spin_iter=prev_inter

            prev_res=words[0]
            prev_lines=[]
            prev_inter=[]

        else:
            if words[0] != prev_res and generic != prev_generic:# and "R:eaac" not in prev_res and "R:eaca" not in words[0]:
                # Next result is different from previous, so close off block

                if prev_res != '#####':
                    # Add generic index to residual tensor name
                    if ("R[" in prev_res):
                        prev_res = prev_res.replace("R[", "R:" + "".join(generic_index(prev_res)) + "[")
                    elif ("G[" in prev_res):
                        prev_res = prev_res.replace("G[", "G:" + "".join(generic_index(prev_res)) + "[")

                    # dont_store obsolete, assume False always
                    # Don't print out store for speical code blocks (ie. not residuals)
                    if not dont_store:

                        # store the special case of R:eaac and R:eaca
                        if "R:eaac" in prev_res or "R:eaca" in prev_res:
                            prev_res = "R:eaac[apqi], R:eaca[apiq]"

                        if (not switched_block):
                            print("store", prev_res.replace('.',''), file=out)
                            print(file=out)
                        else:
                            # write to previous temp file if code block switched
                            print("store", prev_res.replace('.',''), file=out_prev)
                            print(file=out_prev)
                            switched_block = False # reset flag

                        if (init_alloc and initalise):
                            print("store", prev_res.replace('.',''), file=init_res_temp)
                            print(file=init_res_temp)
                            init_alloc = False
                    # Finished special block, so we can start storeing again
                    dont_store=False

                # Check whether to load previously allocated tensor
                # To load, the tensor name and generic index associated with it must be equal
                loaded=False
                for i in range(0, len(declare_name)):
                    if words[0].split('[',1)[0].replace('.','') == declare_name[i]:
                        if generic == declare_index[i]:
                            # Generic index must be at same position as name it belongs to - dangerous!
                            # Load previous tensor
                            tmp_res=words[0]

                            if ("R[" in words[0]):
                                tmp_res = words[0].replace("R[", "R:" + "".join(generic_index(words[0])) + "[")

                            elif ("G[" in words[0]):
                                tmp_res = words[0].replace("G[", "G:" + "".join(generic_index(words[0])) + "[")

                            # load the special case of R:eaac and R:eaca
                            if "R:eaac" in tmp_res or "R:eaca" in tmp_res:
                                tmp_res = "R:eaca[apiq], R:eaac[apqi]"

                            print("load", tmp_res.replace('.',''), file=out)
                            loaded=True
                            break
                    else:
                        continue

                if not loaded:
                    # Add result to global list and alloc new result
                    declare_res.append(words[0].split('.',1)[1])
                    declare_index.append(generic)
                    declare_name.append(words[0].split('[',1)[0].replace('.',''))

                    # If allocating R:eaac OR R:eaca, then also add its partner to
                    # declare_res (both will be allocated at the same time)
                    if "R:eaac" in words[0]:
                        declare_res.append("R:eaca[apiq]")
                        declare_index.append(['e','a','c','a'])
                        declare_name.append("R:eaca")

                    elif "R:eaca" in words[0]:
                        declare_res.append("R:eaac[apqi]")
                        declare_index.append(['e','a','a','c'])
                        declare_name.append("R:eaac")

                    # Add generic index to residual tensor name, ie. R:eecc
                    tmp_res=words[0]
                    if ("R[" in words[0]):
                        tmp_res = words[0].replace("R[", "R:" + "".join(generic_index(words[0])) + "[")
                    elif ("G[" in words[0]):
                        tmp_res = words[0].replace("G[", "G:" + "".join(generic_index(words[0])) + "[")

                    # allocate the special case of R:eaac and R:eaca
                    if "R:eaac" in tmp_res or "R:eaca" in tmp_res:
                        tmp_res = "R:eaca[apiq], R:eaac[apqi]"

                    print("alloc", tmp_res.replace('.',''), file=out)
                    if (init_res and initalise):
                        print("alloc", tmp_res.replace('.',''), file=init_res_temp)
                        init_alloc = True

                # Alloc intermediates if needed
                if prev_inter:
                    print("alloc ", end="", flush=True, file=out)
                    print(*prev_inter,sep=", ", file=out)

                # Print intermediates and load/drop relavant tensors
                print_inter(prev_lines, init_inter_res)

                # Print result line
                #print_result(line)
                print_result(line, False, init_res)

                # Drop intermediates if needed
                if prev_inter and not begin:
                    print("drop ", end="", flush=True, file=out)
                    print(*list(reversed(prev_inter)),sep=", ", file=out)

                # Store intermediate so as to drop at the end of spin summed family
                if begin and not end:
                    if not old_spin_iter:
                        old_spin_iter=prev_inter

                prev_res=words[0]
                prev_lines=[]
                prev_inter=[]

            #elif words[0] == prev_res:
            else:
                # Still within the same result block (ie. still with the same residual)

                # Alloc intermediates if needed
                if prev_inter:
                    print("alloc ", end="", flush=True, file=out)
                    print(*prev_inter,sep=", ", file=out)

                # Print intermediates and load/drop relavant tensors
                print_inter(prev_lines, init_inter_res)

                # Print loop for 3-external integrals
                if (begin and not end):
                    print_loop(line, words)

                # Print result line
                #print_result(line, tab)
                print_result(line, tab, init_res)
                tab = False

                # Drop intermediates if needed, don't drop if needed again in
                # the next contraction which is part of the same spin summed family
                if prev_inter and not begin:
                    print("drop ", end="", flush=True, file=out)
                    print(*list(reversed(prev_inter)),sep=", ", file=out)

                # Store intermediate so as to drop at the end of spin summed family
                if begin and not end:
                    if not old_spin_iter:
                        old_spin_iter=prev_inter

                prev_res=words[0]
                prev_lines=[]
                prev_inter=[]

        # Update generic index for next line
        prev_generic=generic


# Close off final result block
print("store", prev_res.replace('.',''), file=out)

#Close input file
f.close()

#output.close()

# should not be needed as now everything is on temp files
# Open and write file again so as to prepend the declaration of tensors
#f2=open(outp, "r")
#tmp=f2.read()
#f2.close()

f2=open(outp, "w")
gecco = os.environ["GECCO_DIR"]


print_code_block('header', gecco, f2)
now = datetime.datetime.now()
print("// Created on:", now.strftime("%d-%m-%Y %H:%M"), file=f2)
print(file=f2)

# Print index-spaces
if (not multi):
    print_code_block('single_ref/index_spaces', gecco, f2)
else:
    print_code_block('multi_ref/index_spaces', gecco, f2)

# Print out tensors used in contractions
# For now, delcare them on disk
declare_ten.sort()
declare_ten.sort(key=len)
for i in range(0, len(declare_ten)):
    if ("T:" in declare_ten[i] or "K:" in declare_ten[i] or "K4E" in declare_ten[i] or "f:" in declare_ten[i]\
        or "Dm" in declare_ten[i] or "J:" in declare_ten[i] or "INTpp" in declare_ten[i] or "Ym" in declare_ten[i]): continue
    if ("[]" in declare_ten[i]):
        print("tensor:", declare_ten[i] + ",  !Create{type:scalar}", file=f2)
    else:
        print("tensor:", declare_ten[i] + ",  !Create{type:disk}", file=f2)

# Print already existing tensor, ie. don't need !Create{type:disk}
declare_existing_tensors(declare_ten, "K-integral tensors", "K")
declare_existing_tensors(declare_ten, "J-integral tensors", "J")

# If not declared already, declare integrals needed in calculating fock matricies etc.
combined = '\t'.join(declare_ten)
if multi:
    if ("INTkx:eecc" not in combined):
        print("tensor: INTkx:eecc[abij], !Create{type:disk}", file=f2)
    if ("INTkx:eeac" not in combined):
        print("tensor: INTkx:eeac[abpj], !Create{type:disk}", file=f2)
    if ("INTkx:eeaa" not in combined):
        print("tensor: INTkx:eeaa[abpq], !Create{type:disk}", file=f2)
    if "f:ac" not in combined:
        print("tensor: f:ac[pi], f:ac", file=f2)
    if "f:ea" not in combined:
        print("tensor: f:ea[ap], f:ea", file=f2)
    if "f:aa" not in combined:
        print("tensor: f:aa[pq], f:aa", file=f2)
    if "f:ee" not in combined:
        print("tensor: f:ee[ab], f:ee", file=f2)
    if "f:cc" not in combined:
        print("tensor: f:cc[ij], f:cc", file=f2)
    if "f:ec" not in combined:
        print("tensor: f:ec[ai], f:ec", file=f2)
    if "J:eacc" not in combined:
        print("tensor: J:eacc[apij], J:eacc", file=f2)
    if "J:ecca" not in combined:
        print("tensor: J:ecca[aijp], J:ecca", file=f2)
    if "J:ecaa" not in combined:
        print("tensor: J:ecaa[aipq], J:ecaa", file=f2)
    if "J:eaaa" not in combined:
        print("tensor: J:eaaa[apqr], J:eaaa", file=f2)
    if "J:ccaa" not in combined:
        print("tensor: J:ccaa[ijpq], J:ccaa", file=f2)
    if "J:caaa" not in combined:
        print("tensor: J:caaa[ipqr], J:caaa", file=f2)
    if "J:ccca" not in combined:
        print("tensor: J:ccca[ijkp], J:ccca", file=f2)
    if "J:eccc" not in combined:
        print("tensor: J:eccc[aijk], J:eccc", file=f2)
    if "J:eecc" not in combined:
        print("tensor: J:eecc[abij], J:eecc", file=f2)
    if "J:eeaa" not in combined:
        print("tensor: J:eeaa[abpq], J:eeaa", file=f2)
    if "K:eeaa" not in combined:
        print("tensor: K:eeaa[abpq], K:eeaa", file=f2)
    if "K:ccaa" not in combined:
        print("tensor: K:ccaa[ijpq], K:ccaa", file=f2)
    if "K:ecaa" not in combined:
        print("tensor: K:ecaa[aipq], K:ecaa", file=f2)
    if "K:cccc" not in combined:
        print("tensor: K:cccc[ijkl], K:cccc", file=f2)
    if "K:aaaa" not in combined:
        print("tensor: K:aaaa[pqrs], K:aaaa", file=f2)
    print("", file=f2)

declare_existing_tensors(declare_ten, "Special integral tensors", "K4E",True)
if multi:
    if ("K4E:eecc" not in combined):
        print("tensor: K4E:eecc[abij], K4E:eecc", file=f2)
    if ("K4E:eeac" not in combined):
        print("tensor: K4E:eeac[abpj], K4E:eeac", file=f2)
    if ("K4E:eeaa" not in combined):
        print("tensor: K4E:eeaa[abpq], K4E:eeaa", file=f2)
    print("tensor: K4C[abmn], K4C", file=f2)

# if -kext flag was set
if (kext):
    declare_existing_tensors(declare_ten, "Tensor to send to Kext", "INTpp",True)
    if multi:
        print("tensor: INTpp[abmn], INTpp", file=f2)
        print("tensor: INTpp1[abmi], !Create{type:disk}, INTpp1", file=f2)
        print("tensor: INTpp2[abmq], !Create{type:disk}, INTpp2", file=f2)
        print("", file=f2)
        print("tensor: deltaai[pm], DeltaActInt", file=f2)
        print("tensor: deltaci[im], DeltaCloInt", file=f2)
else:
    print(file=f2)
    if multi:
        print("// Tensor to send to Kext", file=f2)
        print("tensor: INTpp[abmn], INTpp", file=f2)
        print("tensor: INTpp1[abmi], !Create{type:disk}, INTpp1", file=f2)
        print("", file=f2)
        print("tensor: deltaai[pm], DeltaActInt", file=f2)
        print("tensor: deltaci[im], DeltaCloInt", file=f2)
    else:
        print("// Tensor to send to Kext", file=f2)
        print("tensor: INTpp[abij], INTpp", file=f2)
declare_existing_tensors(declare_ten, "Fock tensors", "f")
declare_existing_tensors(declare_ten, "Amplitude tensors", "T")
#if (multi): print("tensor: R[I],  R:I", file=f2)
declare_existing_tensors(declare_res, "Residual tensors", "R")

g_residual = False
if "G:eecc[abij]" in declare_res:
    g_residual = True
    print("tensor: R:eecc[abij], R:eecc", file=f2)

if (any('R:eeeccc' in s for s in declare_res)):
   triples = True

if (multi):
    print(file=f2)
    print("// Residuals not used in the code", file=f2)
    for i in range(0, len(declare_res_multi)):
        if any(declare_res_multi[i].split('[',1)[0]+"[" in s for s in declare_res):
            continue
        else:
            tmp_ten = declare_res_multi[i] + ", " + declare_res_multi[i].split('[',1)[0]
            print("tensor:", tmp_ten, file=f2)

    print(file=f2)
    print("// Amplitudes not used in the code", file=f2)
    for i in range(0, len(declare_amp_multi)):
        if any(declare_amp_multi[i].split('[',1)[0]+"[" in s for s in declare_ten):
            continue
        else:
            tmp_ten = declare_amp_multi[i] + ", " + declare_amp_multi[i].split('[',1)[0]
            print("tensor:", tmp_ten, file=f2)


declare_existing_tensors(declare_res, "Energy and DIIS scalars", "ECC", True)

# Check if we have singles amplitudes
singles = False
for i in range(0, len(declare_ten)):
    if ("T:ec" in declare_ten[i]):
        singles = True
        break

if singles:
    if multi:
        print_code_block('multi_ref/declare_tensors_singles', gecco, f2)
    else:
        print_code_block('single_ref/declare_tensors_singles', gecco, f2)

if multi:
    print_code_block('multi_ref/declare_tensors', gecco, f2)

    #if "Ym1" in combined:
    print("tensor: Ym1[pp],      !Create{type:disk} Ym1",file=f2)
    #if "Ym2" in combined:
    print("tensor: Ym2[pppp],    !Create{type:disk} Ym2",file=f2)
    #if "Ym3" in combined:
    print("tensor: Ym3[pppppp],  !Create{type:disk} Ym3",file=f2)
    if "Ym4" in combined:
        print("tensor: Ym4[pppppppp],!Create{type:disk} Ym4",file=f2)
    print(file=f2)

else:
    print_code_block('single_ref/declare_tensors', gecco, f2)


# Print out intermediates
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    if ("[]" in declare_inter[i]):
        print("tensor: %-18s !Create{type:scalar}" % (declare_inter[i] + ","), file=f2)
    else:
        print("tensor: %-18s !Create{type:plain}" % (declare_inter[i] + ","), file=f2)
if g_residual:
    print("tensor: G:eecc[abij], !Create{type:disk}", file=f2)


# Print out code blocks
# Need to initalise the amplitudes first
print(file=f2)
print(file=f2)
print('---- code("Init_Amplitudes")',file=f2)
if (multi):
    for i in range(0, len(declare_ten)):
        if ("T:" in declare_ten[i]):
            print("alloc", declare_ten[i], file=f2)
            print("store", declare_ten[i], file=f2)
else:
    # Initalise amplitudes using MP2
    if singles:
        print_code_block('single_ref/init_amplitudes_singles', gecco, f2)
    else:
        print_code_block('single_ref/init_amplitudes', gecco, f2)


# Calculate the reference energy for single-reference methods
if not multi:
    print_code_block('single_ref/ref_energy', gecco, f2)


# Print out INTpp update
if (not kext):
    print(file=f2)
    print(file=f2)
    print('---- code("Update_INTkx")', file=f2)
    print("// Intermediate to pass to Kext", file=f2)
    if multi:
        print("alloc INTpp1[abmi]", file=f2)
        print("alloc INTpp[abmn]", file=f2)
        print("load deltaci[im], deltaai[pm]", file=f2)
        print("load T:eecc[abij], T:eeac[abpi]", file=f2)
        print(".INTpp1[abmj] += T:eecc[abij] deltaci[im]", file=f2)
        print(".INTpp1[abmi] += T:eeac[abpi] deltaai[pm]", file=f2)
        print(".INTpp[abmn] += INTpp1[abmj] deltaci[jn]", file=f2)
        print("drop T:eeac[abpi], T:eecc[abij]", file=f2)
        print("drop deltaai[pm], deltaci[im]", file=f2)
        print("store INTpp[abmn]", file=f2)
        print("store INTpp1[abmi]", file=f2)
    else:
        print("alloc INTpp[abij]",file=f2)
        print("load T:eecc[abij]",file=f2)
        print(".INTpp[abij] := T:eecc[abij]",file=f2)
        print("drop T:eecc[abij]",file=f2)
        print("store INTpp[abij]",file=f2)
# Now treated via CODE_BLOCK:
#else:
#    kext_temp.seek(0)
#    for line in kext_temp:
#        print(line.strip(), file=f2)

#    print("store INTpp[abij]",file=f2)


# Transform K ext from internal indicies to closed and active and back
if multi:
    print(file=f2)
    print_code_block('multi_ref/transform_k', gecco, f2)
    print_code_block('multi_ref/transform_intk', gecco, f2)


# Print out Init_Residual
#if (initalise):
#    print(file=f2)
#    print(file=f2)
#    print('---- code("Init_Residual")', file=f2)
#    init_res_temp.seek(0)
#    for line in init_res_temp:
#        print(line.strip(), file=f2)

#if multi:
#    print(file=f2)
#    print(file=f2)
#    print('---- code("Init_Residual")', file=f2)
#    combined = '\t'.join(declare_res)
#    if "R:eecc" in combined:
#        print("alloc R:eecc[abij]", file=f2)
#        print("load K:eecc[abij]", file=f2)
#        print(".R:eecc[abij] += K:eecc[abij]", file=f2)
#        print("drop K:eecc[abij]", file=f2)
#        print("store R:eecc[abji]", file=f2)
#        print("", file=f2)
#    if "R:eeac" in combined:
#        print("alloc R:eeac[bapi]", file=f2)
#        print("load K:eeac[baqi], Ym1[qp]", file=f2)
#        print(".R:eeac[bapi] -= K:eeac[baqi] Ym1[qp]", file=f2)
#        print("drop Ym1[qp], K:eeac[baqi]", file=f2)
#        print("store R:eeac[bapi]", file=f2)
#        print("", file=f2)
#    if "R:eacc" in combined:
#        print("alloc R:eacc[apij]", file=f2)
#        print("load K:eacc[apij]", file=f2)
#        print(".R:eacc[apij] += K:eacc[apij]", file=f2)
#        print("drop K:eacc[apij]", file=f2)
#        print("load Ym1[pq], K:eacc[aqij]", file=f2)
#        print(".R:eacc[apij] -= Ym1[pq] K:eacc[aqij]", file=f2)
#        print("drop K:eacc[aqij], Ym1[pq]", file=f2)
#        print("store R:eacc[apij]", file=f2)
#    if "R:ec" in combined:
#        print("", file=f2)
#        print("alloc R:ec[ai]", file=f2)
#        print("load f:ec[ai]", file=f2)
#        print(".R:ec[ai] += f:ec[ai]", file=f2)
#        print("drop f:ec[ai]", file=f2)
#        print("load Ym1[pq], K:eaca[aqip], K:eaac[aqpi]", file=f2)
#        print(".R:ec[ai] += Ym1[pq] (K:eaca[aqip] - K:eaac[aqpi])", file=f2)
#        print("drop K:eaac[aqpi], K:eaca[aqip], Ym1[pq]", file=f2)
#        print("load Ym1[pq], K:eaca[aqip]", file=f2)
#        print(".R:ec[ai] += Ym1[pq] K:eaca[aqip]", file=f2)
#        print("drop K:eaca[aqip], Ym1[pq]", file=f2)
#        print("store R:ec[ai]", file=f2)
#    if "R:ea" in combined:
#        print("", file=f2)
#        print("alloc R:ea[ap]", file=f2)
#        print("load f:ea[aq], Ym1[qp]", file=f2)
#        print(".R:ea[ap] += f:ea[aq] Ym1[qp]", file=f2)
#        print("drop Ym1[qp], f:ea[aq]", file=f2)
#        print("store R:ea[ap]", file=f2)
#    if "R:ac" in combined:
#        print("", file=f2)
#        print("alloc R:ac[pi]", file=f2)
#        print("load f:ac[pi]", file=f2)
#        print(".R:ac[pi] += f:ac[pi]", file=f2)
#        print("drop f:ac[pi]", file=f2)
#        print("load Ym1[pq], f:ac[qi]", file=f2)
#        print(".R:ac[pi] -= Ym1[pq] f:ac[qi]", file=f2)
#        print("drop f:ac[qi], Ym1[pq]", file=f2)
#        print("load Ym1[qr], K:aaac[rpqi]", file=f2)
#        print(".R:ac[pi] += Ym1[qr] (K:aaac[rpqi] - K:aaac[rqpi])", file=f2)
#        print(".R:ac[pi] += Ym1[qr] K:aaac[rpqi]", file=f2)
#        print("drop K:aaac[rpqi], Ym1[qr]", file=f2)
#        print("store R:ac[pi]", file=f2)
#    if "R:eaca" in combined:
#        print("", file=f2)
#        print("alloc R:eaca[apiq]", file=f2)
#        print("load Ym1[pq], f:ec[ai]", file=f2)
#        print(".R:eaca[apiq] -= Ym1[pq] f:ec[ai]", file=f2)
#        print("drop f:ec[ai], Ym1[pq]", file=f2)
#        print("load K:eaca[apir], Ym1[rq]", file=f2)
#        print(".R:eaca[apiq] -= K:eaca[apir] Ym1[rq]", file=f2)
#        print("drop Ym1[rq], K:eaca[apir]", file=f2)
#        print("store R:eaca[apiq]", file=f2)
#    if "R:eaac" in combined:
#        print("", file=f2)
#        print("alloc R:eaac[apqi]", file=f2)
#        print("load K:eaac[apqi], Ym1[rq]", file=f2)
#        print(".R:eaac[apqi] += K:eaac[apri] Ym1[rq]", file=f2)
#        print("drop Ym1[rq], K:eaac[apqi]", file=f2)
#        print("store R:eaac[apqi]", file=f2)

# Print out all code blocks generated from input
for block_name in code_blocks:
    code_block_tmp[block_name].seek(0)
    print('---- code("'+block_name+'")', file=f2)
    for line in code_block_tmp[block_name]:
        print(line.rstrip(), file=f2)

    code_block_tmp[block_name].close()

## Print out residual equations
#if (not initalise): print(file=f2)
#print(file=f2)
#print('---- code("Residual")', file=f2)
#f2.write(tmp)

# JOSH: the following seems to belong to a code block
# Symmetrise tensors
if "G:eecc[abij]" in declare_res:
    print("", file=f2)
    print("alloc R:eecc[abij]", file=f2)
    print("load G:eecc[abij]", file=f2)
    print(".R:eecc[abij] += G:eecc[abij]", file=f2)
    print(".R:eecc[abij] += G:eecc[baji]", file=f2)
    print("drop G:eecc[abij]", file=f2)
    print("store R:eecc[abij]", file=f2)

#if "R[pqij]" in declare_res:
#    print("load R:aacc[pqij]", file=f2)
#    print(".R:aacc[pqij] += R:aacc[qpji]", file=f2)
#    print("store R:aacc[pqij]", file=f2)
#    print(file=f2)
#
#if "R[abpq]" in declare_res:
#    print("load R:eeaa[abpq]", file=f2)
#    print(".R:eeaa[abpq] += R:eeaa[baqp]", file=f2)
#    print("store R:eeaa[abpq]", file=f2)
#    print(file=f2)

print(file=f2)
print(file=f2)


# Print out Generate_Fock_Matricies
if multi:
    print_code_block('multi_ref/generate_fock_matricies', gecco, f2)


# Print out Update_Amplitudes
if multi:
    print_code_block('multi_ref/form_dm3', gecco, f2)
    # OLD:
    print_code_block('multi_ref/transform_residual', gecco, f2)
    print_code_block('multi_ref/create_amplitude_update', gecco, f2)
    print_code_block('multi_ref/update_amplitudes', gecco, f2)
    # NEW:
    print_code_block('multi_ref/amplitude_update',gecco, f2)
    print_code_block('multi_ref/transform_to_pair_index',gecco, f2)
    #
    print_code_block('multi_ref/construct_gs_overlap', gecco, f2)
    print_code_block('multi_ref/construct_s2', gecco, f2)
    print_code_block('multi_ref/construct_projected_s2', gecco, f2)
    print_code_block('multi_ref/construct_offdiag_x', gecco, f2)
else:
    if singles:
        print_code_block('single_ref/update_amplitudes_singles', gecco, f2)
    else:
        print_code_block('single_ref/update_amplitudes', gecco, f2)


print("---- end", file=f2)

# Close off files (including temporary ones)
f2.close()
kext_temp.close()
init_res_temp.close()
