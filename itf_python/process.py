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
        # with there names defined in C++, ie. R:eecc
        names = ["R", "G", "T", "K", "f"]
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
        else:
            load_ten="load "

        # Construct load line
        for tensor in result[:-1]:
            load_ten = load_ten + tensor + ', '
        load_ten = load_ten + result[-1]

        # Print out load, line and drop
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
    if ("K:eeec" in line or "K:eeea" in line):
        for i in range(0, len(words)):
            if ("K:eeec" in words[i] or "K:eeea" in words[i]):
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
    valence  = ['p','q','r','s','t','u','v','w']

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


#def rename_integrals(line):
#    # Rename K -> J and permute indicies to [eecc]
#    # Molpro stores integrals in two sets of arrays: K and J. Not to be confused with the exchange
#    # and Coulomb integrals of HF theory; there is no relationship.
#    # For integrals like K_{ak}^{ic} = K_{ac}^{ik}, these are store on the K arrays.
#    # For integrals K_{ka}^{ic}, these are store on the J arrays.
#
#    second = False
#
#    words = line.split()
#    for i in range(0, len(words)):
#        if ("K:ceec" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[2:3] + tmp[1:2] + tmp[0:1] + tmp[3:4]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[2:3] + tmp2[1:2] + tmp2[0:1] + tmp2[3:4]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:ccee" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[2:3] + tmp[3:4] + tmp[0:1] + tmp[1:2]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[2:3] + tmp2[3:4] + tmp2[0:1] + tmp2[1:2]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:ecec" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[0:1] + tmp[2:3] + tmp[1:2] + tmp[3:4]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[0:1] + tmp2[2:3] + tmp2[1:2] + tmp2[3:4]
#            words[i] = words[i].split(':',1)[0].replace("K", "J") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:ecce" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[0:1] + tmp[3:4] + tmp[2:3] + tmp[1:2]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[0:1] + tmp2[3:4] + tmp2[2:3] + tmp2[1:2]
#            #words[i] = "(" + words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1] + " + J:" + tmp2 + "[" + tmp + "])"
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        # 3 internal integrals
#        # Pairing as above (1+3, 2+4)
#        elif ("K:ccec" in words[i]):
#            if (second):
#                # Second integral in the line, so must be the 'transpose', ie. J:eccc
#                tmp = words[i].split('[',1)[1].split(']',1)[0]
#                tmp = tmp[2:3] + tmp[0:1] + tmp[1:2] + tmp[3:4]
#                tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#                tmp2 = tmp2[2:3] + tmp2[1:2] + tmp2[0:1] + tmp2[3:4]
#                words[i] = words[i].split(':',1)[0].replace("K", "J") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#                second = False
#            else:
#                tmp = words[i].split('[',1)[1].split(']',1)[0]
#                tmp = tmp[2:3] + tmp[1:2] + tmp[0:1] + tmp[3:4]
#                tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#                tmp2 = tmp2[2:3] + tmp2[1:2] + tmp2[0:1] + tmp2[3:4]
#                words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#                second = True
#
#        elif ("K:ccce" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[3:4] + tmp[0:1] + tmp[1:2] + tmp[2:3]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[3:4] + tmp2[0:1] + tmp2[1:2] + tmp2[2:3]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:cecc" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[1:2] + tmp[0:1] + tmp[3:4] + tmp[2:3]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[1:2] + tmp2[0:1] + tmp2[3:4] + tmp2[2:3]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        # 3 internal integrals which arise from a permuation, therefore they become J:eccc
#        elif ("KP:ccce" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[3:4] + tmp[1:2] + tmp[0:1] + tmp[2:3]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[3:4] + tmp2[1:2] + tmp2[0:1] + tmp2[2:3]
#            words[i] = words[i].split(':',1)[0].replace("KP", "J") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("KP:ccec" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[2:3] + tmp[0:1] + tmp[1:2] + tmp[3:4]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[2:3] + tmp2[1:2] + tmp2[0:1] + tmp2[3:4]
#            words[i] = words[i].split(':',1)[0].replace("KP", "J") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#
#        # 3 external integrals
#        # Index pariing is different compared to normal 1+3, 2+4;
#        # Indicies are paired between the 1+2, 3+4, so need to convert index which is pair above to this convention
#        elif ("K:ceee" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[1:2] + tmp[3:4] + tmp[2:3] + tmp[0:1]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[1:2] + tmp2[3:4] + tmp2[2:3] + tmp2[0:1]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:ecee" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[0:1] + tmp[2:3] + tmp[3:4] + tmp[1:2]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[0:1] + tmp2[2:3] + tmp2[3:4] + tmp2[1:2]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("K:eece" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[1:2] + tmp[3:4] + tmp[0:1] + tmp[2:3]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[1:2] + tmp2[3:4] + tmp2[0:1] + tmp2[2:3]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        # 3 external integrals which arise from a permuation, therefore they become J:eeec
#        elif ("KP:ceee" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[1:2] + tmp[3:4] + tmp[2:3] + tmp[0:1]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[1:2] + tmp2[3:4] + tmp2[2:3] + tmp2[0:1]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#            words[i] = words[i].split(':',1)[0].replace("KP", "K") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#        elif ("KP:ecee" in words[i]):
#            tmp = words[i].split('[',1)[1].split(']',1)[0]
#            tmp = tmp[0:1] + tmp[2:3] + tmp[3:4] + tmp[1:2]
#            tmp2 = words[i].split(':',1)[1].split('[',1)[0]
#            tmp2 = tmp2[0:1] + tmp2[2:3] + tmp2[3:4] + tmp2[1:2]
#            words[i] = words[i].split(':',1)[0] + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#            words[i] = words[i].split(':',1)[0].replace("KP", "K") + ":" + tmp2 + "[" + tmp + "]" + words[i].split(']',1)[1]
#
#    line = " ".join(words)
#    return line


# =========================================================================================
# Main program starts here
# =========================================================================================
import argparse     # Parse arguments
import datetime     # Get tima and date
import time         # For timings
import tempfile     # For temporary files

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
        "R:eacc[apij]","R:eaac[apqi]","R:eaaa[apqr]",
        "R:eecc[abij]","R:eeac[abpi]","R:eeaa[abpq]"]      # List of all residuals in multireference case
declare_amp_multi=[
        "T:I[I]",
        "T:ac[pi]","T:ec[ai]","T:ea[ap]",
        "T:aacc[pqij]","T:aaac[pqri]",
        "T:eacc[apij]","T:eaac[apqi]","T:eaaa[apqr]",
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
# Don't store tensors from these code blocks; instead these are printed out at the end
dont_store=False

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

    # Decide which file to write to
    if (special_begin and not special_end):
        out = kext_temp
    else:
        out = output


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


    # Catch if K[ceec] or K[ecec] is on the line. Need to permute the index to [eecc],
    # but also change name to J[eecc]
    #if ("K:" in line or "KP:"):
    #    if ("K:eecc" not in line or "K:eccc" not in line or "K:eecc" not in line or "K:cccc" not in line):
    #        line = rename_integrals(line)
    #        words=line.split()


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
    if "TIN" in words[0]:
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

        if words[0] != prev_res and generic != prev_generic:
            # Next result is different from previous, so close off block

            if prev_res != '#####':
                # Add generic index to residual tensor name
                if ("R[" in prev_res):
                    prev_res = prev_res.replace("R[", "R:" + "".join(generic_index(prev_res)) + "[")
                elif ("G[" in prev_res):
                    prev_res = prev_res.replace("G[", "G:" + "".join(generic_index(prev_res)) + "[")

                # Don't print out store for speical code blocks (ie. not residuals)
                if (not dont_store):
                    print("store", prev_res.replace('.',''), file=out)
                    print(file=out)

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

                # Add generic index to residual tensor name, ie. R:eecc
                tmp_res=words[0]
                if ("R[" in words[0]):
                    tmp_res = words[0].replace("R[", "R:" + "".join(generic_index(words[0])) + "[")
                elif ("G[" in words[0]):
                    tmp_res = words[0].replace("G[", "G:" + "".join(generic_index(words[0])) + "[")

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
print("store", prev_res.replace('.',''), file=output)
output.close()
f.close()

# Open and write file again so as to prepend the declaration of tensors
f2=open(outp, "r")
tmp=f2.read()
f2.close()

f2=open(outp, "w")

print("// This ITF algo file was created using the GeCCo ITF translator", file=f2)
print("// Author: J.A. Black", file=f2)
print(file=f2)
now = datetime.datetime.now()
print("// Created on:", now.strftime("%d-%m-%Y %H:%M"), file=f2)
print(file=f2)

# Declare tensors and index-spaces
print("---- decl", file=f2)
if (not multi):
    print("index-space: ijklmn, Closed  , c", file=f2)
    print("index-space: abcdef, External, e", file=f2)
    print("index-space: CD, Core, C", file=f2)
else:
    print("index-space: pqrstuvw, Active  , a", file=f2)
    print("index-space: ijkl    , Closed  , c", file=f2)
    print("index-space: gh      , ClosedF , f", file=f2)
    print("index-space: abcd    , External, e", file=f2)
    print("index-space: mno     , Internal, i", file=f2)
    print(file=f2)
    print("index-space: I       , ConfigI0, I", file=f2)

print(file=f2)

# Print out tensors used in contractions
# For now, delcare them on disk
declare_ten.sort()
declare_ten.sort(key=len)
for i in range(0, len(declare_ten)):
    if ("T:" in declare_ten[i] or "K:" in declare_ten[i] or "K4E" in declare_ten[i] or "f:" in declare_ten[i]\
        or "Dm" in declare_ten[i] or "J:" in declare_ten[i] or "INTpp" in declare_ten[i]): continue
    if ("[]" in declare_ten[i]):
        print("tensor:", declare_ten[i] + ",  !Create{type:scalar}", file=f2)
    else:
        print("tensor:", declare_ten[i] + ",  !Create{type:disk}", file=f2)

# Print already existing tensor, ie. don't need !Create{type:disk}
declare_existing_tensors(declare_ten, "K-integral tensors", "K")
declare_existing_tensors(declare_ten, "J-integral tensors", "J")
declare_existing_tensors(declare_ten, "Special integral tensors", "K4E",True)
if (kext):
    declare_existing_tensors(declare_ten, "Tensor to send to Kext", "INTpp",True)
else:
    print(file=f2)
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

if (not multi):
    # Tensors needed in CCD
    print("tensor: ERef[], ERef     // Reference energy", file=f2)
    if (singles):
        print("tensor: EMp1[], EMp1     // MP2 energy", file=f2)
        print("tensor: EDi1[], EDi1     // Direct 1st order energy", file=f2)
        print("tensor: Nrm1[], Nrm1     // Singles amplitude norm", file=f2)
        print("tensor: Var1[], Var1     // Singles residual norm", file=f2)
    print("tensor: EMp2[], EMp2     // MP2 energy", file=f2)
    print("tensor: EDi2[], EDi2     // Direct 1st order energy", file=f2)
    print("tensor: Nrm2[], Nrm2     // Doubles amplitude norm", file=f2)
    print("tensor: Var2[], Var2     // Doubles residual norm", file=f2)
    if (triples):
        print("tensor: EDi3[], EDi3     // Direct 1st order energy", file=f2)
        print("tensor: Nrm3[], Nrm3     // Singles amplitude norm", file=f2)
        print("tensor: Var3[], Var3     // Singles residual norm", file=f2)
    print(file=f2)
    print("// Tensors needed to calculate the reference energy", file=f2)
    print("tensor: f:CC[CC],   f:CC", file=f2)
    print("tensor: CoreH[ii],  h:cc", file=f2)
    print("tensor: CoreH[CC],  h:CC", file=f2)
    print("tensor: Delta[ii],  Delta", file=f2)
    print("tensor: DeltaC[CC], DeltaC", file=f2)

# Declare density and overlap tensors
if (multi):
    print(file=f2)
    print("// Fock matricies, including core orbitals",file=f2)
    print("tensor: fc:aa[pq], fc:aa", file=f2)
    print("tensor: fc:cc[ki], fc:cc", file=f2)
    print("tensor: fc:ea[aq], fc:ea", file=f2)
    print("tensor: fc:ec[ai], fc:ec", file=f2)
    print("tensor: fc:ee[ab], fc:ee", file=f2)
    print("tensor: fc:ca[ip], fc:ca", file=f2)
    print("", file=f2)
    print("tensor: J:eacc[apij]", file=f2)
    print("tensor: J:ecca[aijp]", file=f2)
    print("tensor: J:ecaa[aipq]", file=f2)
    print("tensor: J:eaaa[apqr]", file=f2)
    print("tensor: J:ccaa[ijpq]", file=f2)
    print("tensor: J:caaa[ipqr]", file=f2)
    print("tensor: J:ccca[ijkq]", file=f2)
    print("tensor: K:ccaa[ijpq]", file=f2)
    print("tensor: K:ecaa[aipq]", file=f2)
    print("", file=f2)
    print("// Effective Fock matricies",file=f2)
    print("tensor: g:aa[pq]", file=f2)
    print("tensor: g:ac[pi]", file=f2)
    print("tensor: g:cc[ij]", file=f2)
    print("tensor: g:ea[ap]", file=f2)
    print("tensor: g:ec[ai]", file=f2)
    print("tensor: g:ee[ab]", file=f2)


    print(file=f2)
    print("// Resuced density tensors (Icc-Icc coupling-coefficients)",file=f2)
    print("// these are created on C++ side in CreateMrciTensors method",file=f2)
    print("tensor: Dm1[pp],         DDm1",file=f2)
    print("tensor: Dm2[pppp],       DDm2",file=f2)
    print("tensor: Dm3[pppppp],     DDm3",file=f2)
    print("tensor: Dm1H[pp],        DHm1",file=f2)
    print("tensor: Dm2H[pppp],      DHm2",file=f2)
    print("tensor: Dm3H[pppppp],    DHm3",file=f2)
    print(file=f2)

    print("// Non-disk density matrix drivers",file=f2)
    print("// Can be loaded, but can not be stored.",file=f2)
    print("// spec:<Ref|+-+-|Ref> means that Dm2X[pqrs] = <Ref|E^p_q R^r_s|Ref>",file=f2)
    print("tensor: Dm2X[pppp],      !Create{type:cc-drv; spec:<Ref|+-+-|Ref>; irrep:0;}",file=f2)
    print("tensor: Dm3X[pppppp],    !Create{type:cc-drv; spec:<Ref|+-+-+-|Ref>; irrep:0;}",file=f2)
    print("tensor: Dm2HX[pppp],     !Create{type:cc-drv; spec:<Ref|/+/-/+/-|Ref>; irrep:0;}",file=f2)
    print("tensor: Dm3HX[pppppp],   !Create{type:cc-drv; spec:<Ref|/+/-/+/-/+/-|Ref>; irrep:0;}",file=f2)
    print(file=f2)

    print("// Delta tensors",file=f2)
    print("tensor: deltacc[ij],      DeltaCloClo", file=f2)
    print("tensor: deltaaa[pq],      DeltaActAct",file=f2)
    print("tensor: delta4[pppp],     !Create{type:plain}, Delta4  // Intermediate rank4 delta tensor",file=f2)
    print(file=f2)

    print("// Overlap tensors labeled by the exciation class",file=f2)
    print("tensor: S1:I1[pp],       S1:I1",file=f2)
    print("tensor: S2:I1[pppppp],   S2:I1",file=f2)
    print("tensor: S3:I1[pppp],     S3:I1",file=f2)
    print("tensor: S1:I2[pppp],     S1:I2",file=f2)
    print("tensor: S1:S0[pp],       S1:S0",file=f2)
    print("tensor: S2:S0[pppppp],   S2:S0",file=f2)
    print("tensor: S3:S0[pppp],     S3:S0",file=f2)
    print("tensor: S2:S1[pppp],     S2:S1",file=f2)
    print("tensor: S3:S1[pp],       S3:S1",file=f2)
    print("tensor: S1:S2[pp],       S1:S2",file=f2)
    print("tensor: S1:P0[pppp],     S1:P0",file=f2)
    print("tensor: S1:P1[pp],       S1:P1",file=f2)
    print(file=f2)
    print("// Orthogonal residuals",file=f2)
    print("tensor: OR:I[I], OR:I", file=f2)
    print("tensor: OR:ac[pi], OR:ac", file=f2)
    print("tensor: OR:ec[ai], OR:ec", file=f2)
    print("tensor: OR:ea[ap], OR:ea", file=f2)
    print("tensor: OR:aacc[pqij], OR:aacc", file=f2)
    print("tensor: OR:aaac[pqri], OR:aaac", file=f2)
    print("tensor: OR:eacc[apij], OR:eacx", file=f2)
    print("tensor: OR:eaac[apqi], OR:eaac", file=f2)
    print("tensor: OR:eaaa[apqr], OR:eaaa", file=f2)
    print("tensor: OR:eecc[abij], OR:eecc", file=f2)
    print("tensor: OR:eeac[abpi], OR:eeac", file=f2)
    print("tensor: OR:eeaa[abpq], OR:eeaa", file=f2)
    print(file=f2)
    print("// Orthogonal amplitudes",file=f2)
    print("tensor: OT:I[I], OT:I", file=f2)
    print("tensor: OT:ac[pi], OT:ac", file=f2)
    print("tensor: OT:ec[ai], OT:ec", file=f2)
    print("tensor: OT:ea[ap], OT:ea", file=f2)
    print("tensor: OT:aacc[pqij], OT:aacc", file=f2)
    print("tensor: OT:aaac[pqri], OT:aaac", file=f2)
    print("tensor: OT:eacc[apij], OT:eacx", file=f2)
    print("tensor: OT:eaac[apqi], OT:eaac", file=f2)
    print("tensor: OT:eaaa[apqr], OT:eaaa", file=f2)
    print("tensor: OT:eecc[abij], OT:eecc", file=f2)
    print("tensor: OT:eeac[abpi], OT:eeac", file=f2)
    print("tensor: OT:eeaa[abpq], OT:eeaa", file=f2)

# Declare denominator shifts
print(file=f2)
if (singles):
    print("tensor: ShiftS[], ShiftS", file=f2)
print("tensor: ShiftP[], ShiftP", file=f2)

# Print out intermediates
print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    if ("[]" in declare_inter[i]):
        print("tensor: %-18s !Create{type:scalar}" % (declare_inter[i] + ","), file=f2)
    else:
        print("tensor: %-18s !Create{type:plain}" % (declare_inter[i] + ","), file=f2)
if g_residual:
    print("tensor: G:eecc[abij], !Create{type:disk}", file=f2)


# Intermediates needed in single-reference amp update
if (not multi):
    print(file=f2)
    if (singles):
        print("tensor: L1[ai],            !Create{type:plain}", file=f2)
    print("tensor: L2[abij],          !Create{type:plain}", file=f2)
    if (triples):
        print("tensor: L3[abcijk],        !Create{type:plain}", file=f2)
    print("tensor: C[abij],           !Create{type:plain}", file=f2)

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
    print("// Using MP2 amplitudes for starting guess", file=f2)
    if (singles):
        print("alloc T:ec[ai], EMp1[], Nrm1[]", file=f2)
        print("load f:ec[ai], f:ee[aa], f:cc[ii]", file=f2)
        print(".T:ec[ai] -= f:ec[ai]", file=f2)
        print("denom-scale T:ec[ai], f:ee[aa] - f:cc[ii]", file=f2)
        print(".EMp1[] += T:ec[ai] f:ec[ai]", file=f2)
        print(".Nrm1[] += 2.0*T:ec[ai] T:ec[ai]", file=f2)
        print("drop f:cc[ii], f:ee[aa], f:ec[ai]", file=f2)
        print("store Nrm1[], EMp1[], T:ec[ai]", file=f2)
        print("", file=f2)
    print("alloc EMp2[], Nrm2[]", file=f2)
    print("for [i,j]:", file=f2)
    print("   alloc T:eecc[abij]", file=f2)
    print("   load K:eecc[**ij], f:ee[aa], f:cc[ii], f:cc[jj]", file=f2)
    print("   .T:eecc[abij] -= K:eecc[abij]", file=f2)
    print("   denom-scale T:eecc[abij], f:ee[aa] + f:ee[bb] - f:cc[ii] - f:cc[jj]", file=f2)
    print("   .EMp2[] += (2.0*T:eecc[abij] - T:eecc[baij]) K:eecc[abij]", file=f2)
    print("   .Nrm2[] += (2.0*T:eecc[abij] - T:eecc[baij]) T:eecc[abij]", file=f2)
    print("   drop f:cc[jj], f:cc[ii], f:ee[aa], K:eecc[**ij]", file=f2)
    print("   store T:eecc[**ij]", file=f2)
    print("store Nrm2[], EMp2[]", file=f2)
    if (triples):
        print("", file=f2)
        print("alloc T:eeeccc[abcijk]", file=f2)
        print("store T:eeeccc[abcijk]", file=f2)

# Calculate the reference energy for single-reference methods
if (not multi):
    print("", file=f2)
    print("", file=f2)
    print('---- code ("Ref_Energy")', file=f2)
    print("alloc ERef[]", file=f2)
    print("// Closed-shell contribution", file=f2)
    print("load f:cc[ii], CoreH[ii], Delta[ii]", file=f2)
    print(".ERef += f:cc[ij] Delta[ij]", file=f2)
    print(".ERef += CoreH[ij] Delta[ij]", file=f2)
    print("drop Delta, CoreH, f:cc", file=f2)
    print("", file=f2)
    print("// Core contribution", file=f2)
    print("load f:CC[CC], CoreH[CC], DeltaC[CC]", file=f2)
    print(".ERef += f:CC[CD] DeltaC[CD]", file=f2)
    print(".ERef += CoreH[CD] DeltaC[CD]", file=f2)
    print("drop DeltaC, CoreH, f:CC", file=f2)
    print("store ERef[]", file=f2)


# Print out INTpp update
print(file=f2)
print(file=f2)
print('---- code("Update_Kext_Tensor")', file=f2)
if (not kext):
    print("// Intermediate to pass to Kext", file=f2)
    print("alloc INTpp[abij]",file=f2)
    print("load T:eecc[abij]",file=f2)
    print(".INTpp[abij] := T:eecc[abij]",file=f2)
    print("drop T:eecc[abij]",file=f2)
    print("store INTpp[abij]",file=f2)
else:
    kext_temp.seek(0)
    for line in kext_temp:
        print(line.strip(), file=f2)

    print("store INTpp[abij]",file=f2)


# Print out Init_Residual
if (initalise):
    print(file=f2)
    print(file=f2)
    print('---- code("Init_Residual")', file=f2)
    init_res_temp.seek(0)
    for line in init_res_temp:
        print(line.strip(), file=f2)


# Print out residual equations
if (not initalise): print(file=f2)
print(file=f2)
print('---- code("Residual")', file=f2)
f2.write(tmp)

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



if (multi):
    print(file=f2)
    print(file=f2)
    print('---- code("Generate_Fock_Matrices")', file=f2)
    print('alloc f:ac[**]', file=f2)
    print('load fc:ca[**]', file=f2)
    print('.f:ac[pi] := fc:ca[ip]', file=f2)
    print('drop fc:ca[**]', file=f2)
    print('load deltacc[**], J:ccca[****]', file=f2)
    print('.f:ac[pi] += (2.*J:ccca[jkip] - J:ccca[ikjp]) deltacc[kj]', file=f2)
    print('drop J:ccca[****], deltacc[**]', file=f2)
    print('store f:ac[**]', file=f2)
    print('', file=f2)
    print('alloc f:ec[**]', file=f2)
    print('load fc:ec[**]', file=f2)
    print('.f:ec[ai] += fc:ec[ai]', file=f2)
    print('drop fc:ec[**]', file=f2)
    print('load deltacc[**], J:eccc[****]', file=f2)
    print('.f:ec[ai] += (2.*J:eccc[aijk] - J:eccc[ajik]) deltacc[kj]', file=f2)
    print('drop J:eccc[****], deltacc[**]', file=f2)
    print('store f:ec[**]', file=f2)
    print('', file=f2)
    print('alloc f:cc[**]', file=f2)
    print('load fc:cc[**]', file=f2)
    print('.f:cc[ij] := fc:cc[ij]', file=f2)
    print('drop fc:cc[**]', file=f2)
    print('load deltacc[**], K:cccc[****]', file=f2)
    print('.f:cc[ij] += (2.*K:cccc[jkil] - K:cccc[jkli]) deltacc[lk]', file=f2)
    print('drop K:cccc[****], deltacc[**]', file=f2)
    print('store f:cc[**]', file=f2)
    print('', file=f2)
    print('alloc f:ea[**]', file=f2)
    print('load fc:ea[**]', file=f2)
    print('.f:ea[ap] := fc:ea[ap]', file=f2)
    print('drop fc:ea[**]', file=f2)
    print('load deltacc[**]', file=f2)
    print('load J:eacc[****], J:ecca[****]', file=f2)
    print('.f:ea[ap] += (2.*J:eacc[apij] - J:ecca[ajip])  deltacc[ji]', file=f2)
    print('drop J:ecca[****], J:eacc[****]', file=f2)
    print('drop deltacc[**]', file=f2)
    print('store f:ea[**]', file=f2)
    print('', file=f2)
    print('alloc f:ee[**]', file=f2)
    print('load fc:ee[**]', file=f2)
    print('.f:ee[ab] := fc:ee[ab]', file=f2)
    print('drop fc:ee[**]', file=f2)
    print('load deltacc[**]', file=f2)
    print('load J:eecc[****], K:eecc[****]', file=f2)
    print('.f:ee[ab] += (2.*J:eecc[baij] - K:eecc[baji]) deltacc[ji]', file=f2)
    print('drop K:eecc[****], J:eecc[****]', file=f2)
    print('drop deltacc[**]', file=f2)
    print('store f:ee[**]', file=f2)
    print('', file=f2)
    print('alloc f:aa[**]', file=f2)
    print('load fc:aa[**]', file=f2)
    print('.f:aa[pq] := fc:aa[pq]', file=f2)
    print('drop fc:aa[**]', file=f2)
    print('load deltacc[**]', file=f2)
    print('load J:ccaa[****], K:ccaa[****]', file=f2)
    print('.f:aa[pq] += (2.*J:ccaa[ijqp] - K:ccaa[jiqp]) deltacc[ji]', file=f2)
    print('drop K:ccaa[****], J:ccaa[****]', file=f2)
    print('drop deltacc[**]', file=f2)
    print('store f:aa[**]', file=f2)
    print('', file=f2)
    print('alloc g:ac[**]', file=f2)
    print('load f:ac[**]', file=f2)
    print('.g:ac[pi] := f:ac[pi]', file=f2)
    print('drop f:ac[**]', file=f2)
    print('load J:caaa[****], Dm1[**]', file=f2)
    print('.g:ac[pi] += (J:caaa[ipqr] - .5*J:caaa[irqp]) Dm1[qr]', file=f2)
    print('drop Dm1[**], J:caaa[****]', file=f2)
    print('store g:ac[**]', file=f2)
    print('', file=f2)
    print('alloc g:ec[**]', file=f2)
    print('load f:ec[**]', file=f2)
    print('.g:ec[ai] := f:ec[ai]', file=f2)
    print('drop f:ec[**]', file=f2)
    print('load Dm1[**]', file=f2)
    print('load J:ecaa[****], K:ecaa[****]', file=f2)
    print('.g:ec[ai] += (J:ecaa[aipq] - .5*K:ecaa[aipq]) Dm1[pq]', file=f2)
    print('drop K:ecaa[****], J:ecaa[****]', file=f2)
    print('drop Dm1[**]', file=f2)
    print('store g:ec[**]', file=f2)
    print('', file=f2)
    print('alloc g:cc[**]', file=f2)
    print('load f:cc[**]', file=f2)
    print('.g:cc[ij] := f:cc[ij]', file=f2)
    print('drop f:cc[**]', file=f2)
    print('load Dm1[**]', file=f2)
    print('load J:ccaa[****], K:ccaa[****]', file=f2)
    print('.g:cc[ij] += (J:ccaa[jipq] - .5*K:ccaa[jiqp]) Dm1[pq]', file=f2)
    print('drop K:ccaa[****], J:ccaa[****]', file=f2)
    print('drop Dm1[**]', file=f2)
    print('store g:cc[**]', file=f2)
    print('', file=f2)
    print('alloc g:ea[**]', file=f2)
    print('load f:ea[**]', file=f2)
    print('.g:ea[ap] := f:ea[ap]', file=f2)
    print('drop f:ea[**]', file=f2)
    print('load J:eaaa[****], Dm1[**]', file=f2)
    print('.g:ea[ap] += (J:eaaa[apqr] - .5*J:eaaa[arqp]) Dm1[qr]', file=f2)
    print('drop Dm1[**], J:eaaa[****]', file=f2)
    print('store g:ea[**]', file=f2)
    print('', file=f2)
    print('alloc g:ee[**]', file=f2)
    print('load f:ee[**]', file=f2)
    print('.g:ee[ab] := f:ee[ab]', file=f2)
    print('drop f:ee[**]', file=f2)
    print('load Dm1[**]', file=f2)
    print('load J:eeaa[****], K:eeaa[****]', file=f2)
    print('.g:ee[ab] += (J:eeaa[bapq] - .5*K:eeaa[baqp]) Dm1[pq]', file=f2)
    print('drop K:eeaa[****], J:eeaa[****]', file=f2)
    print('drop Dm1[**]', file=f2)
    print('store g:ee[**]', file=f2)
    print('', file=f2)
    print('alloc g:aa[**]', file=f2)
    print('load f:aa[**]', file=f2)
    print('.g:aa[pq] := f:aa[pq]', file=f2)
    print('drop f:aa[**]', file=f2)
    print('load Dm1[**]', file=f2)
    print('load K:aaaa[****]', file=f2)
    print('.g:aa[pq] += (K:aaaa[qrps] - .5*K:aaaa[qrsp])  Dm1[rs]', file=f2)
    print('drop K:aaaa[****]', file=f2)
    print('drop Dm1[**]', file=f2)
    print('store g:aa[**]', file=f2)


# Apply preconditioner and update amplitudes
print(file=f2)
print(file=f2)
print('---- code("Update_Amplitudes")',file=f2)
if (not multi):
    if (singles):
        #print("// Update singles", file=f2)
        #print("load R:ec[ai]", file=f2)
        #print("// L1 = R^i_a/D^i_a; D^i_a = e_a-e_i", file=f2)
        #print("alloc L1[ai]", file=f2)
        #print(".L1[ai] += R:ec[ai]", file=f2)
        #print("denom-scale L1[ai], [1,0]", file=f2)
        #print("", file=f2)
        #print("alloc EDi1[], Nrm1[], Var1[]", file=f2)
        #print("load T:ec[ai], f:ec[ai]", file=f2)
        #print("// Update T:ec",file=f2)
        #print("// T^i_a = T^i_a-R^i_a/D^i_a",file=f2)
        #print(".T:ec[ai] -= L1[ai]", file=f2)
        #print(".EDi1[] += 2.0*T:ec[ai] f:ec[ai]", file=f2)
        #print(".Nrm1[] += 2.0*T:ec[ai] T:ec[ai]", file=f2)
        #print(".Var1[] += 2.0*L1[ai] L1[ai]", file=f2)
        #print("drop f:ec[ai]", file=f2)
        #print("store T:ec[ai]", file=f2)
        #print("store Var1[], Nrm1[], EDi1[]", file=f2)
        #print("", file=f2)
        #print("drop L1[ai], R:ec[ai]", file=f2)
        #print("", file=f2)
        #print("load T:ec[ai]", file=f2)
        print("// Update singles", file=f2)
        print("load R:ec[ai]", file=f2)
        print("load f:ee[aa], f:cc[ii], ShiftS[]", file=f2)
        print("denom-scale R:ec[ai], f:ee[aa] - f:cc[ii] + ShiftS[]", file=f2)
        print("drop ShiftS[], f:cc[ii], f:ee[aa]", file=f2)
        print("", file=f2)
        print("alloc EDi1[], Nrm1[], Var1[]", file=f2)
        print("load T:ec[ai], f:ec[ai]", file=f2)
        print(".T:ec[ai] -= R:ec[ai]", file=f2)
        print(".EDi1[] += 2.0*T:ec[ai] f:ec[ai]", file=f2)
        print(".Nrm1[] += 2.0*T:ec[ai] T:ec[ai]", file=f2)
        print(".Var1[] += 2.0*R:ec[ai] R:ec[ai]", file=f2)
        print("drop f:ec[ai]", file=f2)
        print("store T:ec[ai]", file=f2)
        print("store Var1[], Nrm1[], EDi1[]", file=f2)
        print("drop R:ec[ai]", file=f2)
        print("", file=f2)
    print("// Update doubles", file=f2)
    if (singles): print("load T:ec[ai]", file=f2)
    print("alloc EDi2[], Nrm2[], Var2[]", file=f2)
    print("load R:eecc[abij], K:eecc[abij]", file=f2)
    print("load T:eecc[abij]", file=f2)
    print("load f:ee[aa], f:cc[ii], ShiftP[]", file=f2)
    print("denom-scale R:eecc[abij], f:ee[aa] + f:ee[bb] - f:cc[ii] - f:cc[jj] + ShiftP[]", file=f2)
    print("drop ShiftP[], f:cc[ii], f:ee[aa]", file=f2)
    print("", file=f2)
    print(".T:eecc[abij] -= R:eecc[abij]", file=f2)
    print("", file=f2)
    print("alloc C[abij]", file=f2)
    print(".C[abij] += T:eecc[abij]", file=f2)
    if (singles): print(".C[abij] += T:ec[ai] T:ec[bj]", file=f2)
    print(".EDi2 += (2.0*C[abij] - C[baij]) K:eecc[abij]", file=f2)
    print(".Nrm2 += (2.0*C[abij] - C[baij]) C[abij]", file=f2)
    print("drop C[abij]", file=f2)
    print("", file=f2)
    print(".Var2 += (2.0*R:eecc[abij] - R:eecc[baij]) R:eecc[abij]", file=f2)
    print("", file=f2)
    print("store T:eecc[abij]", file=f2)
    print("drop K:eecc[abij], R:eecc[abij]", file=f2)
    print("store Var2[], Nrm2[], EDi2[]", file=f2)
    if (singles): print("drop T:ec[ai]", file=f2)
    #print("// Update doubles",file=f2)
    #print("alloc EDi2[], Nrm2[], Var2[]",file=f2)
    #print("for [i,j]:",file=f2)
    #print("   load R:eecc[**ij]",file=f2)
    #print("   // L2 = R^{ij}_{ab}/D^{ij}_{ab}; D^{ij}_{ab} = e_a+e_b-e_i-e_j",file=f2)
    #print("   alloc L2[**ij]",file=f2)
    #print("   .L2[**ij] += R:eecc[**ij]",file=f2)
    #print("   denom-scale L2[**ij], [1,1,0,0]",file=f2)
    #print("",file=f2)
    #print("   load K:eecc[**ij]",file=f2)
    #print("   .Var2 += (2.0*L2[abij] - L2[baij]) L2[abij]",file=f2)
    #print("",file=f2)
    #print("   // Update T:eecc",file=f2)
    #print("   // T^{ij}_{ab} = T^{ij}_{ab}-R^{ij}_{ab}/D^{ij}_{ab}",file=f2)
    #print("   load T:eecc[**ij]",file=f2)
    #print("   .T:eecc[abij] -= L2[abij]",file=f2)
    #print("",file=f2)
    #print("   alloc C[**ij]",file=f2)
    #print("   .C[**ij] += T:eecc[**ij]",file=f2)
    #if (singles):
    #    print("   .C[abij] += T:ec[ai] T:ec[bj]",file=f2)
    #print("",file=f2)
    #print("   .EDi2 += (2.0*C[abij] - C[baij]) K:eecc[abij]",file=f2)
    #print("   .Nrm2 += (2.0*C[abij] - C[baij]) C[abij]",file=f2)
    #print("   drop C[**ij]",file=f2)
    #print("   store T:eecc[**ij]",file=f2)
    #print("   drop K:eecc[**ij]",file=f2)
    #print("",file=f2)
    #print("   drop L2[**ij]",file=f2)
    #print("   drop R:eecc[**ij]",file=f2)
    #print("store Var2[], Nrm2[], EDi2[]",file=f2)
    #if (singles): print("drop T:ec[ai]", file=f2)
    if (triples):
        # TODO: Do this properly
        print("", file=f2)
        print("// Update triples", file=f2)
        print("load R:eeeccc[abcijk]", file=f2)
        print("alloc L3[abcijk]", file=f2)
        print(".L3[abcijk] += R:eeeccc[abcijk]", file=f2)
        print("denom-scale L3[abcijk], [1,1,1,0,0,0]", file=f2)
        print("", file=f2)
        print("load Var3[]", file=f2)
        print(".Var3[] += (6.0*L3[abcijk] - L3[bcaijk] - L3[cabijk] - L3[bacijk] - L3[acbijk] - L3[cbaijk]) L3[abcijk]", file=f2)
        print("store Var3[]", file=f2)
        print("", file=f2)
        print("load T:eeeccc[abcijk]", file=f2)
        print(".T:eeeccc[abcijk] -= L3[abcijk]", file=f2)
        print("load Nrm3[]", file=f2)
        print(".Nrm3[] += (6.0*T:eeeccc[abcijk] - T:eeeccc[bcaijk] - T:eeeccc[cabijk] - T:eeeccc[bacijk] - T:eeeccc[acbijk] - T:eeeccc[cbaijk]) T:eeeccc[abcijk]", file=f2)
        print("store Nrm3[]", file=f2)
        print("store T:eeeccc[abcijk]", file=f2)
        print("", file=f2)
        print("drop L3[abcijk]", file=f2)
        print("drop R:eeeccc[abcijk]", file=f2)
else:
    print("alloc Var1[]", file=f2)
    print("", file=f2)
    print("load OR:ac[pi]", file=f2)
    print("load f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:ac[pi], f:aa[pp] - f:cc[ii]", file=f2)
    print("drop f:cc[ii], f:aa[pp]", file=f2)
    print("alloc OT:ac[pi]", file=f2)
    print(".OT:ac[pi] -= OR:ac[pi]", file=f2)
    print("store OT:ac[pi]", file=f2)
    print(".Var1[] += 2.0*OR:ac[pi] OR:ac[pi]", file=f2)
    print("drop OR:ac[pi]", file=f2)
    print("", file=f2)
    print("load OR:ec[ai]", file=f2)
    print("load f:ee[aa], f:cc[ii]", file=f2)
    print("denom-scale OR:ec[ai], f:ee[aa] - f:cc[ii]", file=f2)
    print("drop f:cc[ii], f:ee[aa]", file=f2)
    print("alloc OT:ec[ai]", file=f2)
    print(".OT:ec[ai] -= OR:ec[ai]", file=f2)
    print("store OT:ec[ai]", file=f2)
    print(".Var1[] += 2.0*OR:ec[ai] OR:ec[ai]", file=f2)
    print("drop OR:ec[ai]", file=f2)
    print("", file=f2)
    print("load OR:ea[ap]", file=f2)
    print("load f:ee[aa], f:aa[pp]", file=f2)
    print("denom-scale OR:ea[ap], f:ee[aa] - f:aa[pp]", file=f2)
    print("drop f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:ea[ap]", file=f2)
    print(".OT:ea[ap] -= OR:ea[ap]", file=f2)
    print("store OT:ea[ap]", file=f2)
    print(".Var1[] += 2.0*OR:ea[ap] OR:ea[ap]", file=f2)
    print("drop OR:ea[ap]", file=f2)
    print("", file=f2)
    print("store Var1[]", file=f2)
    print("alloc Var2[]", file=f2)
    print("", file=f2)
    print("load OR:aacc[pqij]", file=f2)
    print("load f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:aacc[pqij], f:aa[pp] + f:aa[qq] - f:cc[ii] - f:cc[jj]", file=f2)
    print("drop f:cc[ii], f:aa[pp]", file=f2)
    print("alloc OT:aacc[pqij]", file=f2)
    print(".OT:aacc[pqij] -= OR:aacc[pqij]", file=f2)
    print("store OT:aacc[pqij]", file=f2)
    print(".Var2[] += (2.0*OR:aacc[pqij] - OR:aacc[qpij]) OR:aacc[pqij]", file=f2)
    print("drop OR:aacc[pqij]", file=f2)
    print("", file=f2)
    print("load OR:aaac[pqri]", file=f2)
    print("load f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:aaac[pqri], f:aa[pp] + f:aa[qq] - f:aa[rr] - f:cc[ii]", file=f2)
    print("drop f:cc[ii], f:aa[pp]", file=f2)
    print("alloc OT:aaac[pqri]", file=f2)
    print(".OT:aaac[pqri] -= OR:aaac[pqri]", file=f2)
    print("store OT:aaac[pqri]", file=f2)
    print(".Var2[] += OR:aaac[pqri] OR:aaac[pqri]", file=f2)
    print("drop OR:aaac[pqri]", file=f2)
    print("", file=f2)
    print("load OR:eacc[apij]", file=f2)
    print("load f:ee[aa], f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:eacc[apij], f:ee[aa] + f:aa[pp] - f:cc[ii] - f:cc[jj]", file=f2)
    print("drop f:cc[ii], f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:eacc[apij]", file=f2)
    print(".OT:eacc[apij] -= OR:eacc[apij]", file=f2)
    print("store OT:eacc[apij]", file=f2)
    print(".Var2[] += OR:eacc[apij] OR:eacc[apij]", file=f2)
    print("drop OR:eacc[apij]", file=f2)
    print("", file=f2)
    print("load OR:eaac[apqi]", file=f2)
    print("load f:ee[aa], f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:eaac[apqi], f:ee[aa] + f:aa[pp] - f:aa[qq] - f:cc[ii]", file=f2)
    print("drop f:cc[ii], f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:eaac[apqi]", file=f2)
    print(".OT:eaac[apqi] -= OR:eaac[apqi]", file=f2)
    print("store OT:eaac[apqi]", file=f2)
    print(".Var2[] += OR:eaac[apqi] OR:eaac[apqi]", file=f2)
    print("drop OR:eaac[apqi]", file=f2)
    print("", file=f2)
    print("load OR:eaaa[apqr]", file=f2)
    print("load f:ee[aa], f:aa[pp]", file=f2)
    print("denom-scale OR:eaaa[apqr], f:ee[aa] + f:aa[pp] - f:aa[qq] - f:aa[rr]", file=f2)
    print("drop f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:eaaa[apqr]", file=f2)
    print(".OT:eaaa[apqr] -= OR:eaaa[apqr]", file=f2)
    print("store OT:eaaa[apqr]", file=f2)
    print(".Var2[] += OR:eaaa[apqr] OR:eaaa[apqr]", file=f2)
    print("drop OR:eaaa[apqr]", file=f2)
    print("", file=f2)
    print("load OR:eecc[abij]", file=f2)
    print("load f:ee[aa], f:cc[ii]", file=f2)
    print("denom-scale OR:eecc[abij], f:ee[aa] + f:ee[bb] - f:cc[ii] - f:cc[jj]", file=f2)
    print("drop f:cc[ii], f:ee[aa]", file=f2)
    print("alloc OT:eecc[abij]", file=f2)
    print(".OT:eecc[abij] -= OR:eecc[abij]", file=f2)
    print("store OT:eecc[abij]", file=f2)
    print(".Var2[] += (2.0*OR:eecc[abij] - OR:eecc[baij]) OR:eecc[abij]", file=f2)
    print("drop OR:eecc[abij]", file=f2)
    print("", file=f2)
    print("load OR:eeac[abpi]", file=f2)
    print("load f:ee[aa], f:aa[pp], f:cc[ii]", file=f2)
    print("denom-scale OR:eeac[abpi], f:ee[aa] + f:ee[bb] - f:aa[pp] - f:cc[ii]", file=f2)
    print("drop f:cc[ii], f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:eeac[abpi]", file=f2)
    print(".OT:eeac[abpi] -= OR:eeac[abpi]", file=f2)
    print("store OT:eeac[abpi]", file=f2)
    print(".Var2[] += OR:eeac[abpi] OR:eeac[abpi]", file=f2)
    print("drop OR:eeac[abpi]", file=f2)
    print("", file=f2)
    print("load OR:eeaa[abpq]", file=f2)
    print("load f:ee[aa], f:aa[pp]", file=f2)
    print("denom-scale OR:eeaa[abpq], f:ee[aa] + f:ee[bb] - f:aa[pp] - f:aa[qq]", file=f2)
    print("drop f:aa[pp], f:ee[aa]", file=f2)
    print("alloc OT:eeaa[abpq]", file=f2)
    print(".OT:eeaa[abpq] -= OR:eeaa[abpq]", file=f2)
    print("store OT:eeaa[abpq]", file=f2)
    print(".Var2[] += (2.0*OR:eeaa[abpq] - OR:eeaa[bapq]) OR:eeaa[abpq]", file=f2)
    print("drop OR:eeaa[abpq]", file=f2)
    print("", file=f2)
    print("store Var2[]", file=f2)


# Print out code needed to evaluate the overlap matrix
if (multi):
    print(file=f2)
    print("// Set up 3rd order denisty and hole tensors",file=f2)
    print("// This is taken from the cic code",file=f2)
    print('---- code("FormDm3OnDisk")',file=f2)
    print("// <E^pqr_stu> += <E^p_s E^q_t R^r_u>",file=f2)
    print("//             -= delta_rs <E^pq_ut>",file=f2)
    print("//             -= delta_rt <E^pq_su>",file=f2)
    print("//             -= delta_qs <E^p_t E^r_u>",file=f2)
    print("",file=f2)
    print("alloc Dm3[******]",file=f2)
    print("load deltaaa[**]",file=f2)
    print("load Dm3X[******]",file=f2)
    print(".Dm3[pqruts] += Dm3X[psqtru]",file=f2)
    print("drop Dm3X",file=f2)
    print("load Dm2[****]",file=f2)
    print(".Dm3[pqruts] -= deltaaa[rs] Dm2[pqtu]",file=f2)
    print(".Dm3[pqruts] -= deltaaa[rt] Dm2[pqus]",file=f2)
    print("drop Dm2",file=f2)
    print("load Dm2X[****]",file=f2)
    print(".Dm3[pqruts] -= deltaaa[qs] Dm2X[ptru]",file=f2)
    print("drop Dm2X",file=f2)
    print("drop deltaaa[**]",file=f2)
    print("store Dm3",file=f2)
    print("",file=f2)
    print("// also hole dm3",file=f2)
    print("// <E^k_s E^l_t E^m_u E^r_m E^q_l E^p_k> = HDm3[sturqp]",file=f2)
    print("//          + <E^k_s E^p_k E^l_t E^q_l E^m_u E^r_m>",file=f2)
    print("//          - delta_pu <E^m_s E^l_t E^q_l E^r_m>",file=f2)
    print("//          - delta_pt <E^l_s E^m_u E^r_m E^q_l>",file=f2)
    print("//          - delta_qu <E^k_s E^p_k E^m_t E^r_m>",file=f2)
    print("",file=f2)
    print("alloc Dm3H[******]",file=f2)
    print("load deltaaa[**]",file=f2)
    print("load Dm3HX[******]",file=f2)
    print(".Dm3H[sturqp] += Dm3HX[sptqur]",file=f2)
    print("drop Dm3HX",file=f2)
    print("load Dm2H[****]",file=f2)
    print(".Dm3H[sturqp] -= deltaaa[pu] Dm2H[stqr]",file=f2)
    print(".Dm3H[sturqp] -= deltaaa[pt] Dm2H[surq]",file=f2)
    print("drop Dm2H",file=f2)
    print("load Dm2HX[****]",file=f2)
    print(".Dm3H[sturqp] -= deltaaa[qu] Dm2HX[sptr]",file=f2)
    print("drop Dm2HX",file=f2)
    print("drop deltaaa[**]",file=f2)
    print("store Dm3H",file=f2)

    print(file=f2)
    print(file=f2)
    print('---- code("MRCC_SBlock")',file=f2)
    print("// Set up overlap metric, ready to construct X used",file=f2)
    print("// in sequential orthogonalisation",file=f2)
    print("",file=f2)
    print("// I1",file=f2)
    print("alloc S1:I1[pq]",file=f2)
    print("load Dm1H[pq]",file=f2)
    print(".S1:I1[pq] := Dm1H[pq]",file=f2)
    print("drop Dm1H[pq]",file=f2)
    print("store S1:I1[pq]",file=f2)
    print("",file=f2)
    print("alloc S2:I1[pqrstu]",file=f2)
    print("load Dm3[pppppp], Dm2[pppp], Dm1[pp], deltaaa[pp]",file=f2)
    print(".S2:I1[pqrstu] := Dm3[pqrstu]",file=f2)
    print(".S2:I1[pqrstu] -= deltaaa[pt] Dm2[rqsu]",file=f2)
    print(".S2:I1[pqrstu] += deltaaa[pu] Dm2[rqst]",file=f2)
    print(".S2:I1[pqrstu] += deltaaa[qt] Dm2[rpsu]",file=f2)
    print(".S2:I1[pqrstu] -= deltaaa[qu] Dm2[rpst]",file=f2)
    print("alloc delta4[pppp]",file=f2)
    print(".delta4[pqtu] += deltaaa[pt] deltaaa[qu]",file=f2)
    print(".S2:I1[pqrstu] += delta4[pqtu] Dm1[rs]",file=f2)
    print(".S2:I1[pqrstu] -= delta4[qptu] Dm1[rs]",file=f2)
    print("drop delta4[pppp]",file=f2)
    print("drop deltaaa[pp], Dm1[pp], Dm2[pppp], Dm3[pppppp]",file=f2)
    print("store S2:I1[pqrstu]",file=f2)
    print("",file=f2)
    print("alloc S3:I1[pqrs]",file=f2)
    print("load Dm2[pppp], Dm1[pp], deltaaa[pp]",file=f2)
    print(".S3:I1[pqrs] -= Dm2[qprs]",file=f2)
    print(".S3:I1[pqrs] -= deltaaa[pr] Dm1[qs]",file=f2)
    print(".S3:I1[pqrs] += deltaaa[ps] Dm1[qr]",file=f2)
    print("drop deltaaa[pp], Dm1[pp], Dm2[pppp]",file=f2)
    print("store S3:I1[pqrs]",file=f2)
    print("",file=f2)
    print("// I2",file=f2)
    print("alloc S1:I2[pqrs]",file=f2)
    print("load Dm2H[pppp]",file=f2)
    print(".S1:I2[pqrs] := Dm2H[pqrs]",file=f2)
    print("drop Dm2H[pppp]",file=f2)
    print("store S1:I2[pqrs]",file=f2)
    print("",file=f2)
    print("// S0",file=f2)
    print("alloc S1:S0[pq]",file=f2)
    print("load Dm1[pp]",file=f2)
    print(".S1:S0[pq] := Dm1[pq]",file=f2)
    print("drop Dm1[pp]",file=f2)
    print("store S1:S0[pq]",file=f2)
    print("",file=f2)
    print("alloc S2:S0[pqrstu]",file=f2)
    print("load Dm3[pppppp], Dm2[pppp], deltaaa[pp]",file=f2)
    print(".S2:S0[pqrstu] -= Dm3[pqrstu]",file=f2)
    print(".S2:S0[pqrstu] += deltaaa[pu] Dm2[qrst]",file=f2)
    print("drop deltaaa[pp], Dm2[pppp], Dm3[pppppp]",file=f2)
    print("store S2:S0[pqrstu]",file=f2)
    print("",file=f2)
    print("alloc S3:S0[pqrs]",file=f2)
    print("load Dm2[pppp]",file=f2)
    print(".S3:S0[pqrs] := Dm2[pqsr]",file=f2)
    print("drop Dm2[pppp]",file=f2)
    print("store S3:S0[pqrs]",file=f2)
    print("",file=f2)
    print("// S1",file=f2)
    print("alloc S2:S1[pqrs]",file=f2)
    print("load Dm2[pppp], Dm1[pp], deltaaa[pp]",file=f2)
    print(".S2:S1[pqrs] -= Dm2[qprs]",file=f2)
    print(".S2:S1[pqrs] += deltaaa[ps] Dm1[qr]",file=f2)
    print("drop deltaaa[pp], Dm1[pp], Dm2[pppp]",file=f2)
    print("store S2:S1[pqrs]",file=f2)
    print("",file=f2)
    print("alloc S3:S1[pq]",file=f2)
    print("load Dm1[pp]",file=f2)
    print(".S3:S1[pq] := Dm1[qp]",file=f2)
    print("drop Dm1[pp]",file=f2)
    print("store S3:S1[pq]",file=f2)
    print("",file=f2)
    print("// S2",file=f2)
    print("alloc S1:S2[pq]",file=f2)
    print("load Dm1H[pq]",file=f2)
    print(".S1:S2[pq] := Dm1H[pq]",file=f2)
    print("drop Dm1H[pq]",file=f2)
    print("store S1:S2[pq]",file=f2)
    print("",file=f2)
    print("// P0",file=f2)
    print("alloc S1:P0[pqrs]",file=f2)
    print("load Dm2[pppp]",file=f2)
    print(".S1:P0[pqrs] := Dm2[pqrs]",file=f2)
    print("drop Dm2[pppp]",file=f2)
    print("store S1:P0[pqrs]",file=f2)
    print("",file=f2)
    print("// P1",file=f2)
    print("alloc S1:P1[pq]",file=f2)
    print("load Dm1[pp]",file=f2)
    print(".S1:P1[pq] := Dm1[pq]",file=f2)
    print("drop Dm1[pp]",file=f2)
    print("store S1:P1[pq]",file=f2)
    print(file=f2)

print("---- end", file=f2)

# Close off files (including temporary ones)
f2.close()
kext_temp.close()
init_res_temp.close()
