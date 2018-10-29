class itf_line:
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
        names = ["R", "T", "K", "f"]
        for i in range(0, len(names)):
            self.rename_line_names(names[i])
        self.line = " ".join(self.parts)

    def rename_line_names(self, name):
        """Add generic index to tensor name"""
        if (name+"[" in self.line):
            for i in range(0, len(self.parts)):
                if (name+"[" in self.parts[i]):
                    self.parts[i] = self.parts[i].replace(name, name + ":" + "".join(generic_index(self.parts[i])))


def print_inter(prev_lines):
    # Load, contract, drop tensors involved with intermediates
    for i in range(0, len(prev_lines)):
        inter_words=prev_lines[i].split()

        if len(inter_words)==3:
            if "TIN" not in inter_words[2]:
                load_ten="load " + inter_words[2].split('*',1)[-1]

                print(load_ten.strip(), file=out)

            print(prev_lines[i].strip(), file=out)

            if "TIN" not in inter_words[2]:
                drop_ten="drop " + inter_words[2].split('*',1)[-1]

                print(drop_ten.strip(), file=out)

        else:
            if "TIN" not in inter_words[2] or "TIN" not in inter_words[3]:
                load_ten="load "
                if "TIN" not in inter_words[2]:
                    load_ten=load_ten + inter_words[2].split('*',1)[-1]
                if "TIN" not in inter_words[2] and inter_words[3].split('*',1)[-1] != inter_words[2].split('*',1)[-1] and "TIN" not in inter_words[3]:
                    load_ten=load_ten + ", "
                if "TIN" not in inter_words[3] and inter_words[3].split('*',1)[-1] != inter_words[2].split('*',1)[-1]:
                    # Do not load if an intermediate or if the same as previous loaded tensor
                    load_ten=load_ten + inter_words[3].split('*',1)[-1]

                print(load_ten.strip(), file=out)

            print(prev_lines[i].strip(), file=out)

            if "TIN" not in inter_words[2] or "TIN" not in inter_words[3]:
                drop_ten="drop "
                if "TIN" not in inter_words[3]:
                    drop_ten=drop_ten + inter_words[3].split('*',1)[-1]
                if "TIN" not in inter_words[3] and inter_words[2].split('*',1)[-1] != inter_words[3].split('*',1)[-1] and "TIN" not in inter_words[2]:
                    drop_ten=drop_ten + ", "
                if "TIN" not in inter_words[2] and inter_words[2].split('*',1)[-1] != inter_words[3].split('*',1)[-1]:
                    drop_ten=drop_ten + inter_words[2].split('*',1)[-1]

                print(drop_ten.strip(), file=out)


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
    change_line_names("T", line_o, words)
    change_line_names("K", line_o, words)
    change_line_names("f", line_o, words)
    line = " ".join(words)
    return line


def print_result(line):
    # Load, contract, drop tensors involved with result tensors

    # Change tensor names within the line
    #line = change_line(line_o) 
    words=line.split()

    # Load tensors, cases depend on how many tensors are on the right
    if len(words)==3:

        # This line adds the reference energy to the correlation energy
        # So we can skip it
        if ("ECC[]" in words[0] and "K:[]" in words[2]): return

        # Either a simple adding or copying case
        if "TIN" not in words[2]:
            load_ten="load " + words[2].split('*',1)[-1]

            print(load_ten.strip(), file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "TIN" not in words[2]:
            drop_ten="drop " +  words[2].split('*',1)[-1]

            print(drop_ten.strip(), file=out)

    elif len(words)==4:

        if "TIN" not in words[2] or "TIN" not in words[3]:
            load_ten="load "
            if "TIN" not in words[2]:
                load_ten=load_ten + words[2].split('*',1)[-1]
            if "TIN" not in words[2] and words[2].split('*',1)[-1] != words[3].split('*',1)[-1] and "TIN" not in words[3]:
                load_ten=load_ten + ", "
            if "TIN" not in words[3] and words[2].split('*',1)[-1] != words[3].split('*',1)[-1]:
                load_ten=load_ten + words[3].split('*',1)[-1]

            print(load_ten.strip(), file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "TIN" not in words[2] or "TIN" not in words[3]:
            drop_ten="drop "
            if "TIN" not in words[3]:
                drop_ten=drop_ten + words[3].split('*',1)[-1]
            if "TIN" not in words[3] and words[3].split('*',1)[-1] != words[2].split('*',1)[-1] and "TIN" not in words[2]:
                drop_ten=drop_ten + ", "
            if "TIN" not in words[2] and words[3].split('*',1)[-1] != words[2].split('*',1)[-1]:
                drop_ten=drop_ten + words[2].split('*',1)[-1]

            print(drop_ten.strip(), file=out)

    elif len(words)==6:

        # Line contains brackets, may have to load more than two tensors
        if "TIN" not in words:
            load_ten="load "
            if '(' in words[2]:
                t1=words[2].split('*',1)[-1].replace('(','')
                t2=words[4].split('*',1)[-1].replace(')','')
                t3=words[5].split('*',1)[-1]


                # Brackets for first tensor
                # Better way to do this???
                if "TIN" not in t1:
                    load_ten=load_ten + t1
                if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                    # Need to compare generic index as well
                    if generic_index(words[2])!=generic_index(words[4]):
                        load_ten=load_ten + ", "
                if "TIN" not in t2 and t1 != t2:
                    if generic_index(words[2])!=generic_index(words[4]):
                        load_ten=load_ten + t2
                if "TIN" not in t2 and t2 != t3 and \
                   "TIN" not in t1 and t1 != t3 and "TIN" not in t3:
                    load_ten=load_ten + ", "
                if "TIN" not in t3 and t2 != t3 and t1 != t3:
                    load_ten=load_ten + t3

            elif '(' in words[3]:
                t1=words[2].split('*',1)[-1]
                t2=words[3].split('*',1)[-1].replace('(','')
                t3=words[5].split('*',1)[-1].replace(')','')

                # Brackets for second tensor
                # Better way to do this???
                if "TIN" not in t1:
                    load_ten=load_ten + t1
                if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                    load_ten=load_ten + ", "
                if "TIN" not in t2 and t1 != t2:
                    load_ten=load_ten + t2
                #if "TIN" not in t2 and t2 != t3 and \
                #   "TIN" not in t1 and t1 != t3 and "TIN" not in t3:
                if "TIN" not in t2 and t2 != t3 and \
                   t1 != t3 and "TIN" not in t3:
                    # Need to compare generic index
                    if generic_index(words[3])!=generic_index(words[5]):
                        load_ten=load_ten + ", "
                if "TIN" not in t3 and t2 != t3 and t1 != t3:
                    if generic_index(words[3])!=generic_index(words[5]):
                        load_ten=load_ten + t3
            else:
                print("Error in bracket determination")
                exit(1)

            print(load_ten.strip(), file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "TIN" not in words:
            drop_ten="drop "
            if '(' in words[2]:
                t3=words[2].split('*',1)[-1].replace('(','')
                t2=words[4].split('*',1)[-1].replace(')','')
                t1=words[5].split('*',1)[-1]

                # Brackets for first tensor
                # Better way to do this???
                if "TIN" not in t1:
                    drop_ten=drop_ten + t1
                if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                    # Need to compare generic index as well
                    drop_ten=drop_ten + ", "
                if "TIN" not in t2 and t1 != t2:
                    drop_ten=drop_ten + t2
                if "TIN" not in t2 and t2 != t3 and "TIN" not in t3 \
                   and t1 != t3:
                    if generic_index(words[2])!=generic_index(words[4]):
                        drop_ten=drop_ten + ", "
                if "TIN" not in t3 and t2 != t3 and t1 != t3:
                    if generic_index(words[2])!=generic_index(words[4]):
                        drop_ten=drop_ten + t3

            elif '(' in words[3]:
                t3=words[2].split('*',1)[-1]
                t2=words[3].split('*',1)[-1].replace('(','')
                t1=words[5].split('*',1)[-1].replace(')','')

                # Brackets for second tensor
                if "TIN" not in t1:
                    drop_ten=drop_ten + t1
                if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                    # Need to compare generic index as well
                    if generic_index(words[3])!=generic_index(words[5]):
                        drop_ten=drop_ten + ", "
                if "TIN" not in t2 and t1 != t2:
                    if generic_index(words[3])!=generic_index(words[5]):
                        drop_ten=drop_ten + t2
                if "TIN" not in t2 and t2 != t3 and \
                   t1 != t3 and "TIN" not in t3:
                #if "TIN" not in t2 and t2 != t3 and \
                #   "TIN" not in t1 and t1 != t3 and "TIN" not in t3:
                    drop_ten=drop_ten + ", "
                if "TIN" not in t3 and t2 != t3 and t1 != t3:
                    drop_ten=drop_ten + t3
            else:
                print("Error in bracket determination")
                exit(1)

            print(drop_ten.strip(), file=out)

    elif len(words)>6:
        # Line contains two brackets, may have to load more than two tensors
        if "TIN" not in words:
            t1=words[2].split('*',1)[-1].replace('(','')
            t2=words[4].split('*',1)[-1].replace(')','')
            t3=words[5].split('*',1)[-1].replace('(','')
            t4=words[7].split('*',1)[-1].replace(')','')

            load_ten="load "

            if "TIN" not in t1:
                load_ten=load_ten + t1
            if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                # Need to compare generic index as well
                if generic_index(words[2])!=generic_index(words[4]):
                    load_ten=load_ten + ", "
            if "TIN" not in t2 and t1 != t2:
                if generic_index(words[2])!=generic_index(words[4]):
                    load_ten=load_ten + t2
            if "TIN" not in t2 and t2 != t3 and \
               "TIN" not in t1 and t1 != t3 and "TIN" not in t3:
                load_ten=load_ten + ", "
            if "TIN" not in t3 and t2 != t3 and t1 != t3:
                load_ten=load_ten + t3
            if "TIN" not in t1 and "TIN" not in t2 and "TIN" not in t3 and "TIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(words[5])!=generic_index(words[7]):
                    load_ten=load_ten + ", "
            if "TIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(words[5])!=generic_index(words[7]):
                    load_ten=load_ten + t4

            print(load_ten.strip(), file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "TIN" not in words:
            t4=words[2].split('*',1)[-1].replace('(','')
            t3=words[4].split('*',1)[-1].replace(')','')
            t2=words[5].split('*',1)[-1].replace('(','')
            t1=words[7].split('*',1)[-1].replace(')','')

            drop_ten="drop "

            if "TIN" not in t1:
                drop_ten=drop_ten + t1
            if "TIN" not in t1 and t1 != t2 and "TIN" not in t2:
                # Need to compare generic index as well
                if generic_index(words[7])!=generic_index(words[5]):
                    drop_ten=drop_ten + ", "
            if "TIN" not in t2 and t1 != t2:
                if generic_index(words[7])!=generic_index(words[5]):
                    drop_ten=drop_ten + t2
            if "TIN" not in t2 and t2 != t3 and \
               "TIN" not in t1 and t1 != t3 and "TIN" not in t3:
                drop_ten=drop_ten + ", "
            if "TIN" not in t3 and t2 != t3 and t1 != t3:
                drop_ten=drop_ten + t3
            if "TIN" not in t1 and "TIN" not in t2 and "TIN" not in t3 and "TIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(words[4])!=generic_index(words[2]):
                    drop_ten=drop_ten + ", "
            if "TIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(words[4])!=generic_index(words[2]):
                    drop_ten=drop_ten + t4

            print(drop_ten.strip(), file=out)



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


def generic_index(tensor):
    # Construct generic index representation of specific tensor index
    # Used to check if tensors index the same space, even if they have
    # different indices
        
    index=list(tensor[tensor.find("[")+1:tensor.find("]")])

    hole=['i','j','k','l']
    particle=['a','b','c','d']
    valence=['p','q','r','s','t','u','v','w']

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


# =========================================================================================
# Main program starts here
# =========================================================================================
import argparse     # Parse arguments
import datetime     # Get tima and date
import time         # For timings

# Parse arguments from gecco
parser = argparse.ArgumentParser(
                description="""Process ITF binary contraction file and output ITF algo file
                """,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF binary contraction file')
parser.add_argument('-o','--output',default=None,help='ITF algo file')
parser.add_argument('-s','--overlap',default=1,help='0: single-reference; >0: multireference')
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
outp = args.output
olap = int(args.overlap)

# Open bcontr.tmp file to read from
f=open(inp,"r")
# Open first output file
out=open(outp, "w+")

# Declare lists needed in program
prev_lines=[]       # Previous intermediate lines which belong to next result block
prev_inter=[]       # List of previous intemediates, used to alloc/drop
prev_res='#####'    # Result of previous line
prev_generic=[]     # Previous generic result

declare_res=[]      # Global list of result tensors
declare_index=[]    # Global list of result indices
declare_name=[]     # Global of just result tensor names

declare_inter=[]        # Global list of intermediates
declare_inter_index=[]  # Global list of intermediates
declare_inter_name=[]   # Global list of intermediates

declare_ten=[]          # Global list of tensors involved in binary contractions
declare_ten_index=[]    # Global list of tensor indicies
declare_ten_name=[]     # Global list of tensor names

prev_K4E_lines={0:"start"}
K4E_count=1

# Spin summed family = One or more equations the arise from the spin summation
# proccedure on one result tensor contraction line
begin=False         # Marks the start of a spin summed family of contractions
end=False           # Marks the end of a spin summed family of contractions
old_spin_iter=[]    # Stores list of intermediates used throughout the spin summed family

# Read each line of bcontr.tmp and process it
for line_o in f:

    tensor_line = itf_line(line_o, out)
    tensor_line.rename_line()

    # Change names of external tensors (add : + generic index to name)
    line = change_line(line_o)
    words=line.split()

    # Check for spin summed block
    if (words[0]=='BEGIN'):
        # BEGIN marks the begining of a spin summed series of contractions

        begin=True
        end=False
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


    # Check if brackets in the binary contraction
    if (len(words)>=4):
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
            for i in range(0, len(prev_inter)):
                if words[0].split('[',1)[0].replace('.','') not in prev_inter[i].split('[',1)[0].replace('.',''):
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

                print("store", prev_res.replace('.',''), file=out)
                print(file=out)

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
                print("alloc", tmp_res.replace('.',''), file=out)
                

            # Alloc intermediates if needed
            if prev_inter:
                print("alloc ", end="", flush=True, file=out)
                print(*prev_inter,sep=", ", file=out)

            # Print intermediates and load/drop relavant tensors
            print_inter(prev_lines)

            # Print result line
            print_result(line)

            # Drop intermediates if needed
            if prev_inter:
                print("drop ", end="", flush=True, file=out)
                print(*list(reversed(prev_inter)),sep=", ", file=out)

            prev_res=words[0]
            prev_lines=[]
            prev_inter=[]

        #elif words[0] == prev_res:
        else:
            # Still within the same result block

            # Alloc intermediates if needed
            if prev_inter:
                print("alloc ", end="", flush=True, file=out)
                print(*prev_inter,sep=", ", file=out)

            # Print intermediates and load/drop relavant tensors
            print_inter(prev_lines)

            # Print result line
            print_result(line)

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
out.close()
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
if (olap==0):
    print("index-space: ijkl, Closed  , c", file=f2)
    print("index-space: abcd, External, e", file=f2)
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
        or "Dm" in declare_ten[i]): continue
    if ("[]" in declare_ten[i]):
        print("tensor:", declare_ten[i] + ",  !Create{type:scalar}", file=f2)
    else:
        print("tensor:", declare_ten[i] + ",  !Create{type:disk}", file=f2)

# Print already existing tensor, ie. don't need !Create{type:disk}
declare_existing_tensors(declare_ten, "Integral tensors", "K")
declare_existing_tensors(declare_ten, "Special integral tensors", "K4E",True)
declare_existing_tensors(declare_ten, "Fock tensors", "f")
declare_existing_tensors(declare_ten, "Amplitude tensors", "T")
if (olap): print("tensor: R[I],  R:I", file=f2)
declare_existing_tensors(declare_res, "Residual tensors", "R")
declare_existing_tensors(declare_res, "Energy and DIIS scalars", "ECC", True)

if (olap==0):
    # Tensors needed in CCD
    print("tensor: ERef[], ERef     // Reference energy", file=f2)
    print("tensor: EMp2[], EMp2     // MP2 energy", file=f2)
    print("tensor: EDi2[], EDi2     // Direct 2nd order energy", file=f2)
    print("tensor: Nrm2[], Nrm2     // Doubles amplitude norm", file=f2)
    print("tensor: Var2[], Var2     // Doubles residual norm", file=f2)
    print(file=f2)
    print("// Tensors needed to calculate the reference energy", file=f2)
    print("tensor: f:CC[CC],   f:CC", file=f2)
    print("tensor: CoreH[ii],  h:cc", file=f2)
    print("tensor: CoreH[CC],  h:CC", file=f2)
    print("tensor: Delta[ii],  Delta", file=f2)
    print("tensor: DeltaC[CC], DeltaC", file=f2)

# Declare density and overlap tensors
if (olap>0):
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

# Print out intermediates
print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    if ("[]" in declare_inter[i]):
        print("tensor: %-18s !Create{type:scalar}" % (declare_inter[i] + ","), file=f2)
    else:
        print("tensor: %-18s !Create{type:plain}" % (declare_inter[i] + ","), file=f2)
if (olap==0):
    print(file=f2)
    print("tensor: L1[abij],          !Create{type:plain}", file=f2)

# Print out code blocks
# Need to initalise the amplitudes first
print(file=f2)
print(file=f2)
print('---- code("Init_Amplitudes")',file=f2)
if (olap>0):
    for i in range(0, len(declare_ten)):
        if ("T[" in declare_ten[i]):
            generic=generic_index(declare_ten[i])
            declare_ten[i] = declare_ten[i][:1] + ":" + "".join(generic) + declare_ten[i][1:]

            print("alloc", declare_ten[i], file=f2)
            print("store", declare_ten[i], file=f2)
else:
    # Initalise amplitudes using MP2
    print("// Using MP2 amplitudes for starting guess", file=f2)
    print("alloc EMp2[], Nrm2[]", file=f2)
    print("for [i,j]:", file=f2)
    print("   alloc T:eecc[abij]", file=f2)
    print("   load K:eecc[**ij]", file=f2)
    print("   .T:eecc[abij] -= K:eecc[abij]", file=f2)
    print("   denom-scale T:eecc[abij], [1,1,0,0]", file=f2)
    print("   .EMp2[] += (2.0*T:eecc[abij] - T:eecc[baij]) K:eecc[abij]", file=f2)
    print("   .Nrm2 += (2.0*T:eecc[abij] - T:eecc[baij]) T:eecc[abij]", file=f2)
    print("   drop K:eecc[**ij]", file=f2)
    print("   store T:eecc[**ij]", file=f2)
    print("store Nrm2[], EMp2[]", file=f2)

# Calculate the reference energy for single-reference methods
if (olap==0):
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

# Print out amplitude update
# For single-reference methods, this also calculates the energy
if (olap==0):
    # Update for CCD
    print(file=f2)
    print(file=f2)
    print('---- code("Update_Amplitudes")',file=f2)
    print("// Update doubles",file=f2)
    print("alloc ECC[], EDi2[], Nrm2[], Var2[]",file=f2)
    print("for [i,j]:",file=f2)
    print("   load R:eecc[**ij]",file=f2)
    print("   // L1 = R^{ij}_{ab}/D^{ij}_{ab}; D^{ij}_{ab} = e_a+e_b-e_i-e_j",file=f2)
    print("   alloc L1[**ij]",file=f2)
    print("   .L1[**ij] += R:eecc[**ij]",file=f2)
    print("   denom-scale L1[**ij], [1,1,0,0]",file=f2)
    print("",file=f2)
    print("   load K:eecc[**ij]",file=f2)
    print("   // Compute Hylleraas-functional-like energy",file=f2)
    print("   // \tilde{T}^{ij}_{ab} K{ij}_{ab} + \tilde{(T^{ij}_{ab}-R^{ij}_{ab}/D^{ij}_{ab})} R{ij}_{ab}",file=f2)
    print("   load T:eecc[**ij]",file=f2)
    print("   .ECC += (2.0*T:eecc[abij] - T:eecc[baij]) (R:eecc[abij] + K:eecc[abij])",file=f2)
    print("   drop T:eecc[**ij]",file=f2)
    print("   .ECC -= (2.0*L1[abij] - L1[baij]) R:eecc[abij]",file=f2)
    print("   .Var2 += (2.0*L1[abij] - L1[baij]) L1[abij]",file=f2)
    print("",file=f2)
    print("   // Update T:eecc",file=f2)
    print("   // T^{ij}_{ab} = T^{ij}_{ab}-R^{ij}_{ab}/D^{ij}_{ab}",file=f2)
    print("   load T:eecc[**ij]",file=f2)
    print("   .T:eecc[abij] -= L1[abij]",file=f2)
    print("",file=f2)
    print("   .EDi2 += (2.0*T:eecc[abij] - T:eecc[baij]) K:eecc[abij]",file=f2)
    print("   .Nrm2 += (2.0*T:eecc[abij] - T:eecc[baij]) T:eecc[abij]",file=f2)
    print("   store T:eecc[**ij]",file=f2)
    print("   drop K:eecc[**ij]",file=f2)
    print("",file=f2)
    print("   drop L1[**ij]",file=f2)
    print("   drop R:eecc[**ij]",file=f2)
    print("store Var2[], Nrm2[], EDi2[], ECC[]",file=f2)

# Print out residual equations
print(file=f2)
print(file=f2)
print('---- code("Residual")', file=f2)
f2.write(tmp)
print(file=f2)

# Symmetrise tensors
if "R[abij]" in declare_res:
    print("load R:eecc[abij]", file=f2)
    print(".R:eecc[abij] += R:eecc[baji]", file=f2)
    print("store R:eecc[abij]", file=f2)
    print(file=f2)

if "R[pqij]" in declare_res:
    print("load R:aacc[pqij]", file=f2)
    print(".R:aacc[pqij] += R:aacc[qpji]", file=f2)
    print("store R:aacc[pqij]", file=f2)
    print(file=f2)

if "R[abpq]" in declare_res:
    print("load R:eeaa[abpq]", file=f2)
    print(".R:eeaa[abpq] += R:eeaa[baqp]", file=f2)
    print("store R:eeaa[abpq]", file=f2)
    print(file=f2)

# Print out code needed to evaluate the overlap matrix
if (olap>0):
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

f2.close()
