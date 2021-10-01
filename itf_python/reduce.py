#!/usr/bin/env python3

# =========================================================================================
# This file is used to post-proccess the ITF algo file and reduce the size of the file.
# So far, this removes repeated load/drop lines which invlove the same tensors
# =========================================================================================
import argparse     # Parse arguments

# =========================================================================================
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

# =========================================================================================
# =========================================================================================
#   main starts here:
# =========================================================================================
# =========================================================================================

# Parse arguments from gecco
parser = argparse.ArgumentParser(
                description="""Reduce the size of the ITF algo file""",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF algo file')
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)

inp = args.input
out = 'tmp.itfaa'

# Open bcontr.tmp file to read from
f1 = open(inp,"r")
# Open bcontr2.tmp file to read from
f2 = open(out,"w+")

# Each line will be saved a compared to the previous one
old_line = ""

# Read each line of bcontr.tmp
for line in f1:

    loop = False
    # Check if inside a loop (loops begin with 4x whitspace)
    if line.startswith(('   ')):
        loop = True

    # Look for consecutive load and drop line and see if they do the same thing
    if ('load' in line and 'drop' in old_line):

        words1 = line.split()
        words2 = old_line.split()

        if (len(words1) == len(words2)):
            # Compare tensor names (not index)
            tmp1=[]
            tmp2=[]
            for i in range(1, len(words1)):
                if not loop:
                    try:
                        tmp1.append(str(words1[i].split('[')[0])+str(generic_index(words1[i].split('[')[1])))
                        tmp2.append(str(words2[len(words2)-i].split('[',1)[0])+str(generic_index(words2[len(words2)-i].split('[')[1])))
                    except:
                        tmp1.append(words1[i])
                        tmp2.append(words2[len(words2)-i])
                else:
                    # Inside a loop, so should compare indices
                    tmp1.append(words1[i].split(']')[0])
                    tmp2.append(words2[len(words2)-i].split(']',1)[0])
            if (tmp1 == tmp2):
                #print("OLD2 ", old_line, end="", flush=True, file=f2)
                #print("OLD1 ", line, end="", flush=True, file=f2)
                #print("DELETE this", file=f2)
                old_line = ""
                line = ""
            else:
                print(old_line, end="", flush=True, file=f2)
                old_line = line
        else:
            print(old_line, end="", flush=True, file=f2)
            old_line = line
    else:
        print(old_line, end="", flush=True, file=f2)
        old_line = line


# Print the last reamining line of the file
print(old_line, end="", flush=True, file=f2)

# Clsoe the files
f1.close()
f2.close()
