def print_inter(prev_lines):
    # Load, contract, drop tensors involved with intermediates
    for i in range(0, len(prev_lines)):
        inter_words=prev_lines[i].split()
        if "STIN" not in inter_words[2] or "STIN" not in inter_words[3]:
            load_ten="load "
            if "STIN" not in inter_words[2]:
                load_ten=load_ten + inter_words[2].split('*',1)[-1]
            if "STIN" not in inter_words[2] and inter_words[3].split('*',1)[-1] != inter_words[2].split('*',1)[-1] and "STIN" not in inter_words[3]:
                load_ten=load_ten + ", "
            if "STIN" not in inter_words[3] and inter_words[3].split('*',1)[-1] != inter_words[2].split('*',1)[-1]:
                # Do not load if an intermediate or if the same as previous loaded tensor
                load_ten=load_ten + inter_words[3].split('*',1)[-1] 

            print(load_ten, file=out)

        print(prev_lines[i].strip(), file=out)

        if "STIN" not in inter_words[2] or "STIN" not in inter_words[3]:
            drop_ten="drop "
            if "STIN" not in inter_words[3]:
                drop_ten=drop_ten + inter_words[3].split('*',1)[-1]
            if "STIN" not in inter_words[3] and inter_words[2].split('*',1)[-1] != inter_words[3].split('*',1)[-1] and "STIN" not in inter_words[2]:
                drop_ten=drop_ten + ", "
            if "STIN" not in inter_words[2] and inter_words[2].split('*',1)[-1] != inter_words[3].split('*',1)[-1]:
                drop_ten=drop_ten + inter_words[2].split('*',1)[-1]

            print(drop_ten, file=out)

def print_result(line, words):
    # Load, contract, drop tensors involved with result tensors

    # Load tensors, cases depend on how many tensors are on the right
    if len(words)<=4:

        if "STIN" not in words[2] or "STIN" not in words[3]:
            load_ten="load "
            if "STIN" not in words[2]:
                load_ten=load_ten + words[2].split('*',1)[-1]
            if "STIN" not in words[2] and words[2].split('*',1)[-1] != words[3].split('*',1)[-1] and "STIN" not in words[3]:
                load_ten=load_ten + ", "
            if "STIN" not in words[3] and words[2].split('*',1)[-1] != words[3].split('*',1)[-1]:
                load_ten=load_ten + words[3].split('*',1)[-1]

            print(load_ten, file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "STIN" not in words[2] or "STIN" not in words[3]:
            drop_ten="drop "
            if "STIN" not in words[3]:
                drop_ten=drop_ten + words[3].split('*',1)[-1]
            if "STIN" not in words[3] and words[3].split('*',1)[-1] != words[2].split('*',1)[-1] and "STIN" not in words[2]:
                drop_ten=drop_ten + ", "
            if "STIN" not in words[2] and words[3].split('*',1)[-1] != words[2].split('*',1)[-1]:
                drop_ten=drop_ten + words[2].split('*',1)[-1]

            print(drop_ten, file=out)

    elif len(words)==6:

        # Line contains brackets, may have to load more than two tensors
        if "STIN" not in words:
            load_ten="load "
            if '(' in words[2]:
                t1=words[2].split('*',1)[-1].replace('(','')
                t2=words[4].split('*',1)[-1].replace(')','')
                t3=words[5].split('*',1)[-1]

                # Will need to compare generic index
                index2=list(words[2][words[2].find("[")+1:words[2].find("]")])
                index4=list(words[4][words[4].find("[")+1:words[4].find("]")])


                # Brackets for first tensor
                # Better way to do this???
                if "STIN" not in t1:
                    load_ten=load_ten + t1
                if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                    # Need to compare generic index as well
                    if generic_index(index2)!=generic_index(index4):
                        load_ten=load_ten + ", "
                if "STIN" not in t2 and t1 != t2:
                    if generic_index(index2)!=generic_index(index4):
                        load_ten=load_ten + t2
                if "STIN" not in t2 and t2 != t3 and \
                   "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                    load_ten=load_ten + ", "
                if "STIN" not in t3 and t2 != t3 and t1 != t3:
                    load_ten=load_ten + t3

            elif '(' in words[3]:
                t1=words[2].split('*',1)[-1]
                t2=words[3].split('*',1)[-1].replace('(','')
                t3=words[5].split('*',1)[-1].replace(')','')

                index3=list(words[3][words[3].find("[")+1:words[3].find("]")])
                index5=list(words[5][words[5].find("[")+1:words[5].find("]")])

                # Brackets for second tensor
                # Better way to do this???
                if "STIN" not in t1:
                    load_ten=load_ten + t1
                if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                    load_ten=load_ten + ", "
                if "STIN" not in t2 and t1 != t2:
                    load_ten=load_ten + t2
                #if "STIN" not in t2 and t2 != t3 and \
                #   "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                if "STIN" not in t2 and t2 != t3 and \
                   t1 != t3 and "STIN" not in t3:
                    # Need to compare generic index
                    if generic_index(index3)!=generic_index(index5):
                        load_ten=load_ten + ", "
                if "STIN" not in t3 and t2 != t3 and t1 != t3:
                    if generic_index(index3)!=generic_index(index5):
                        load_ten=load_ten + t3
            else:
                print("Error in bracket determination")
                exit(1)

            print(load_ten, file=out)

        print(line.strip(), file=out)

        # Drop tensors
        if "STIN" not in words:
            drop_ten="drop "
            if '(' in words[2]:
                t3=words[2].split('*',1)[-1].replace('(','')
                t2=words[4].split('*',1)[-1].replace(')','')
                t1=words[5].split('*',1)[-1]

                # Will need to compare generic index
                index2=list(words[2][words[2].find("[")+1:words[2].find("]")])
                index4=list(words[4][words[4].find("[")+1:words[4].find("]")])


                # Brackets for first tensor
                # Better way to do this???
                if "STIN" not in t1:
                    drop_ten=drop_ten + t1
                if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                    # Need to compare generic index as well
                    drop_ten=drop_ten + ", "
                if "STIN" not in t2 and t1 != t2:
                    drop_ten=drop_ten + t2
                if "STIN" not in t2 and t2 != t3 and "STIN" not in t3 \
                   and t1 != t3:
                    if generic_index(index2)!=generic_index(index4):
                        drop_ten=drop_ten + ", "
                if "STIN" not in t3 and t2 != t3 and t1 != t3:
                    if generic_index(index2)!=generic_index(index4):
                        drop_ten=drop_ten + t3

            elif '(' in words[3]:
                t3=words[2].split('*',1)[-1]
                t2=words[3].split('*',1)[-1].replace('(','')
                t1=words[5].split('*',1)[-1].replace(')','')

                index3=list(words[3][words[3].find("[")+1:words[3].find("]")])
                index5=list(words[5][words[5].find("[")+1:words[5].find("]")])

                # Brackets for second tensor
                if "STIN" not in t1:
                    drop_ten=drop_ten + t1
                if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                    # Need to compare generic index as well
                    if generic_index(index3)!=generic_index(index5):
                        drop_ten=drop_ten + ", "
                if "STIN" not in t2 and t1 != t2:
                    if generic_index(index3)!=generic_index(index5):
                        drop_ten=drop_ten + t2
                if "STIN" not in t2 and t2 != t3 and \
                   t1 != t3 and "STIN" not in t3:
                #if "STIN" not in t2 and t2 != t3 and \
                #   "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                    drop_ten=drop_ten + ", "
                if "STIN" not in t3 and t2 != t3 and t1 != t3:
                    drop_ten=drop_ten + t3
            else:
                print("Error in bracket determination")
                exit(1)

            print(drop_ten, file=out)

    elif len(words)>6:
        # Line contains two brackets, may have to load more than two tensors
        if "STIN" not in words:
            t1=words[2].split('*',1)[-1].replace('(','')
            t2=words[4].split('*',1)[-1].replace(')','')
            t3=words[5].split('*',1)[-1].replace('(','')
            t4=words[7].split('*',1)[-1].replace(')','')

            load_ten="load "

            # Will need to compare generic index
            index1=list(words[2][words[2].find("[")+1:words[2].find("]")])
            index2=list(words[4][words[4].find("[")+1:words[4].find("]")])
            index3=list(words[5][words[5].find("[")+1:words[5].find("]")])
            index4=list(words[7][words[7].find("[")+1:words[7].find("]")])

            if "STIN" not in t1:
                load_ten=load_ten + t1
            if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                # Need to compare generic index as well
                if generic_index(index1)!=generic_index(index2):
                    load_ten=load_ten + ", "
            if "STIN" not in t2 and t1 != t2:
                if generic_index(index1)!=generic_index(index2):
                    load_ten=load_ten + t2
            if "STIN" not in t2 and t2 != t3 and \
               "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                load_ten=load_ten + ", "
            if "STIN" not in t3 and t2 != t3 and t1 != t3:
                load_ten=load_ten + t3
            if "STIN" not in t1 and "STIN" not in t2 and "STIN" not in t3 and "STIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(index3)!=generic_index(index4):
                    load_ten=load_ten + ", "
            if "STIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(index3)!=generic_index(index4):
                    load_ten=load_ten + t4

            print(load_ten, file=out)

        print(line.strip(), file=out)


        # Drop tensors
        if "STIN" not in words:
            t4=words[2].split('*',1)[-1].replace('(','')
            t3=words[4].split('*',1)[-1].replace(')','')
            t2=words[5].split('*',1)[-1].replace('(','')
            t1=words[7].split('*',1)[-1].replace(')','')

            drop_ten="drop "

            # Will need to compare generic index
            index4=list(words[2][words[2].find("[")+1:words[2].find("]")])
            index3=list(words[4][words[4].find("[")+1:words[4].find("]")])
            index2=list(words[5][words[5].find("[")+1:words[5].find("]")])
            index1=list(words[7][words[7].find("[")+1:words[7].find("]")])

            if "STIN" not in t1:
                drop_ten=drop_ten + t1
            if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                # Need to compare generic index as well
                if generic_index(index1)!=generic_index(index2):
                    drop_ten=drop_ten + ", "
            if "STIN" not in t2 and t1 != t2:
                if generic_index(index1)!=generic_index(index2):
                    drop_ten=drop_ten + t2
            if "STIN" not in t2 and t2 != t3 and \
               "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                drop_ten=drop_ten + ", "
            if "STIN" not in t3 and t2 != t3 and t1 != t3:
                drop_ten=drop_ten + t3
            if "STIN" not in t1 and "STIN" not in t2 and "STIN" not in t3 and "STIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(index3)!=generic_index(index4):
                    drop_ten=drop_ten + ", "
            if "STIN" not in t4 and \
               t1 != t4 and t2 != t4 and t3 != t4:
                if generic_index(index3)!=generic_index(index4):
                    drop_ten=drop_ten + t4

            print(drop_ten, file=out)



def add_to_global(word,declare_ten,declare_ten_index,declare_ten_name):
    if word.split('*',1)[-1] not in declare_ten and "STIN" not in word:
        index=list(word[word.find("[")+1:word.find("]")])

        # Construct generic index
        generic=generic_index(index)

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


def generic_index(index):
    # Construct generic index representation of specific tensor index
    # Used to check if tensors index the same space, even if they have
    # different indices

    hole=['i','j','k','l']
    particle=['a','b','c','d']
    valence=['p','q','r','s','t','u','v','w']

    gen=[]
    for i in range (0,len(index)):
        if index[i] in particle:
            gen.append('e')
        if index[i] in valence:
            gen.append('a')
        if index[i] in hole:
            gen.append('c')
    return gen




import argparse

parser = argparse.ArgumentParser(
                description="""Process ITF binary contraction file and output ITF algo file
                """,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF binary contraction file')
parser.add_argument('-o','--output',default=None,help='ITF algo file')
# s = 0: Single-reference, no need for overlap tensors in algo file
#     1: Multireference, print out all overlap tensors
parser.add_argument('-s','--overlap',default=1,help='Overlap tensors needed for a specific mrcc calculation')
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
outp = args.output
olap = args.overlap

f=open(inp,"r")
out=open(outp, "w+")

prev_lines=[]       # Previous intermediate lines which belong to next result block
prev_inter=[]       # List of previous intemediates, used to alloc/drop
prev_res='#####'    # Result of previous line
prev_generic=[]     # Previous generic result

declare_res=[]      # Global list of result tensors
declare_name=[]     # Global of just result tensor names
declare_index=[]    # Global list of result indices

declare_inter=[]        # Global list of intermediates
declare_inter_index=[]  # Global list of intermediates
declare_inter_name=[]   # Global list of intermediates

declare_ten=[]
declare_ten_index=[]
declare_ten_name=[]

# Spin summed family = One or more equations the arise from the spin summation
# proccedure on one result tensor contraction line
begin=False         # Marks the start of a spin summed family of contractions
end=False           # Marks the end of a spin summed family of contractions
old_spin_iter=[]    # Stores list of intermediates used throughout the spin summed family

for line in f:
    words=line.split()

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

    # Check if brackets in the binary contraction
    if '(' in words[2] and '(' in words[5]:
        # Check if first tensor needs to be declared
        add_to_global(words[2].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[4].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
        # Check if second tensor needs to be declared
        add_to_global(words[5].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[7].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)

        # Arrange tensors so they are in positions 2 and 3
#        words[2]=words[2].replace('(','')
#        words[3]=words[5].replace('(','')
    elif '(' in words[3]:
        # Check if first tensor needs to be declared
        add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[3].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[5].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)

        # Brackets for second tensor
#        words[3]=words[3].replace('(','')
    elif '(' in words[2] :
        # Check if first tensor needs to be declared
        add_to_global(words[2].replace('(',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[4].replace(')',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[5],declare_ten,declare_ten_index,declare_ten_name)

        # Brackets for first tensor
#        words[2]=words[2].replace('(','')
#        words[3]=words[5]
    else:
        # No brackets
        # Check if first tensor needs to be declared
        add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)
        # Check if second tensor needs to be declared
        add_to_global(words[3],declare_ten,declare_ten_index,declare_ten_name)


    if "STIN" in words[0]:
        # Check if contraction forms and intermediate
        prev_lines.append(line)
        prev_inter.append(words[0].replace('.',''))
        if words[0].replace('.','') not in declare_inter:
            # Add intermedate to global list
            index=list(words[0][words[0].find("[")+1:words[0].find("]")])
    
            generic=generic_index(index)
    
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
        index=list(words[0][words[0].find("[")+1:words[0].find("]")])

        generic=generic_index(index)

        if words[0] != prev_res and generic != prev_generic:
            # Next result is different from previous, so close off block

            if prev_res != '#####':
                print("store", prev_res.replace('.',''), file=out)
                print(file=out)

            # Check whether to load previously allocated tensor, or load it back from memory
            # To load, the tensor name and generic index associated with it must be equal
            loaded=False
            for i in range(0, len(declare_name)):
                if words[0].split('[',1)[0].replace('.','') == declare_name[i]:
                    if generic == declare_index[i]:
                        # Generic index must be at same position as name it belongs to - dangerous! 
                        # Load previous tensor
                        print("load", words[0].replace('.',''), file=out)
                        loaded=True
                        break
                else:
                    continue

            if not loaded:
                # Add result to global list and alloc new result
                #declare_res.append(words[0].replace('.',''))
                declare_res.append(words[0].split('.',1)[1])
                declare_index.append(generic)
                declare_name.append(words[0].split('[',1)[0].replace('.',''))
                print("alloc ", words[0].replace('.',''), file=out)
                

            # Alloc intermediates if needed
            if prev_inter:
                print("alloc ", end="", flush=True, file=out)
                print(*prev_inter,sep=", ", file=out)

            # Print intermediates and load/drop relavant tensors
            print_inter(prev_lines)

            # Print result line
            print_result(line, words)

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
            print_result(line, words)

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
# Declare tensros and index-spaces
print("---- decl", file=f2)
print("index-space: pqrstuvw , Active    , a", file=f2)
print("index-space: ijkl     , Closed    , c", file=f2)
print("index-space: gh       , ClosedF   , f", file=f2)
print("index-space: abcd     , External  , e", file=f2)
print("index-space: mno      , Internal  , i", file=f2)
print(file=f2)

# Print out tensors used in contractions
for i in range(0, len(declare_ten)):
    print("tensor: ", declare_ten[i], file=f2)

# Print of result tensors
print(file=f2)
for i in range(0, len(declare_res)):
    print("tensor: ", declare_res[i], file=f2)

# Declare density and overlap tensors
if (olap>0):
    print(file=f2)
    print("// Resuced density tensors (Icc-Icc coupling-coefficients)",file=f2)
    print("// these are created on C++ side in CreateMrciTensors method",file=f2)
    print("tensor: Dm1[pp]       , DDm1",file=f2)
    print("tensor: Dm2[pppp]     , DDm2",file=f2)
    print("tensor: Dm3[pppppp]   , DDm3",file=f2)
    print("tensor: Dm1H[pp]      , DHm1",file=f2)
    print("tensor: Dm2H[pppp]    , DHm2",file=f2)
    print("tensor: Dm3H[pppppp]  , DHm3",file=f2)
    print(file=f2)

    print("// Delta tensors",file=f2)
    print("tensor: deltaaa[pq]   , DeltaActAct",file=f2)
    print("tensor: delta4[pppp]  , !Create{type:plain}, Delta4  // Intermediate rank4 delta tensor",file=f2)
    print(file=f2)

    print("// Overlap tensors labeled by the exciation class",file=f2)
    print("tensor: S1:I1[pp]     , S1:I1",file=f2)
    print("tensor: S2:I1[pppppp] , S2:I1",file=f2)
    print("tensor: S3:I1[pppp]   , S3:I1",file=f2)
    print("tensor: S1:I2[pppp]   , S1:I2",file=f2)
    print("tensor: S1:S0[pp]     , S1:S0",file=f2)
    print("tensor: S2:S0[pppppp] , S2:S0",file=f2)
    print("tensor: S3:S0[pppp]   , S3:S0",file=f2)
    print("tensor: S2:S1[pppp]   , S2:S1",file=f2)
    print("tensor: S3:S1[pp]     , S3:S1",file=f2)
    print("tensor: S1:S2[pp]     , S1:S2",file=f2)
    print("tensor: S1:P0[pppp]   , S1:P0",file=f2)
    print("tensor: S1:P1[pp]     , S1:P1",file=f2)
    print(file=f2)

# Print out intermediates
print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    print("tensor: %-20s, !Create{type:disk}" % (declare_inter[i]), file=f2)

# Print out code blocks
print(file=f2)
print(file=f2)
print('---- code("MRCC_Residual")', file=f2)
f2.write(tmp)
print(file=f2)

# Symmetrise tensors
if "O2g[abij]" in declare_res:
    print("load O2g[abij]", file=f2)
    print(".O2g[abij] += O2g[baji]", file=f2)
    print("store O2g[abij]", file=f2)
    print(file=f2)

if "O2g[pqij]" in declare_res:
    print("load O2g[pqij]", file=f2)
    print(".O2g[pqij] += O2g[qpji]", file=f2)
    print("store O2g[pqij]", file=f2)
    print(file=f2)

if "O2g[abpq]" in declare_res:
    print("load O2g[abpq]", file=f2)
    print(".O2g[abpq] += O2g[baqp]", file=f2)
    print("store O2g[abpq]", file=f2)
    print(file=f2)

# Print out code needed to evaluate the overlap matrix
if (olap>0):
    print(file=f2)
    print("// Set up 3rd order denisty and hole tensors")
    print("// This is taken from the cic code")
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
    print('---- code("icMRCC_SBlock")',file=f2)
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
