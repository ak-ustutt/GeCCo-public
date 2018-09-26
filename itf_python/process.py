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
                #t3=words[2].split('*',1)[-1]
                #t2=words[3].split('*',1)[-1].replace('(','')
                #t1=words[5].split('*',1)[-1].replace(')','')

                #index3=list(words[3][words[3].find("[")+1:words[3].find("]")])
                #index5=list(words[5][words[5].find("[")+1:words[5].find("]")])

                ## Brackets for second tensor
                ## Better way to do this???
                #if "STIN" not in t1:
                #    drop_ten=drop_ten + t1
                #if "STIN" not in t1 and t1 != t2 and "STIN" not in t2:
                #    drop_ten=drop_ten + ", "
                #if "STIN" not in t2 and t1 != t2:
                #    drop_ten=drop_ten + t2
                ##if "STIN" not in t2 and t2 != t3 and \
                ##   "STIN" not in t1 and t1 != t3 and "STIN" not in t3:
                #if "STIN" not in t2 and t2 != t3 and \
                #   t1 != t3 and "STIN" not in t3:
                #    # Need to compare generic index
                #    if generic_index(index3)!=generic_index(index5):
                #        drop_ten=drop_ten + ", "
                #if "STIN" not in t3 and t2 != t3 and t1 != t3:
                #    if generic_index(index3)!=generic_index(index5):
                #        drop_ten=drop_ten + t3

                index3=list(words[3][words[3].find("[")+1:words[3].find("]")])
                index5=list(words[5][words[5].find("[")+1:words[5].find("]")])

                # Brackets for second tensor
                # Better way to do this???
                if "STIN" not in words[5]:
                    drop_ten=drop_ten + words[5].split('*',1)[-1].replace(')','')
                if "STIN" not in words[5] and words[5].split('*',1)[-1].replace(')','') != words[3].split('*',1)[-1].replace('(','') and "STIN" not in words[3]:
                    # Need to compare generic index as well
                    if generic_index(index3)!=generic_index(index5):
                        drop_ten=drop_ten + ", "
                if "STIN" not in words[3] and words[5].split('*',1)[-1].replace(')','') != words[3].split('*',1)[-1].replace('(',''):
                    if generic_index(index3)!=generic_index(index5):
                        drop_ten=drop_ten + words[3].split('*',1)[-1].replace('(','')
                if "STIN" not in words[3] and words[3].split('*',1)[-1].replace('(','') != words[2].split('*',1)[-1] and \
                   "STIN" not in words[5] and words[5].split('*',1)[-1].replace(')','') != words[2].split('*',1)[-1] and "STIN" not in words[5]:
                    drop_ten=drop_ten + ", "
                if "STIN" not in words[2] and words[3].split('*',1)[-1].replace('(','') != words[2].split('*',1)[-1] and \
                    words[5].split('*',1)[-1].replace(')','') != words[2].split('*',1)[-1]:
                    drop_ten=drop_ten + words[2].split('*',1)[-1]
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
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
outp = args.output

f=open(inp,"r")
out=open(outp, "w+")

prev_lines=[]       # Previous intermediate lines which belong to next result block
prev_inter=[]       # List of previous intemediates, used to alloc/drop
prev_res='#####'    # Result of previous line
prev_generic=[]     # Previous generic result

declare_res=[]      # Global list of result tensors
declare_name=[]     # Global of just result tensor names
declare_index=[]    # Global list of result indices

declare_inter=[]    # Global list of intermediates
declare_inter_index=[]    # Global list of intermediates
declare_inter_name=[]    # Global list of intermediates

declare_ten=[]
declare_ten_index=[]
declare_ten_name=[]

begin=False
end=False
old_spin_iter=[]

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
                        # TODO: use dict
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

# Print out intermediates
print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    print("tensor: %-20s, !Create{type:disk}" % (declare_inter[i]), file=f2)

# Print out code blocks
print(file=f2)
print('---- code("Test")', file=f2)
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

print("---- end", file=f2)

f2.close()
