def print_inter(prev_lines):
    # Load, contract, drop tensors involved with intermediates
    for i in range(0, len(prev_lines)):
        inter_words=prev_lines[i].split()
        if "STIN" not in inter_words[2] or "STIN" not in inter_words[3]:
            load_ten="load "
            if "STIN" not in inter_words[2]:
                load_ten=load_ten + inter_words[2].split('*',1)[-1]
            if "STIN" not in inter_words[3] and inter_words[3].split('*',1)[-1] != inter_words[2].split('*',1)[-1]:
                # Do not load if an intermediate or if the same as previous loaded tensor
                load_ten=load_ten + ", " + inter_words[3].split('*',1)[-1] 
            print(load_ten, file=out)
        print(prev_lines[i].strip(), file=out)
        if "STIN" not in inter_words[2] or "STIN" not in inter_words[3]:
            drop_ten="drop "
            if "STIN" not in inter_words[3]:
                drop_ten=drop_ten + inter_words[3].split('*',1)[-1]
            if "STIN" not in inter_words[2] and inter_words[2].split('*',1)[-1] != inter_words[3].split('*',1)[-1]:
                drop_ten=drop_ten + ", " + inter_words[2].split('*',1)[-1]
            print(drop_ten, file=out)

def print_result(line, words):
    # Load tensors
    if "STIN" not in words[2]:
        print("load ", words[2].split('*',1)[-1], file=out)
    if "STIN" not in words[3]:
        print("load ", words[3].split('*',1)[-1], file=out)

    print(line.strip(), file=out)

    # Drop tensors
    if "STIN" not in words[3]:
        print("drop ", words[3].split('*',1)[-1], file=out)
    if "STIN" not in words[2]:
        print("drop ", words[2].split('*',1)[-1], file=out)



f=open("gecco.itfaa","r")
out=open("code.itfaa", "w+")

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

hole=['i','j','k','l']
particle=['a','b','c','d']
valence=['p','q','r','s','t','u','v','w']

for line in f:
    words=line.split()

    # Add contracted tensors to global list
    if words[2].split('*',1)[-1] not in declare_ten and "STIN" not in words[2]:
        index=list(words[2][words[2].find("[")+1:words[2].find("]")])

        generic=[]
        for i in range (0,len(index)):
            for j in range(0,len(particle)):
                if index[i] in particle[j]:
                    generic.append('e')
        for i in range (0,len(index)):
            for j in range(0,len(valence)):
                if index[i] in valence[j]:
                    generic.append('a')
        for i in range (0,len(index)):
            for j in range(0,len(hole)):
                if index[i] in hole[j]:
                    generic.append('c')

        declared=False
        for i in range(0, len(declare_ten)):
            if words[2].split('[',1)[0].split('*',1)[-1] == declare_ten_name[i]:
                if generic == declare_ten_index[i]:
                    # Generic index must be at same position as name it belongs to - dangerous! 
                    # Load previous tensor
                    declared=True
                    break
            else:
                continue

        if not declared:
            # Add result to global list
            declare_ten.append(words[2].split('*',1)[-1])
            declare_ten_index.append(generic)
            declare_ten_name.append(words[2].split('[',1)[0].split('*',1)[-1])

#    elif words[3].split('*',1)[-1] not in declare_ten and "STIN" not in words[3]:
#        declare_ten.append(words[3].split('*',1)[-1])

    if words[3].split('*',1)[-1] not in declare_ten and "STIN" not in words[3]:
        index=list(words[3][words[3].find("[")+1:words[3].find("]")])

        generic=[]
        for i in range (0,len(index)):
            for j in range(0,len(particle)):
                if index[i] in particle[j]:
                    generic.append('e')
        for i in range (0,len(index)):
            for j in range(0,len(valence)):
                if index[i] in valence[j]:
                    generic.append('a')
        for i in range (0,len(index)):
            for j in range(0,len(hole)):
                if index[i] in hole[j]:
                    generic.append('c')

        declared=False
        for i in range(0, len(declare_ten)):
            if words[3].split('[',1)[0].split('*',1)[-1] == declare_ten_name[i]:
                if generic == declare_ten_index[i]:
                    # Generic index must be at same position as name it belongs to - dangerous! 
                    # Load previous tensor
                    declared=True
                    break
            else:
                continue

        if not declared:
            # Add result to global list
            declare_ten.append(words[3].split('*',1)[-1])
            declare_ten_index.append(generic)
            declare_ten_name.append(words[3].split('[',1)[0].split('*',1)[-1])

    if "STIN" in words[0]:
        # Check if contraction forms and intermediate
        prev_lines.append(line)
        prev_inter.append(words[0].replace('.',''))
        if words[0].replace('.','') not in declare_inter:
            # Add intermedate to global list
            index=list(words[0][words[0].find("[")+1:words[0].find("]")])
    
            generic=[]
            for i in range (0,len(index)):
                for j in range(0,len(particle)):
                    if index[i] in particle[j]:
                        generic.append('e')
            for i in range (0,len(index)):
                for j in range(0,len(valence)):
                    if index[i] in valence[j]:
                        generic.append('a')
            for i in range (0,len(index)):
                for j in range(0,len(hole)):
                    if index[i] in hole[j]:
                        generic.append('c')
    
            declared=False
            for i in range(0, len(declare_inter)):
                if words[0].split('[',1)[0].replace('.','') == declare_inter_name[i]:
                    if generic == declare_inter_index[i]:
                        # Generic index must be at same position as name it belongs to - dangerous! 
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
        # Compare gneeric indicies to see if belongs in same block
        # + Load instead of alloc
        index=list(words[0][words[0].find("[")+1:words[0].find("]")])

        generic=[]
        for i in range (0,len(index)):
            for j in range(0,len(particle)):
                if index[i] in particle[j]:
                    generic.append('e')
        for i in range (0,len(index)):
            for j in range(0,len(valence)):
                if index[i] in valence[j]:
                    generic.append('a')
        for i in range (0,len(index)):
            for j in range(0,len(hole)):
                if index[i] in hole[j]:
                    generic.append('c')


        if words[0] != prev_res and generic != prev_generic:
            # Next result is different from previous, so close off block


            if prev_res != '#####':
                print("store", prev_res.replace('.',''), file=out)
                print(file=out)

#            if words[0].replace('.','') in declare_res:
#                # Load previous tensor
#                print("load", words[0].replace('.',''), file=out)
#            else:
#                # Add result to global list and alloc new result
#                declare_res.append(words[0].replace('.',''))
#                print("alloc ", words[0].replace('.',''), file=out)


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
                declare_res.append(words[0].replace('.',''))
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

            # Drop intermediates if needed
            if prev_inter:
                print("drop ", end="", flush=True, file=out)
                print(*list(reversed(prev_inter)),sep=", ", file=out)

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
f2=open("code.itfaa", "r")
tmp=f2.read()
f2.close()

f2=open("code.itfaa", "w")
# Declare tensros and index-spaces
print("---- decl", file=f2)
print("index-space: pqrstuvw , Active    , a", file=f2)
print("index-space: ijkl     , Closed    , c", file=f2)
print("index-space: gh       , ClosedF   , f", file=f2)
print("index-space: abcd     , External  , e", file=f2)
print("index-space: mno      , Internal  , i", file=f2)
print(file=f2)

for i in range(0, len(declare_ten)):
    print("tensor: ", declare_ten[i], file=f2)

print(file=f2)
for i in range(0, len(declare_res)):
    print("tensor: ", declare_res[i], file=f2)

print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_inter)):
    print("tensor: %-20s, !Create{type:disk}" % (declare_inter[i]), file=f2)

print(file=f2)
print('---- code("Test")', file=f2)
f2.write(tmp)
print(file=f2)
print("---- end", file=f2)

f2.close()
