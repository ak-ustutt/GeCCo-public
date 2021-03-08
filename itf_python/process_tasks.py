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


def add_to_global(word,declare_ten,declare_ten_index,declare_ten_name):
    # Add tensor to global lists, so it can be declared at the start of the algo file

    if word.split('*',1)[-1] not in declare_ten:

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
import os

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
output=open(outp, "w+")

# Open tempfile to catch INTpp
kext_temp = tempfile.TemporaryFile(mode='w+t')

# Open tempfile for Initalist_Residual
init_res_temp = tempfile.TemporaryFile(mode='w+t')
init_res = False
init_alloc = False

# TODO: Now all residuals are declared for the mutli cases, don't need all of them, all of the time
declare_res_multi=[
        "R:I[I]",
        "R:ac[pi]","R:ec[ai]","R:ea[ap]",
        "R:aacc[pqij]","R:aaac[pqri]",
        "R:eacc[apij]","R:eaac[apqi]","R:eaaa[apqr]",
        "R:eeac[abpi]","R:eeaa[abpq]"]      # List of all residuals in multireference case
declare_amp_multi=[
        "T:I[I]",
        "T:ac[pi]","T:ec[ai]","T:ea[ap]",
        "T:aacc[pqij]","T:aaac[pqri]",
        "T:eacc[apij]","T:eaac[apqi]","T:eaaa[apqr]",
        "T:eeac[abpi]","T:eeaa[abpq]"]      # List of all residuals in multireference case

declare_ten=[]          # Global list of tensors involved in binary contractions
declare_ten_index=[]    # Global list of tensor indicies
declare_ten_name=[]     # Global list of tensor names


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


# Flag if dealing with triples
triples = False


# Read each line of the file and add each tensor to a global list
# This list will be used to declare tensors in the final itf file
for line_o in f:

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

    # Decide which file to write to
    #if (special_begin and not special_end):
    #    out = kext_temp
    #else:
    #    out = output

    # Check if brackets in the binary contraction
    if (len(words)>=4):
        # Check if result tensor needs to be declared
        add_to_global(words[0].replace('.',''),declare_ten,declare_ten_index,declare_ten_name)
        # Check if first tensor needs to be declared
        add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)
        # Check if second tensor needs to be declared
        add_to_global(words[3],declare_ten,declare_ten_index,declare_ten_name)
    elif (len(words)==3):
        # Simple add or assign line
        add_to_global(words[0].replace('.',''),declare_ten,declare_ten_index,declare_ten_name)
        add_to_global(words[2],declare_ten,declare_ten_index,declare_ten_name)

    # Print lines to a tmp file. Each file represents a code/task block
    print(line, file=out)

output.close()
f.close()



# =========================================================================================
# Print final .itfaa file
# =========================================================================================

# Open and write file again so as to prepend the declaration of tensors
f2=open(outp, "r")
tmp=f2.read()
f2.close()

f2=open(outp, "w")
gecco = os.environ["GECCO_DIR"]

# init and save arrays, used in the task blocks
init = []
save = []


print_code_block('header', gecco, f2)
now = datetime.datetime.now()
print("// Created on:", now.strftime("%d-%m-%Y %H:%M"), file=f2)
print(file=f2)

# Print index-spaces
if (not multi):
    print_code_block('single_ref/index_spaces', gecco, f2)
else:
    print_code_block('multi_ref/index_spaces', gecco, f2)


# Print already existing tensor, ie. don't need !Create{type:disk}
declare_existing_tensors(declare_ten, "K-integral tensors", "K")

# If not declared already, declare integrals needed in calculating fock matricies etc.
combined = '\t'.join(declare_ten)
if multi:
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

declare_existing_tensors(declare_ten, "J-integral tensors", "J")
declare_existing_tensors(declare_ten, "Special integral tensors", "K4E",True)
if multi:
    print("tensor: K4C[abmn], K4C", file=f2)

#if (kext):
#    declare_existing_tensors(declare_ten, "Tensor to send to Kext", "INTpp",True)
#else:
#    print(file=f2)
#    print("// Tensor to send to Kext", file=f2)
#    print("tensor: INTpp[abij], INTpp", file=f2)

if (kext):
    declare_existing_tensors(declare_ten, "Tensor to send to Kext", "INTpp",True)
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
print("",file=f2)


print("// Residual tensors", file=f2)
for i in range(0, len(declare_ten)):
    if "ITIN" in declare_ten[i]:
        gindex = declare_ten[i].split('[')[0]+":"+"".join(generic_index(declare_ten[i]))
        gindex = gindex.split(':')[1]
        #print("tensor: R%-18s" % (":"+gindex+"["+ declare_ten[i].split('[')[1] + ", R:" + gindex), file=f2)
        print("tensor: R:" + gindex + "[" + declare_ten[i].split('[')[1] + ", " + "R:" + gindex, file=f2)
        init.append("R:" + gindex + "[" + declare_ten[i].split('[')[1])
        save.append("R:" + gindex + "[" + declare_ten[i].split('[')[1])

        # Add R to declare tensor
        declare_ten.append("R:"+gindex+"[]")

    elif "R" in declare_ten[i]:
        print("tensor: " + declare_ten[i] + ", " + declare_ten[i].split('[')[0], file=f2)
        init.append(declare_ten[i])
        save.append(declare_ten[i])
print("",file=f2)


if (any('R:eeeccc' in s for s in declare_ten)):
   triples = True

if (multi):
    print(file=f2)
    print("// Residuals not used in the code", file=f2)
    for i in range(0, len(declare_res_multi)):
        if any(declare_res_multi[i].split('[',1)[0]+"[" in s for s in declare_ten):
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


print("// Energy and DIIS scalars", file=f2)
for i in range(0, len(declare_ten)):
    if "ECC" in declare_ten[i]:
        print("tensor: %-10s" % (declare_ten[i] + ", " + declare_ten[i].split('[')[0]), file=f2)
        init.append(declare_ten[i])
        save.append(declare_ten[i])
        break
print("",file=f2)


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
    print("tensor: Ym1[pp], Ym1",file=f2)
    #if "Ym2" in combined:
    print("tensor: Ym2[pppp], Ym2",file=f2)

else:
    print_code_block('single_ref/declare_tensors', gecco, f2)


# Print out intermediates
print(file=f2)
print("// Intermediates", file=f2)
for i in range(0, len(declare_ten)):
    if "TIN" in declare_ten[i] or "X" in declare_ten[i]:
        if ("[]" in declare_ten[i]):
            print("tensor: %-18s !Create{type:scalar}" % (declare_ten[i] + ","), file=f2)
        else:
            # TODO: possible bug is ItfTasks.cpp, it stores a plain tensor. For now
            # declare everything as disk.
            #print("tensor: %-18s !Create{type:plain}" % (declare_ten[i] + ","), file=f2)
            print("tensor: %-18s !Create{type:disk}" % (declare_ten[i] + ","), file=f2)
        init.append(declare_ten[i])


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
print(file=f2)
print(file=f2)
print('---- code("Update_INTkx")', file=f2)
print('load T:eecc[abij]', file=f2)
print('drop T:eecc[abij]', file=f2)
print('', file=f2)
print('', file=f2)

print('---- task("Update_INTkx")', file=f2)
if (not kext):
    if multi:
        print("init INTpp1[abmi], INTpp[abmn]",file=f2)
        print("save INTpp1[abmi], INTpp[abmn]",file=f2)
        print(".INTpp1[abmj] += T:eecc[abij] deltaci[im]", file=f2)
        print(".INTpp1[abmi] += T:eeac[abpi] deltaai[pm]", file=f2)
        print(".INTpp[abmn] += INTpp1[abmj] deltaci[jn]", file=f2)
    else:
        print("init INTpp",file=f2)
        print("save INTpp",file=f2)
        print(".INTpp[abij] += T:eecc[abij]",file=f2)
else:
    print("init INTpp, ITIN[abij]",file=f2)
    print("save INTpp",file=f2)
    kext_temp.seek(0)
    for line in kext_temp:
        print(line.strip(), file=f2)

    print(".INTpp[abij] += ITIN[abij]",file=f2)
    print(".INTpp[abij] += ITIN[baji]",file=f2)

print('', file=f2)

if multi:
    print_code_block('multi_ref/transform_k', gecco, f2)
    print_code_block('multi_ref/transform_intk', gecco, f2)


## Print out Init_Residual
#if (initalise):
#    print(file=f2)
#    print(file=f2)
#    print('---- code("Init_Residual")', file=f2)
#    init_res_temp.seek(0)
#    for line in init_res_temp:
#        print(line.strip(), file=f2)


# Print out all code blocks generated from input
for block_name in code_blocks:
    code_block_tmp[block_name].seek(0)
    print('---- task("'+block_name+'")', file=f2)

    # TODO: hack, every code_block should have it's own init and save
    if block_name == 'Residual':
        print("init ", end="", flush=True, file=f2)
        print(*init, sep=", ", file=f2)
        print("save ", end="", flush=True, file=f2)
        print(*save, sep=", ", file=f2)

    for line in code_block_tmp[block_name]:
        print(line.rstrip(), file=f2)

    if block_name == 'Residual':
        # Symmetrise tensors
        for i in range(0, len(declare_ten)):
            if "ITIN" in declare_ten[i]:
                # Pick out specific tensor we wish to symmetrise
                if "".join(generic_index(declare_ten[i])) == 'eecc':
                    gindex = declare_ten[i].split('[')[0]+":"+"".join(generic_index(declare_ten[i]))
                    gindex = gindex.split(':')[1]
                    res_ten = "R:" + gindex + "[" + declare_ten[i].split('[')[1]

                    index = declare_ten[i].split('[')[1]
                    index = index.replace("]","")
                    index2 = ""
                    index2 += index[1:2]
                    index2 += index[0:1]
                    index2 += index[3:4]
                    index2 += index[2:3]

                    print("." + res_ten + " += " + declare_ten[i], file=f2)
                    print("." + res_ten + " += " + declare_ten[i].split('[')[0] + "[" + index2 + "]", file=f2)


    code_block_tmp[block_name].close()


# Print out residual equations
print(file=f2)
print('---- code("Residual")', file=f2)
print('load T:eecc[abij]', file=f2)
print('drop T:eecc[abij]', file=f2)
print(file=f2)

f2.write(tmp)



# Print out Generate_Fock_Matricies
if multi:
    print_code_block('multi_ref/generate_fock_matricies', gecco, f2)


# Print out Update_Amplitudes
if multi:
    print_code_block('multi_ref/form_dm3', gecco, f2)
    print_code_block('multi_ref/sblock', gecco, f2)
    print_code_block('multi_ref/transform_residual', gecco, f2)
    print_code_block('multi_ref/create_amplitude_update', gecco, f2)
    print_code_block('multi_ref/update_amplitudes', gecco, f2)
    print_code_block('multi_ref/construct_gs_overlap', gecco, f2)
    print_code_block('multi_ref/construct_s2', gecco, f2)
    print_code_block('multi_ref/construct_projected_s2', gecco, f2)
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
