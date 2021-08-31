#!/usr/bin/env python3

# =========================================================================================
# This file is used to pre-proccess the bcontr.tmp file and find simple simplifications.
# So far this includes collecting lines together which are the same and introducing
# a factor.
#
#  largely extended version: Aug. 2021 by A.K.
#     now extensive parsing of lines summation of identical terms
#
# =========================================================================================
import argparse     # Parse arguments
import copy


def get_cntidx(cl_item):
    # entries 2 and 4 are indices, indentify those that are common
    cntidx = []
    for ch in cl_item[2]:
        if ch in cl_item[4]:
            cntidx.append(ch)
    return cntidx


def make_typelist(cntidx):
    # analyse cntidx and sort indices of same contraction space to alphabet. order
    hole     = ['i','j','k','l','m','n','o']
    particle = ['a','b','c','d','e','f','g','h']
    valence  = ['p','q','r','s','t','u','v','w','x','y','z']
    
    cntidxT = [] # type list
    for ch in cntidx:
        if ch in hole:
            cntidxT.append('c')
        if ch in particle:
            cntidxT.append('e')
        if ch in valence:
            cntidxT.append('a')

    return cntidxT
            

def get_canidx(cntidx,f2):

    cntidxN = copy.deepcopy(cntidx)
    cntidxT = make_typelist(cntidx)
    cntidxTT = make_typelist(cntidx)
    ##print("cntidxN (b): ",cntidxN,file=f2)
    ##print("cntidxT (b): ",cntidxT,file=f2)

    for idx_typ in ['c','e','a']:
        nidx = cntidxT.count(idx_typ)
        if nidx <= 1:
            continue
        
        idxmap = []
        for idx in range(len(cntidxT)):
            if cntidxT[idx] == idx_typ:
                idxmap.append(idx)

        ##print(idx_typ," idxmap: ",idxmap,file=f2)

        idxlst = []
        for idx in idxmap:
            idxlst.append(cntidxN[idx])

        ##print("idxlst:     ",idxlst,file=f2)

        idxlst.sort()

        ##print("idxlst (s): ",idxlst,file=f2)

        jdx = 0
        for idx in idxmap:
            cntidxN[idx] = idxlst[jdx]
            jdx += 1

    
    ##print("cntidxN (a): ",cntidxN,file=f2)
    ##print("cntidxT (a): ",cntidxT,file=f2)

    return cntidxN


def make_canidx(ctb,f2):
    # entries 2 and 4 are indices ... change contraction indices to canonical (= alphabet.) order

    ctbN = copy.deepcopy(ctb)
    cntidx = get_cntidx(ctb)
    cntidxC = get_canidx(cntidx,f2)

    ##print("cntidx:  ",cntidx,file=f2)
    ##print("cntidxC: ",cntidxC,file=f2)

    idxstr = []
    for ch in ctbN[2]:
        idxstr.append(ch)
    for kk in range(len(idxstr)):
        ch = idxstr[kk]
        if ch in cntidx:
            ii = cntidx.index(ch)
            idxstr[kk] = cntidxC[ii]
    ctbN[2] = "".join(idxstr)
    idxstr = []
    for ch in ctbN[4]:
        idxstr.append(ch)
    for kk in range(len(idxstr)):
        ch = idxstr[kk]
        if ch in cntidx:
            ii = cntidx.index(ch)
            idxstr[kk] = cntidxC[ii]
    ctbN[4] = "".join(idxstr)
    return ctbN


def same_index(tensor_name,idx1,idx2):

    # may have 01/23 symmetry (in case of J due to Hermiticity)
    sym_ten_list0123 = ["T2","Ym2","K","J","INTpp"]
    
    # quick decision?
    same = (idx1 == idx2)

    # special tensors?  
    if not same and tensor_name in sym_ten_list0123:
        # check if idx1 allows for this interchange
        typelist1 = make_typelist(idx1)
        ####print (tensor_name,": ",typelist1)
        if typelist1[0]==typelist1[1] and typelist1[2]==typelist1[3]:
            idx1ch = "".join([idx1[1],idx1[0],idx1[3],idx1[2]])
            same = (idx1ch == idx2)
            # no need to also test the index types of list2, if not compatibel, the above will be false anyway

        if not same and tensor_name == "J":
            if typelist1[0]==typelist1[1]:
                idx1ch = "".join([idx1[1],idx1[0],idx1[2],idx1[3]])
                same = (idx1ch == idx2)
            if not same and typelist1[2]==typelist1[3]:
                idx1ch = "".join([idx1[0],idx1[1],idx1[3],idx1[2]])
                same = (idx1ch == idx2)
                
    return same


def sum_up(clist):
    # sum up duplicate entries

    clistS = copy.deepcopy(clist)
    
    for ii in range(len(clistS)):
        if abs(float(clistS[ii][0]))<1e-12:
            continue
        
        for jj in range(ii+1,len(clistS)):

            # tensor names should be the same
            if clistS[ii][1] == clistS[jj][1] and clistS[ii][3] == clistS[jj][3]:
                # indices should be the same (including 1/2 permutation, which we test in same index for certain known Tensors)
                if same_index(clistS[ii][1],clistS[ii][2],clistS[jj][2]) and same_index(clistS[ii][3],clistS[ii][4],clistS[jj][4]):
                    clistS[ii][0] += float(clistS[jj][0])
                    clistS[jj][0] = 0.0

    clistN = []
    for ctb in clistS:
        if abs(float(ctb[0]))>1e-12:
            clistN.append(ctb)

    return clistN


def factor_out(clist,f2):
    # factorize expression ... we assume that we do not need to interchange the sequence of the factors
    symbols = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W']
    
    # compare factors ... assign symbolic labels A, B, ... for further analysis
    symlist = []
    for ctb in clist:
        symlist.append([float(ctb[0]),"","",True])
    trans_dict = {}  # translation dict
        
    idx = 0
    # start with left tensors
    for ii in range(len(clist)):
        equal = False
        for jj in range(ii):
            if (clist[ii][1]==clist[jj][1]):
                if same_index(clist[ii][1],clist[ii][2],clist[jj][2]):
                    symlist[ii][1] = symlist[jj][1]
                    equal = True
        if not equal:
            if (idx>=len(symbols)):
                print("Programming restriction: Extend symbol list!")
                quit()
            symlist[ii][1] = symbols[idx]
            trans_dict[symbols[idx]] = [clist[ii][1],clist[ii][2]]
            idx = idx+1
    # and now the right ones
    for ii in range(len(clist)):
        equal = False
        for jj in range(ii):
            if (clist[ii][3]==clist[jj][3]):
                if same_index(clist[ii][3],clist[ii][4],clist[jj][4]):
                    symlist[ii][2] = symlist[jj][2]
                    equal = True
        if not equal:
            if (idx>=len(symbols)):
                print("Programming restriction: Extend symbol list!")
                quit()
            symlist[ii][2] = symbols[idx]
            trans_dict[symbols[idx]] = [clist[ii][3],clist[ii][4]]
            idx = idx+1

    ###print("symbolic list: ",symlist,file=f2)
    ###print("transl. dict:  ",trans_dict,file=f2)

    # try a factorization
    # new list will have the structure:
    #  [ terms ] ;   term =  {factor1, factor2} ;  factor = [ summands ] ; summand = {fac, ten} 
    newlist = []
    for ii in range(len(symlist)):
        # skip terms that have been processed
        if not symlist[ii][3]:
            continue
        # take current term and find further ones
        term = {}
        fac1 = []
        fac2 = []
        fact = float(symlist[ii][0])  # let us see where to better put the factor
        s1 = {'fac': fact, 'ten': symlist[ii][1]}
        s2 = {'fac': 1.0, 'ten': symlist[ii][2]}  
        fac1.append(s1)
        fac2.append(s2)
        term["factor1"] = fac1
        term["factor2"] = fac2
        symlist[ii][3] = False
        do_left = False
        do_right = False
        init = True
        for jj in range(ii+1,len(symlist)):
            if not symlist[jj][3]:
                continue
            if (init or do_right) and symlist[ii][1]==symlist[jj][1]:       # same first tensor? then collect to right factor
                do_right = True
                if init:
                    s1["fac"] = 1.0
                    s2["fac"] = fact  #  put the factor now here (still refers to s2 above)
                    init = False

                # collect new term to s2
                s2 = {}  # new dict
                s2 = {'fac': symlist[jj][0], 'ten': symlist[jj][2]}
                fac2.append(s2)
                symlist[jj][3] = False

            if (init or do_left) and symlist[ii][2]==symlist[jj][2]:       # same second tensor? then collect to left factor
                do_left = True
                if init:
                    s1["fac"] = fact  #  put the factor now here (still refers to s1 above)
                    s2["fac"] = 1.0
                    init = False

                # collect new term to s1
                s1 = {}  # new dict
                s1 = {'fac': symlist[jj][0], 'ten': symlist[jj][1]}
                fac1.append(s1)
                symlist[jj][3] = False
                
        newlist.append(term)

    ###print("Factored list",newlist,file=f2)
    # backtranslate
    clistF = []
    for term in newlist:
        termF = {}
        fac1F = []
        fac2F = []
        # avoid that second terms start with "-"
        tmp = float(term["factor2"][0]["fac"])
        if (tmp<0.0):
            gfac = -1.0
        else:
            gfac = 1.0
        for summand1 in term["factor1"]:
            s1F = {}
            s1F["fac"] = summand1["fac"]*gfac
            s1F["name"] = trans_dict[summand1["ten"]][0]
            s1F["index"] = trans_dict[summand1["ten"]][1]
            fac1F.append(s1F)
        for summand2 in term["factor2"]:
            s2F = {}
            s2F["fac"] = summand2["fac"]*gfac
            s2F["name"] = trans_dict[summand2["ten"]][0]
            s2F["index"] = trans_dict[summand2["ten"]][1]
            fac2F.append(s2F)
        termF["factor1"] = fac1F
        termF["factor2"] = fac2F

        clistF.append(termF)
        
    return clistF


def make_line(r_ten,r_idx,clist):
    # process input into ITF syntac line
    lines = []

    for term in clist:
        line = "."+r_ten+"["+r_idx+"] "
        first = True
        sfac = 1.0
        for factor in ["factor1", "factor2"]:
            nsum = len(term[factor])
            init = True
            if factor == "factor2":
                line += " "
                sfac = 1.0  # modify only pre-factors of first factor
            for summand in term[factor]:
                if summand["name"]=="_One":
                    continue
                
                if first:
                    if (float(summand["fac"])<0.0):
                        sfac = -1.0
                        line += "-= "
                    else:
                        sfac = +1.0
                        line += "+= "
                    first = False
                fac = float(summand["fac"])*sfac
                if init and nsum > 1:
                    line += "("
                if not init and fac>=0:
                    line += " + "
                if not init and fac<0:
                    line += " - "
                if not init and (abs(fac)!=1.0):
                    line += str(abs(fac))+"*"
                if init and ((abs(fac)!=1.0) or fac<0.0):
                    line += str(fac)+"*"
                if init:
                    init = False
                line += summand["name"]+"["+summand["index"]+"]"

            if nsum > 1:
                line += ')'

        lines.append(line)

    return lines



def process_contribs(results,contribs,f2):

    lines = []

    for r_ten_idx in results:
        clist = contribs[r_ten_idx]
        r_ten,r_idx = r_ten_idx.split('-',1)

        ##print("ori: \n",clist,file=f2)

        for ii in range(len(clist)):
            # canonicalize contraction index
            ctbc = make_canidx(clist[ii],f2)
            clist[ii] = ctbc
            
        ##print("can: \n",clist,file=f2)

        clistS = sum_up(clist)

        ##print("summed: \n",clistS,file=f2)

        clistF = factor_out(clistS,f2)
            
        ##print("factored: \n",clistF,file=f2)

        lines_res = make_line(r_ten,r_idx,clistF)

        for line in lines_res:
            lines.append(line)
        
    return lines



# Parse arguments from gecco
parser = argparse.ArgumentParser(
                description="""Simplify ITF algo code, but collecting lines together which
                               are the same and proceed one another""",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input',default=None,help='ITF binary contraction file')
parser.add_argument('-o','--output',default=None,help='ITF binary contraction file')
parser.add_argument('-t','--tasks',dest='tasks',action='store_true',help='toggle task mode')
parser.set_defaults(tasks=False)
args = parser.parse_args()

if args.input is None:
    print("Error in python ITF processor: Must provide input file")
    exit(1)
if args.output is None:
    print("Error in python ITF processor: Must provide output file")
    exit(1)

inp = args.input
out = args.output
tasks = args.tasks

# Open bcontr.tmp file to read from
f1 = open(inp,"r")
# Open bcontr2.tmp file to read from
f2 = open(out,"w")

# Each line will be saved a compared to the previous one
old_line = ""
factor = 1.0

debug = False

# whether we currently are within a BEGIN/END section
in_section = False

# Read each line of bcontr.tmp
for line in f1:

    if ("Error" in line):
        print("Error in translating GeCCo code to ITF, check bcontr.tmp file for more details", file=f2)
        quit()

    if not in_section:
        if "CODE_BLOCK" in line:
            # just print this line and continue
            print(line, end="", flush=True, file=f2)
            continue
        if "BEGIN" in line:
            # print this line and change flag
            in_section = True
            print(line, end="", flush=True, file=f2)
            # initialize lists and dictionaries:
            results  = []
            contribs = {}
            prev_r_tensor = ""
            prev_r_index  = ""
            
            continue

        if not tasks:
            # if we arrive here, something is wrong
            print("Error: Found unexpected line:", file=f2)
            print(line, end="", flush=True, file=f2)
            quit()

    # at this point, in_section must be True
    if "END" in line:
        # block is over, process all collected data and print new line(s)
        if (False):
            print("final contribs:\n",contribs,file=f2)

        new_lines = process_contribs(results,contribs,f2)

        for nline in new_lines:
            print(nline, file=f2)
        
        print(line, end="", flush=True, file=f2)
        in_section = False
        continue

    # collect contents of BEGIN/END block

    # these are some lines from the old simplify.py, we keep them for the time being:
    # Remove unwanted lines before process.py
    # Removes the reference energy
    if ("ECC[]" in line and "K[]" in line): continue
    # Removes the CASSCF energy
    #if ("ECC[]" in line and "Dm[]" in line): continue
    if ("ECC[]" in line and "Ym1" in line and "f" in line): continue
    if ("ECC[]" in line and "Ym2" in line and "K" in line): continue

    if (debug):
        print("O: ",line,end="",file=f2)

    # for tasks just echo the file
    if tasks:
        print(line, end="", flush=True, file=f2)
        continue
        
    words = line.split()

    r_tensor = words[0].split('[',1)[0].replace('.','')
    r_index  = words[0].split('[',1)[1].replace(']','')
    # check if this is still the same as before:
    # for the way, we currently autogenerate the code, this check is sufficient
    if  r_tensor != prev_r_tensor or r_index != prev_r_index:
        current_tensor = r_tensor+'-'+r_index
        results.append(current_tensor)
        contribs[current_tensor] = []
        prev_r_tensor = r_tensor
        prev_r_index  = r_index
    
    # I do not think that ":" is generated ... better stop to see
    #sign   = words[1].replace('=','').replace(':','+')
    sign   = words[1].replace('=','')
    if sign == '+':  # global sign from += or -= statement
        gfac = 1.0
    elif sign == '-':
        gfac = -1.0
    else:
        print("Error: unexpected ':' or even more unexpected assigment line:", file=f2)
        print(line, end="", flush=True, file=f2)
        quit()
        

    paren_open = False      # whether a "(" has been found, but still no ")"
    look_for_sign = False   # whether to look next for "+" or "-"
    have_sfac = False       # whether a sign was found (to modifiy the factor)
    no_second_ten = True          # flag to False if there is a second tensor
    idx = 0
    sfac = +1.0
    in_tensor_index0 = []
    in_tensor_index1 = []
    for term in words[2:]:
        if look_for_sign:
            if term == "+":
                sfac = +1.0
                have_sfac = True
            elif term == "-":
                sfac = -1.0
                have_sfac = True
            else:
                print("Error: unexpected symbole (expected '+' or '-')",file=f2)
                print(line, end="", flush=True, file=f2)
                quit()
            look_for_sign = False
            continue
        if "*(" in term: # factor in front of parentheses?
            gfac = gfac*float(term.split('*',1)[0])  # update global factor
            term = term.split('*',1)[1]
        if "(" in term:
            if paren_open:
                print("Error: unexpected double parenthesis opening, inspect this line:", file=f2)
                print(line, end="", flush=True, file=f2)
                quit()
            paren_open = True
            look_for_sign = True
            term = term.replace('(','')
        if "*" in term:
            fac = float(term.split('*',1)[0])
            term = term.split('*',1)[1]
        else:
            fac = 1.0
        if have_sfac:
            fac *= sfac
            have_sfac = False
        if ")" in term:
            if not paren_open:
                print("Error: unexpected parenthesis closing, inspect this line:", file=f2)
                print(line, end="", flush=True, file=f2)
                quit()
            paren_open = False
            term = term.replace(')','')
        tensor = term.split('[',1)[0]
        index  = term.split('[',1)[1].replace(']','')
        if idx==0:
            in_tensor_index0.append([fac,tensor,index])
        else:
            no_second_ten = False
            in_tensor_index1.append([fac,tensor,index])
        if not paren_open:
            idx += 1

    # provide dummy contribution if there was no second tensor
    if no_second_ten:
        in_tensor_index1.append(['1.0','_One',''])
            
    # go over in_tensor_indexX lists an assemble terms contributing to current r_tensor
    for fac1,ten1,idx1 in in_tensor_index0:
        for fac2,ten2,idx2 in in_tensor_index1:
            contribs[current_tensor].append([gfac*float(fac1)*float(fac2),ten1,idx1,ten2,idx2])

    if (False):
        print("found so far:",file=f2)
        print("r_tensor: ",r_tensor,file=f2)
        print("r_index:  ",r_index,file=f2)
        print("contribs:\n",contribs,file=f2)
        


# Print the last reamining line of the file
print(old_line, end="", flush=True, file=f2)

# Clsoe the files
f1.close()
f2.close()
