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
print("// Intermediate to pass to Kext", file=f2)
if (not kext):
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
else:
    kext_temp.seek(0)
    for line in kext_temp:
        print(line.strip(), file=f2)

    print("store INTpp[abij]",file=f2)

# Transform K ext from internal indicies to closed and active
if multi:
    print("", file=f2)
    print("", file=f2)
    print('---- code("Transform_K")', file=f2)
    print("alloc K4E1[abij], K4E2[abpi]", file=f2)
    print("alloc INTpp1[abmi]", file=f2)
    print("load deltaci[im], deltaai[pm]", file=f2)
    print("load K4C[abmn]", file=f2)
    print(".INTpp1[abmj] += K4C[abmn] deltaci[jn]", file=f2)
    print(".K4E1[abij] += INTpp1[abmj] deltaci[im]", file=f2)
    print(".K4E2[abpi] += INTpp1[abmi] deltaai[pm]", file=f2)
    print("drop K4C[abmn]", file=f2)
    print("drop deltaai[pm], deltaci[im]", file=f2)
    print("drop INTpp1[abmi]", file=f2)
    print("store K4E2[abpi], K4E1[abij]", file=f2)

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

print(file=f2)
print(file=f2)


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
