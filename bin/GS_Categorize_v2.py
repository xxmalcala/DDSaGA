#!/usr/bin/python3

'''Categorizes the Germline-Soma architecture information from the BLASTN
report. This is not necessarily the best method, but it provides a lot
of sensitivity.

The loci are categorized as scrambled, non-scrambled, single-MDS,
multi-germline-loci, or possessing bad (inconclusive) pointer sequences. For
more complete definitions see Maurer-Alcalá et al. 2018 (mBio).'''

##__Updated__: 2020-12-16
##__Author__: Xyrus X. Maurer-Alcalá
##__e-mail__: maurerax@gmail.com; xmaurer-alcala@amnh.org

import sys

# Recursively removes hits that have 100% overlap with a longer alignment
def collapse_repetitive(hit_info, run_loop):
    bad_hits = []
    for n in range(len(hit_info)-1):
        n_fp = int(hit_info[n].split('\t')[3])
        n_tp = int(hit_info[n].split('\t')[4])
        n1_fp = int(hit_info[n+1].split('\t')[3])
        n1_tp = int(hit_info[n+1].split('\t')[4])
        if n_fp <= n1_fp and n_tp >= n1_tp:
            bad_hits.append(hit_info[n+1])
    if not bad_hits:
        return [i for i in hit_info if i not in bad_hits], False
    else:
        return [i for i in hit_info if i not in bad_hits], True


# Categorizes Germ-Soma hits as potentially scrambled, nonscrambled or
# inconclusive if given poor pointers
def checkIES(hit_info):
    # Simple "catch" for MDS-MDS orientations
    sign = []
    for n in range(len(hit_info)-1):
        sn_tp = int(hit_info[n].split('\t')[4])
        sn1_fp = int(hit_info[n+1].split('\t')[3])
        gn_fp = int(hit_info[n].split('\t')[-2])
        gn_tp = int(hit_info[n].split('\t')[-1])
        gn1_fp = int(hit_info[n+1].split('\t')[-2])

        # Pointer sequence length – minimum of 1bp and maximum of 25bp
        pointer_len = sn_tp - sn1_fp
        if 25 >= pointer_len >= 1:
            same_MDS = gn_tp - gn_fp
            between_MDS = gn1_fp - gn_tp
            sign += '+' if same_MDS > 0 else '-'

    ### Minimum IES size is currently set to 20bp - XXMA 2020-12-17
    ### Feel free to change it! Make sure it is the desired value - 1
            if abs(between_MDS) > 19: ###
                sign += '+' if between_MDS > 0 else '-'

    # Checks for a single orientation, multiple, or a "failure"
    if len(set(sign)) == 1:
        return 'nonscrambled'
    elif len(set(sign)) > 1:
        return 'scrambled'
    else:
        return 'bad-pointer'

# Cleans the data and categorizes it according to criteria described in
# Maurer-Alcalá, Knight and Katz (2018, mBIO)
def categorizeLoci(tsv_file, perc_len=0.6):
    # Will store the soma-germ information for quick access
    tsvDict = {}
    tsvDict_filt = {}

    # Will store the categorized the germ-soma hits based on their
    # germline-loci organization
    germ_arch = {'single':{}, 'scrambled':{}, 'nonscrambled':{}, 'multiloc':{},
        'bad-pointer':{}}

    # Only keep the relevant information:
    # germline locus, percent ID, soma start/end, germline start/end
    for i in open(tsv_file).readlines():
        h_pid = '\t'.join(i.split('\t')[1:4])
        summary_hit = '\t'.join(i.split('\t')[6:10])
        tsvDict.setdefault(i.split('\t')[0],[]).append(h_pid+'\t'+summary_hit)

    # Remove somatic sequences that have hits that represent < 60% of
    # their length
    poor_hit = [k for k, v in tsvDict.items() if sum([int(i.split('\t')[2])
        for i in v]) < perc_len*int(k.split('_')[3])]

    for k in poor_hit:
        del tsvDict[k]

    # Sort the hits based on their somatic start positions
    for k, v in tsvDict.items():
        v.sort(key=lambda x: int(x.split('\t')[-4]))
        temp = [i for i in v]
        if len(v) > 1:
            looping = True
            while looping == True:
                temp, looping = collapse_repetitive(temp, looping)
        tsvDict_filt[k] = temp

    # Store and remove the somatic sequences that hit a single germline
    # locus once
    for k, v in tsvDict_filt.items():
        if len(v) == 1:
            germ_arch['single'][k] = v[0]

    for k in germ_arch['single']:
        del tsvDict_filt[k]

    # Identify, store, and remove mulit-locus hits
    for k, v in tsvDict_filt.items():
        if len(set([i.split('\t')[0] for i in v])) > 1:
            germ_arch['multiloc'][k] = v

    for k in germ_arch['multiloc']:
        del tsvDict_filt[k]

    # Identify non-scrambled germline loci
    for k, v in tsvDict_filt.items():
        MDS_info = checkIES(v)
        germ_arch[MDS_info][k] = v

    return germ_arch


# Generates .tsv files with the temporarily categorized germ-soma loci
# These are templates for the alignment step (if the flag is set to run!),
# which will clean these data further!

# BUG-FIX (to-do) --- XXMA 09-22
# Other issues to fix: issue with output coordinates and indels, NSC outputs are
# not easy for inexperienced, merge multi-loc with scrambled (add notes column).
def convertToDF(germ_arch, tsv_file):
    # Header for the tables
    header = f'Soma\tGerm\tID\tAln\tS_start\tS_end\tG_start\tG_end\n'

    # Quick formatting prep for each of the spreadsheets to generate
    temp_sgl = [f'{k}\t{v}' for k, v in germ_arch['single'].items()]
    temp_nsc = [f'{k}\t{i}' for k, v in germ_arch['nonscrambled'].items() for i in v]
    temp_sc = [f'{k}\t{i}' for k, v in germ_arch['scrambled'].items() for i in v]
    temp_ml = [f'{k}\t{i}' for k, v in germ_arch['multiloc'].items() for i in v]
    temp_bp = [f'{k}\t{i}' for k, v in germ_arch['bad-pointer'].items() for i in v]

    sgl_tsv_out = tsv_file.replace('Master','Categorized').split('MASTER')[0]\
        +'SingleMDS.tsv'
    nsc_tsv_out = sgl_tsv_out.replace('SingleMDS','NonScrambled')
    sc_tsv_out = sgl_tsv_out.replace('SingleMDS','Scrambled')
    ml_tsv_out = sgl_tsv_out.replace('SingleMDS','MultiLoci')
    bp_tsv_out = sgl_tsv_out.replace('SingleMDS','BadPointers')

    with open(sgl_tsv_out,'w+') as v:
        v.write(header)
        v.write('\n'.join(temp_sgl))

    with open(nsc_tsv_out,'w+') as w:
        w.write(header)
        w.write('\n'.join(temp_nsc))

    with open(sc_tsv_out,'w+') as x:
        x.write(header)
        x.write('\n'.join(temp_sc))

    with open(ml_tsv_out,'w+') as y:
        y.write(header)
        y.write('\n'.join(temp_ml))

    with open(bp_tsv_out,'w+') as z:
        z.write(header)
        z.write('\n'.join(temp_bp))
