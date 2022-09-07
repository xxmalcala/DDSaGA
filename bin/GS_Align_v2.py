#!/usr/bin/env python3
# coding: utf-8

'''STUFFFFFFF! I will get to this when I finishe the beast!'''

##__Updated__: 2020-12-17
##__Author__: Xyrus Maurer-Alcalá; maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch

import fnmatch
import os
import re
import subprocess
import more_itertools as mit

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Emboss.Applications import WaterCommandline


def printProgressBar(iteration, total, prefix = '', suffix = '',
	decimals = 1, length = 50, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 *
        (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def parseTSVs(out_folder, tsv_file):
    coord_dict = {}
    for i in open(tsv_file).readlines()[1:]:
        soma_germ = i.split('\t')[0]+'\t'+i.split('\t')[1]
        germ_fp = int(i.split('\t')[-2])
        germ_tp = int(i.split('\t')[-1])
        coord_dict.setdefault(soma_germ,[])
        coord_dict[soma_germ] += [germ_fp, germ_tp]
    return coord_dict


def waterAlign(out_folder, soma, germ):
    aln_out = f'{out_folder}/Alignments/{soma}_XX_{germ}.water'
    water_cmd = WaterCommandline()
    water_cmd.asequence = f'{out_folder}/Alignments/s1.fas'
    water_cmd.bsequence = f'{out_folder}/Alignments/s2.fas'
    water_cmd.gapopen = 10
    water_cmd.gapextend = 0.5
    water_cmd.aformat = 'fasta'
    water_cmd.outfile = aln_out
    # water_cmd = f'water -asequence {out_folder}/Alignments/s1.fas '\
    #     f'-bsequence {out_folder}/Alignments/s2.fas -gapopen 10 '\
    #     f'-gapextend 0.5 -aformat fasta -out {aln_out}'
    stdout, stderr = water_cmd()


def prepTempFasta(out_folder, coord, soma_seqs, germ_seqs):
    total_aln = len(coord.keys())
    aln_num = 0
    for k, v in coord.items():
        soma, germ = k.split('\t')
        # Coordinates to slice up the germline-locus. Intended to help reduce
        # the computational burden of the Waterman-Smith alignment implemented
        # by EMBOSS.
        fp = max([min(v)-200, 0])
        tp = max(v)+200
        # Grabs the somatic sequence
        s1_soma = f'>{soma}\n{soma_seqs[soma]}'
        if v[0] > v[1]:
            # Re-orient the germline locus to be in the same "reading frame"
            # as the corresponding soma sequence for the alignment
            s2_germ = f'>{germ}_ReverseComplement\n'\
                f'{germ_seqs[germ][fp:tp].reverse_complement()}'
        else:
            s2_germ = f'>{germ}\n{germ_seqs[germ][fp:tp]}'
        with open(f'{out_folder}/Alignments/s1.fas', 'w+') as w:
            w.write(s1_soma)
        with open(f'{out_folder}/Alignments/s2.fas', 'w+') as w:
            w.write(s2_germ)
        waterAlign(out_folder, soma, germ)
        aln_num += 1
        printProgressBar(aln_num, total_aln)


def parseWater(aln_fas):
    # Read the water alignment
    aln_in = {i.description:str(i.seq) for i in SeqIO.parse(aln_fas,'fasta')}
    soma = list(aln_in.values())[0]
    germ = list(aln_in.values())[1]
    # Get the positions of every gap character present in the somatic sequence
    gap_pos = [m.start() for m in re.finditer('-',soma)]
    # Group the gap positions into consecutive ranges
    temp_gaps = [list(group) for group in mit.consecutive_groups(gap_pos)]
    # Convert the ranges into discrete tuples for IES boundaries!
    cluster_gap = [(i[0],i[-1]) for i in temp_gaps if len(i) > 20]
    # Get the corresponding MDS boundaries!
    cluster_non_gap = []
    if len(cluster_gap) == 0:
        cluster_non_gap.append((0, len(soma)))
    else:
        for n in range(len(cluster_gap)):
            if n == 0:
                cluster_non_gap.append((0,cluster_gap[n][0]-1))
            if n == len(cluster_gap)-1:
                cluster_non_gap.append((cluster_gap[n][-1]+1,len(soma)))
            else:
                f = cluster_gap[n][-1]+1
                t = cluster_gap[n+1][0]-1
                cluster_non_gap.append((f,t))
    # Check if the alignment suggests a single large MDS or multiple MDSs
    if len(cluster_non_gap) == 1:
        return 'SingleMDS', (list(aln_in.keys())[0], list(aln_in.keys())[1])
    else:
        return 'Nonscrambled', (list(aln_in.keys())[0], list(aln_in.keys())[1])


def double_check_Nsc_BadPs(out_folder, soma_fasta, germ_fasta):

    for filename in os.listdir(f'{out_folder}/BlastReports/Categorized/'):
        if fnmatch.fnmatch(filename, '*.BadPointers.*'):
            bp_tsv = f'{out_folder}/BlastReports/Categorized/{filename}'
        elif fnmatch.fnmatch(filename, '*.NonScrambled.*'):
            nsc_tsv = f'{out_folder}/BlastReports/Categorized/{filename}'

    check_coord = parseTSVs(out_folder, bp_tsv)
    check_coord.update(parseTSVs(out_folder, nsc_tsv))

    soma_seqs = {i.description:i.seq for i in SeqIO.parse(soma_fasta,'fasta')}
    germ_seqs = {i.description:i.seq for i in SeqIO.parse(germ_fasta,'fasta')}

    prepTempFasta(out_folder, check_coord, soma_seqs, germ_seqs)

    os.system(f'rm {out_folder}/Alignments/s1.fas {out_folder}/Alignments/s2.fas')

    re_categorize = {'SingleMDS':[], 'Nonscrambled':[]}

    for f in os.listdir(out_folder+'/Alignments/'):
        if f.endswith('.water'):
            aln_file = out_folder+'/Alignments/'+f
            x = parseWater(aln_file)
            re_categorize[x[0]].append(x[1])

    additional_single_mds = f'{out_folder}/BlastReports/Categorized/'\
        f'Soma.GermHits.AdditionalSingleMDS.tsv'

    additional_non_scram = f'{out_folder}/BlastReports/Categorized/'\
        f'Soma.GermHits.UpdatedNonScrambled.tsv'

    with open(additional_single_mds, 'w+') as w:
        w.write('Soma\tGerm\n')
        for pair in re_categorize['SingleMDS']:
            w.write('\t'.join(pair))

    with open(additional_non_scram, 'w+') as w:
        w.write('Soma\tGerm\n')
        for pair in re_categorize['Nonscrambled']:
            w.write('\t'.join(pair))
