#!/usr/bin/python3

'''Removes telomere-bearing germline sequences (by default), captures
telomere-bearing somatic chomosomes, and can mask telomeric ends as well.'''

##__Updated__: 2020-12-15
##__Author__: Xyrus X. Maurer-Alcalá
##__e-mail__: maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch

import re
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO


class TeloCapture():
    def __init__(self, seq, telo_seq, name):
        self.seq = seq
        self.telomere = telo_seq
        self.description = name

    # Identifies presence/absence of telomeres and if telos are present at
    # both ends of the contig/scaffold
    def telo_eval(self):
        # telo_regex = []
        f_telo, r_telo = 0, 0
        self.ft_pos, self.rt_pos = [],[]

        # Prepares the telomere sequence
        telo_repeat = self.telomere.upper()*2
        step = int(float(len(self.telomere))/2)

        # telo_regex += [telo_repeat[:len(self.telomere)+step], telo_repeat[step:]]
        telo_regex = [telo_repeat[:len(self.telomere)+step], telo_repeat[step:]]

        # Returns True if telomeres are present within the first/last 70bp
        for t_re in list(set(telo_regex)):
            self.ft_pos += [m.start()+len(t_re) for m in re.finditer('(?='+t_re+')', str(self.seq)[:70])]
            self.rt_pos += [m.start()+len(t_re) for m in re.finditer('(?='+t_re+')', str(self.seq.reverse_complement())[:70])]
            if self.ft_pos:
                f_telo += 1
            if self.rt_pos:
                r_telo += 1
        ### must have two telo-counts per end to be safe! -- XXMA 2020-12-15
        if f_telo > 1 and r_telo > 1:
            self.telo_count = 2

        elif (f_telo > 1 and r_telo <= 1) or (f_telo <= 1 and r_telo > 1):
            self.telo_count = 1

        else:
            self.telo_count = 0

    # Masks telomeres if they are present in the soma and masking is on
    def telo_mask(self):
        temp_seq = str(self.seq)

        # Finds the telomere end position(s) if there is one at the specified end
        mask_fpos = max(self.ft_pos) if self.ft_pos else 0
        mask_rpos = max(self.rt_pos) if self.rt_pos else 0

        # Masks each end, for a fully telomere-masked contig/chromosome
        ft_masked = str(Seq(temp_seq.replace(temp_seq[:mask_fpos],'N'*mask_fpos)).reverse_complement())
        self.masked_seq = str(Seq(ft_masked.replace(ft_masked[:mask_rpos],'N'*mask_rpos)).reverse_complement())



def clean_FASTA(out_folder, fasta_file, telo_seq, is_soma = True,
        masking = True, all_soma = False):

    # Set up the telomere-number dictionary (stores data)
    telo_summary = {2:[],1:[],0:[]}
    f_file = fasta_file.rsplit('/')[-1]

    for seq_rec in SeqIO.parse(fasta_file,'fasta'):
        seq_telo = TeloCapture(seq_rec.seq, telo_seq, seq_rec.description)
        seq_telo.telo_eval()
        telo_summary[seq_telo.telo_count].append(seq_telo)

    if is_soma:
        # NoTelo.blahblahblah has the "less good" somatic data - this will be
        # updated for transcriptomes ––XXMA
        no_telo_fasta = f'NoTelo.{f_file}'
        with open(f'{out_folder}/Soma/Telomere_Eval/{no_telo_fasta}','w+') as w:
            for seq_rec in telo_summary[0]:
                w.write(f'>{seq_rec.description}_Len_{len(seq_rec.seq)}_No_Telo'\
                f'\n{seq_rec.seq}\n')

        # WithTelo.blahblahblah has the "good" somatic chromosomes. For the somatic
        # genomes, these are the data that are used by default (no telos are
        # currently ignored) ––XXMA
        with_telo_fasta = f'WithTelo.{f_file}'
        with open(f'{out_folder}/Soma/Telomere_Eval/{with_telo_fasta}','w+') as w:
                for seq_rec in telo_summary[1]:
                    w.write(f'>{seq_rec.description}_Len_{len(seq_rec.seq)}'\
                        f'_Single_Telo\n{seq_rec.seq}\n')

                for seq_rec in telo_summary[2]:
                    w.write(f'>{seq_rec.description}_Len_{len(seq_rec.seq)}'\
                        f'_Both_Telo\n{seq_rec.seq}\n')

        # For somatic chromosomes bearing telomeres, these will be masked with
        # N's. This is HIGHLY recommended ––XXMA
        if masking:
            masked_telo_fasta = f'MaskedTelo.{f_file}'
            with open(f'{out_folder}/Soma/Telomere_Eval/{masked_telo_fasta}',
                    'w+') as w:

                for seq_rec in telo_summary[1]:
                    seq_rec.telo_mask()
                    mask_count = seq_rec.masked_seq.count('N')
                    masked_len = f'{len(seq_rec.seq)-mask_count}'
                    w.write(f'>{seq_rec.description}_Len_{masked_len}_Single_'\
                        f'Masked_Telo\n{seq_rec.masked_seq}\n')

                for seq_rec in telo_summary[2]:
                    seq_rec.telo_mask()
                    mask_count = seq_rec.masked_seq.count('N')
                    masked_len = f'{len(seq_rec.seq)-mask_count}'
                    w.write(f'>{seq_rec.description}_Len_{masked_len}_Both_'\
                        f'Masked_Telo\n{seq_rec.masked_seq}\n')

            merge_fasta = f'AllSoma.WithMask.{f_file}'

            merge_cmd = f'cat {out_folder}/Soma/Telomere_Eval/'\
                f'{masked_telo_fasta} {out_folder}/Soma/Telomere_Eval/'\
                f'{no_telo_fasta} > {out_folder}/Soma/Telomere_Eval/{merge_fasta}'

            if all_soma:
                status = subprocess.call(merge_cmd, shell=True)
                return f'{out_folder}/Soma/Telomere_Eval/{merge_fasta}'
            else:
                return f'{out_folder}/Soma/Telomere_Eval/{masked_telo_fasta}'
        else:

            merge_fasta = f'AllSoma.NoMask.{f_file}'
            merge_cmd = f'cat {out_folder}/Soma/Telomere_Eval/{with_telo_fasta}'\
                f'{out_folder}/Soma/Telomere_Eval/{no_telo_fasta} > '\
                f'{out_folder}/Soma/Telomere_Eval/{merge_fasta}'

            if all_soma:
                status = subprocess.call(merge_cmd, shell=True)
                return f'{out_folder}/Soma/Telomere_Eval/{merge_fasta}'
            else:
                f'{out_folder}/Soma/Telomere_Eval/{with_telo_fasta}'


        #     return (out_folder+'/Soma/Telomere_Eval/NoTelo.'+fasta_file.rsplit('/')[-1],
        #         out_folder+'/Soma/Telomere_Eval/MaskedTelo.'+fasta_file.rsplit('/')[-1])
        #
        # else:
        #     return (out_folder+'/Soma/Telomere_Eval/NoTelo.'+fasta_file.rsplit('/')[-1],
        #         out_folder+'/Soma/Telomere_Eval/WithTelo.'+fasta_file.rsplit('/')[-1])

    # Write data out for the germline genome
    if not is_soma:
        no_telo_fasta = f'NoTelo.{f_file}'

        # NoTelo.blahblahblah has the "good" putative germline loci ––XXMA
        with open(f'{out_folder}/Germ/Telomere_Eval/{no_telo_fasta}','w+') as w:
            for seq_rec in telo_summary[0]:
                w.write(f'>{seq_rec.description}_No_Telo\n{seq_rec.seq}\n')

        # WithTelo.blahblahblah has the "bad" putative germline loci, which harbor
        # telomeres. This is a silly distinction, but as of 2020 there are no
        # complete (telomere-capped) germline genomes. ––XXMA
        with_telo_fasta = f'WithTelo.{f_file}'
        with open(f'{out_folder}/Germ/Telomere_Eval/{with_telo_fasta}','w+') as w:
            for seq_rec in telo_summary[1]:
                w.write(f'>{seq_rec.description}_Single_Telo\n{seq_rec.seq}\n')

            for seq_rec in telo_summary[2]:
                w.write(f'>{seq_rec.description}_Both_Telo\n{seq_rec.seq}\n')

        return f'{out_folder}/Germ/Telomere_Eval/{no_telo_fasta}'

        # return out_folder+'/Germ/Telomere_Eval/NoTelo.'+fasta_file.rsplit('/')[-1]
