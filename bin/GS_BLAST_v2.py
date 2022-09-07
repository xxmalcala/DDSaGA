#!/usr/bin/python3

##__Updated__: 2020-12-17
##__Author__: Xyrus X. Maurer-Alcalá
##__e-mail__: maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch

import subprocess
import sys


def format_germ(germ_fasta, out_folder):
    mk_germ_blastdb_nuc = f'makeblastdb -in {germ_fasta} -dbtype nucl -out '\
        f'{out_folder}/db_Blast/germ_db'
    status = subprocess.call(mk_germ_blastdb_nuc, shell=True,
        stdout=subprocess.DEVNULL)

    if status > 0:
        print(f'Error: Failed to make BLASTN database. Exit code {status}')
        sys.exit()
    else:
        return f'{out_folder}/db_Blast/germ_db'


def blastn_gs(germ_fasta, soma_fasta, threads, pid, wrd, out_folder):

    germ_db = format_germ(germ_fasta, out_folder)

    out_blast_tsv = f'{out_folder}/BlastReports/Master/Soma.GermHits.'\
        f'MASTER.pid{pid}.wrd{wrd}.tsv'

    # Originally used the BLAST archive file, but it apparently is renaming
    # some sequences in the reformatting/conversion and I am (honestly) too
    # lazy to figure out why... So ignore this chunk and use the TSV instead!
    # Lazy XXMA 2020-12-17

    # blastn_cmd = f'blastn -num_threads {threads} -query {soma_fasta} -db '\
    #     f'{germ_db} -perc_identity {pid} -ungapped -word_size {wrd} -outfmt 11 '\
    #     f'-out {out_blast_archive}'

    blastn_cmd = f'blastn -num_threads {threads} -query {soma_fasta} -db '\
        f'{germ_db} -perc_identity {pid} -ungapped -word_size {wrd} -outfmt 6 '\
        f'-out {out_blast_tsv}'

    status = subprocess.call(blastn_cmd, shell=True)

    if status > 0:
        print(f'Error: BLASTN failed. Exit code {status}')
        sys.exit()
    else:
        return out_blast_tsv


### Ignore this function! For some reason, the reformatting of the BLAST
### archive file is introducing artefacts/renaming sequences...
def blast_reformatting(bl_asn):
    outfmt6_txt = bl_asn.replace('asn','fmt6.tsv')

    convertTo6 = f'blast_formatter -archive {bl_asn} -outfmt 6 -out '\
        f'{outfmt6_txt}'

    status = subprocess.call(convertTo6, shell=True)

    if status > 0:
        print(f'Error: Reformatting BLAST archive file to format 6 '\
            f'(tab-separated values) failed. Exit code {status}')
        sys.exit()
    else:
        return outfmt6_txt
