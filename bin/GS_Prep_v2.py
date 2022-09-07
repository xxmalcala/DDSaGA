#!/usr/bin/python3

##__Updated__: 2020-12-16
##__Author__: Xyrus X. Maurer-Alcalá
##__e-mail__: maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch

import errno
import os
import shutil

from Bio import SeqIO
from Bio.SeqUtils import GC

# Creates various directories
class PrepFolders():

    @staticmethod
    # Germline-related folders
    def make_germ_folders(out_folder, telomere):
        germ_folder = ['/Germ/Original', '/Germ/Spreadsheets']
        if telomere:
            germ_folder.append('/Germ/Telomere_Eval')
        for folder_path in germ_folder:
            os.makedirs(out_folder+folder_path, exist_ok=True)

    @staticmethod
    # Soma-related folders - same for transcriptomes AND somatic genomes
    # Currently unable to process transcriptomes AND somatic genomes in the
    # same analyses – XXMA 2020-12-15
    def make_soma_folders(out_folder, telomere):
        soma_folder = ['/Soma/Original', '/Soma/Spreadsheets']
        if telomere:
            soma_folder.append('/Soma/Telomere_Eval')
        for folder_path in soma_folder:
            os.makedirs(out_folder+folder_path, exist_ok=True)

    @staticmethod
    # Folder to harbor the BLAST+ formatted germline database
    def make_BLAST_db_folder(out_folder):
        try:
            os.mkdir(out_folder+'/db_Blast')
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_folder+'/db_Blast'):
                pass

    @staticmethod
    # Folders for the master BLASTN report and its processed/categorized copies
    def make_BLAST_reports_folder(out_folder):
        report_folders = ['/BlastReports/Master', '/BlastReports/Categorized']
        for folder_path in report_folders:
            os.makedirs(out_folder+folder_path, exist_ok=True)

    @staticmethod
    # Folder for the WATER alignments
    def make_alignment_folder(out_folder):
        try:
            os.mkdir(out_folder+'/Alignments')
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(out_folder+'/Alignments'):
                pass

class BackUpFiles():
    @staticmethod
    def backup_germ_fasta(out_folder, germ):
        # os.system(f'cp {germ} {out_folder}/Germ/Original/')
        shutil.copy2(germ, out_folder+'/Germ/Original/')

    @staticmethod
    def backup_soma_fasta(out_folder, soma):
        # os.system(f'cp {soma} {out_folder}/Soma/Original/')
        shutil.copy2(soma, out_folder+'/Soma/Original/')

    @staticmethod
    def backup_tsv(out_folder, tsv):
        # os.system(f'cp {tsv} {out_folder}/BlastReports/Master/')
        shutil.copy2(tsv, out_folder+'/BlastReports/Master/')


def size_filt_FASTA(out_folder, fasta_file, min_size, telomere = None, transcriptome = False,
        soma = True):
    seq_names = {}
    seq_count = 1

    if soma:

    # New output fasta names
        soma_out = f'{out_folder}.SOMA.{min_size}bp.Renamed.fas'

        with open(out_folder+'/Soma/'+soma_out, 'w+') as w:
            temp_seqs = [seq_rec for seq_rec in SeqIO.parse(fasta_file,'fasta')
                if len(seq_rec.seq) >= min_size]

            temp_seqs.sort(key=lambda x: -len(x.seq))

            for seq_rec in temp_seqs:
                if transcriptome:
                    new_name = f'Transcript_{seq_count}_Len_{len(seq_rec.seq)}'

                elif not telomere:
                    new_name = f'Soma_{seq_count}_Len_{len(seq_rec.seq)}'
                else:
                    new_name = f'Soma_{seq_count}'

                seq_names[seq_rec.description] = [
                    new_name,
                    f'{GC(seq_rec.seq):.2f}',
                    f'{len(seq_rec.seq)}']

                w.write(f'>{new_name}\n{seq_rec.seq.upper()}\n')
                seq_count += 1

        with open(out_folder+'/Soma/Spreadsheets/NameConversion'\
                '.SomaInitialFilt.tsv','w+') as x:
            x.write('Original Name\tSimplified Name\tGC Content\tLength (Bp)\n')
            for k, v in seq_names.items():
                x.write(k+'\t'+'\t'.join(v)+'\n')
        return out_folder+'/Soma/'+soma_out

    if not soma:
        germ_out = f'{out_folder}.GERM.{min_size}Kbp.Renamed.fas'
        with open(out_folder+'/Germ/'+germ_out, 'w+') as w:
            temp_seqs = [seq_rec for seq_rec in SeqIO.parse(fasta_file,'fasta')
                if len(seq_rec.seq) >= min_size*1000]
            temp_seqs.sort(key=lambda x: -len(x.seq))
            for seq_rec in temp_seqs:
                new_name = f'Germ_{seq_count}_Len_{len(seq_rec.seq)}'
                seq_names[seq_rec.description] = [
                    new_name,
                    f'{GC(seq_rec.seq):.2f}',
                    f'{len(seq_rec.seq)}']
                w.write(f'>{new_name}\n{seq_rec.seq.upper()}\n')
                seq_count += 1

        with open(out_folder+'/Germ/Spreadsheets/NameConversion'\
                '.GermInitialFilt.tsv','w+') as x:
            x.write('Original Name\tSimplified Name\tGC Content\tLength (Bp)\n')
            for k, v in seq_names.items():
                x.write(k+'\t'+'\t'.join(v)+'\n')
        return out_folder+'/Germ/'+germ_out


def prep_folders_all(out_folder, telomere, align):
    PrepFolders.make_soma_folders(out_folder, telomere)
    PrepFolders.make_germ_folders(out_folder, telomere)
    PrepFolders.make_BLAST_reports_folder(out_folder)
    PrepFolders.make_BLAST_db_folder(out_folder)

    if align:
        PrepFolders.make_alignment_folder(out_folder)


def prep_folders_tsv(out_folder):
    PrepFolders.make_BLAST_reports_folder(out_folder)
    BackUpFiles.backup_tsv(out_folder, tsv)


def prep_fasta(out_folder, germ_fasta, soma_fasta, transcriptome,
    telomere, soma_size, germ_size):

    if transcriptome != None:
        soma_fasta = transcriptome
    BackUpFiles.backup_germ_fasta(out_folder, germ_fasta)
    BackUpFiles.backup_soma_fasta(out_folder, soma_fasta)
    # print(transcriptome)
    filt_soma = size_filt_FASTA(out_folder, soma_fasta, soma_size,
        telomere, transcriptome)
    filt_germ = size_filt_FASTA(out_folder, germ_fasta, germ_size, soma = False)

    return filt_soma, filt_germ
