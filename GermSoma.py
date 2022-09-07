#!/usr/bin/python3

##__Updated__: 2020-12-17
##__Author__: Xyrus X. Maurer-Alcalá
##__e-mail__: maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch

from bin import GS_Prep_v2 as pr
from bin import GS_Telo_v2 as tl
from bin import GS_BLAST_v2 as bl
from bin import GS_Categorize_v2 as ctg
from bin import GS_Align_v2 as aln

# import GS_Prep_v2 as pr
# import GS_Telo_v2 as tl
# import GS_BLAST_v2 as bl
# import GS_Categorize_v2 as ctg
# import GS_Align_v2 as aln


import argparse
from argparse import RawTextHelpFormatter
from argparse import SUPPRESS
import sys

def check_args():

	parser = argparse.ArgumentParser(description='\nThis program will identify '
        'and attempt to characterize\ngermline-soma architecture from an '\
        "organism's somatic\n(transcriptome/genome) and germline (partial/"\
        'complete)\nFASTA files.', formatter_class=RawTextHelpFormatter)

    # Set of required options for the script
	required_arg_group = parser.add_argument_group('Input-Output Options')

	required_arg_group.add_argument('--soma', '-s', action = 'store',
        metavar = 'somatic genome', help = ' Somatic Genome FASTA File\n\n')

	required_arg_group.add_argument('--transcriptome','-trans',
        action = 'store', metavar = 'transcriptome',
        help = ' Transcriptome FASTA File\n\n')

	required_arg_group.add_argument('--germ','-g', action = 'store',
        required = True, metavar = 'germline genome',
        help = ' Germline Genome FASTA File\n\n')

	required_arg_group.add_argument('--out','-o', action = 'store',
        required = True, metavar = '[output folder]',
        help=" Name of output folder\n\n")

	required_arg_group.add_argument('--tsv','-tsv', action = 'store',
        help = " Categorizes the user-provided TSV ONLY\n\n")

    # Set of FASTA-handling/parsing options for the script
	optional_arg_group = parser.add_argument_group('FASTA Handling Options')

	optional_arg_group.add_argument('--all_soma','-as', action = 'store_true',
        help = ' Use all somatic sequences, rather than just '\
        'telomere-bound seqs\n (default with transcriptomes)\n\n')

	optional_arg_group.add_argument('--telomere','-telo', action = 'store',
        help = " Single telomeric repeat (e.g., CCCCAAAA)\n\n")

	optional_arg_group.add_argument('--mask','-m',
        action = 'store_true',
        help = ' Masks telomeric ends in the SOMATIC assembly\n\n')

	optional_arg_group.add_argument('--soma_len','-slen', action = 'store',
        metavar = '[size in bp]', type = int, default = 300,
        help = " Length (in bp) of somatic sequences (default = 300)\n\n")

	optional_arg_group.add_argument('--germ_len','-glen', action = 'store',
        metavar = '[size in Kbp]', default = 10, type = int,
        help=" Length (in Kbp) of germline sequences (default = 10)\n\n")

	optional_arg_group.add_argument('--no_align','-naln', action = 'store_true',
        help=' Skips building alignments between soma\nand corresponding '\
        'germline loci\n\n')

    # Set of BLASTN-handling options for the script
	blast_arg_group = parser.add_argument_group('BLASTN Options')

	blast_arg_group.add_argument('--pid','-id', action = 'store',
        default = 97,
        help = " Percent identity for germline-soma hits (default = 97)\n\n")

	blast_arg_group.add_argument('--word_size','-w', action = 'store',
        default = 25, help = ' "Word-size" for BLASTN alignments '\
        '(default = 25)\n\n')

	blast_arg_group.add_argument('--threads','-p', action = 'store',
        default = 4, help = ' Number of threads for BLASTN (default = 4)\n\n')

    # Ensures that just script description provided if no arguments provided
	if len(sys.argv[1:]) == 0:
		print (parser.description+'\n')
		sys.exit()

	args = parser.parse_args()

	return args



def prep_data(out_folder, germ_fasta, soma_fasta, transcriptome, align, tsv,
        telomere, soma_size = 300, germ_size = 10):
    print(f'\n{"#"*50}\n\nPreparing folders and filtering data (by length)')
    if not tsv:
        pr.prep_folders_all(out_folder, telomere, align)
        return pr.prep_fasta(out_folder, germ_fasta, soma_fasta, transcriptome,
            telomere, soma_size, germ_size)
    if tsv:
        pr.prep_folders_tsv(out_folder, tsv)
        return None, None


def handle_telomeres(out_folder, fasta_file, telomere, is_soma, masking,
        all_soma):
    if is_soma:
        print(f'\n{"#"*50}\n\nIdentifying and handling telomere-bearing\n'\
        f'somatic sequences')

    else:
        print(f'\n{"#"*50}\n\nIdentifying and removing telomere-bearing '\
            f'likely\nsomatic contamination from germline sequences')

    return tl.clean_FASTA(out_folder, fasta_file, telomere, is_soma, masking,
        all_soma)


def blastn_germ_soma_steps(out_folder, germ_fasta, soma_fasta, pid = 97,
        threads = 4, word_size = 25):

    print(f'\n{"#"*50}\n\nCreating BLAST database and running BLASTN with\n'\
        f'the following parameters:\n  perc_identity: {pid}\n'\
        f'  word_size: {word_size}\n  ungapped')

    print(f'\nBLASTN may take a while, so do not be alarmed!')

    return bl.blastn_gs(germ_fasta, soma_fasta, threads, pid, word_size,
        out_folder)


def categorize_germ_soma(blastn_tsv, perc_aln_length=0.6):
    print(f'\n{"#"*50}\n\nCategorizing gemline-soma hits as:\n  "Single-MDS"'\
        f'\n  "Non-scrambled"\n  "Scrambled"\n  "Multi-germline-loci"\n'\
        f'  "Bad pointers"')

    print(f'\n\nNOTE: Minimum IES length (e.g. distance between\npointer '\
        f'sequences) is hard-coded as 20bp.\n\nFeel free to change this value'
        f' as needed\n      (line 53 in GS_Categorize.py)')

    gs_dict = ctg.categorizeLoci(blastn_tsv, perc_aln_length)
    ctg.convertToDF(gs_dict, blastn_tsv)


def check_germ_soma_by_align(out_folder, soma_fasta, germ_fasta):
    print(f'\n{"#"*50}\n\nUsing local alignments to check the validity\nof'\
    f' data in the "Bad pointers" and "Non-scrambled"\ncategories\n')

    aln.double_check_Nsc_BadPs(out_folder, soma_fasta, germ_fasta)

    print(f'\nFinished checking "Bad pointer" and "Non-scrambled" categories\n')


if __name__ == "__main__":

    args = check_args()
    args.align = not args.no_align
    args.out = f'{args.out}_GS'
    # print(args.mask)
    # print(args.all_soma)
    # sys.exit()

    if args.soma and args.transcriptome:
        print(f'Error: You cannot specify both a somatic genome (--soma) and '
        f'a transcriptome (--transcriptome)')
        sys.exit(2)

    if not (args.soma or args.transcriptome):
        print(f'Error: You must specify either a somatic genome (--soma) or '
        f'a transcriptome (--transcriptome)')
        sys.exit(2)

    if not args.tsv:
        size_filt_soma, size_filt_germ = prep_data(args.out, args.germ,
            args.soma, args.transcriptome, args.align, args.tsv,
            args.telomere, args.soma_len, args.germ_len, )

        if args.transcriptome:
            soma_for_blast = size_filt_soma

        if args.telomere:
            germ_for_blast = handle_telomeres(args.out, size_filt_germ,
                args.telomere, False, False, False)
            if not args.transcriptome:
                soma_for_blast = handle_telomeres(args.out, size_filt_soma,
                    args.telomere, True, args.mask, args.all_soma)
        else:
            soma_for_blast = size_filt_soma
            germ_for_blast = size_filt_germ

        args.tsv = blastn_germ_soma_steps(args.out, germ_for_blast, soma_for_blast,
            args.pid, args.threads, args.word_size)

    categorize_germ_soma(args.tsv)
    check_germ_soma_by_align(args.out, soma_for_blast, germ_for_blast)
