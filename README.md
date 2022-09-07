# DSaGA
**D**elineating **S**omatic **a**nd **G**ermline **A**rchitecture 

## Description
DSaGA compares assembled somatic data (transcriptome or genome) against a reference germline genome assembly to determine the organization of somatic-destined sequences in the germline. This is largely directed towards ciliates (Alveolata: Ciliophora), but may be applicable elsewhere. This tool has largely been coupled with single-cell -omics data, but works well with "complete" somatic and germline genomes as well.

It takes FASTA DNA sequence as input, writing out tab-separated value tables summarizing the germline genome architecture.

## Dependencies
[Python 3.6+](https://www.python.org/downloads/)\
[BioPython](https://biopython.org/wiki/Download)\
[BLAST 2.12+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)\
[EMBOSS](http://emboss.sourceforge.net/)\
[CD-HIT](https://github.com/weizhongli/cdhit)

## Usage
Can be run with both somatic genomes or solely transcriptomic data. For detailed usage information and options linked to data type (genome or transcriptome):
```
$ python3 GS.py transcriptome --help
$ python3 GS.py genome --help
```

## Planned Updates
-   GFF3 outputs
-   Options for Visualization of Germ-Soma Architecture
-   Conda packaging!

### When using the tool in published research, please cite:

-   Maurer-Alcalá XX, Knight R, Katz LA. 2018. \"Exploration of the Germline Genome of the Ciliate Chilodonella uncinata through Single-Cell Omics (Transcriptomics and Genomics)\", *_mBio_* **9(1)**: e01836-17. [doi.org/10.1128/mBio.01836-17](https://doi.org/10.1128/mBio.01836-17)
    
 #### Further publications:
 -   Maurer-Alcalá XX, Yan Y, Pilling OA, Knight R, Katz LA. 2018. \"Twisted tales: insights into genome diversity of ciliates using single-cell ‘omics\", *_GBE_* **10(8)**: 1927-1938. [doi.org/10.1093/gbe/evy133](https://doi.org/10.1093/gbe/evy133)
  -   Smith SA, Maurer-Alcalá XX, Yan Y, Katz LA, Santoferra LF, McManus GB. 2020. \"Combined Genome and Transcriptome Analyses of the Ciliate Schmidingerella arcuata (Spirotrichea) Reveal Patterns of DNA Elimination, Scrambling, and Inversion\", *_GBE_* **12(9)**: 1616-1622. [doi.org/10.1093/gbe/evaa185](https://doi.org/10.1093/gbe/evaa185)



