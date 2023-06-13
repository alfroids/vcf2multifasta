# vcf2multifasta

> Haploid genomes only!

Create a FASTA file with multiple sequences from a VCF file and the reference genome to which all genomes in the VCF file were originally aligned.

## How to use

This script requires `PyVCF`, `pysam`, `biopython` and `pandas`, which can be installed via `pip`:
```
pip install pandas biopython pysam PyVCF
```
> `pysam` is only available for UNIX systems

### Arguments

```
(req) --vcf            Path to the VCF file. There must also exist an indexation file with the same name and '.tbi' extension.
(req) --ref            Path to the reference FASTA file.
(req) --chr            Chromosome to be converted to Multi-FASTA. This name must be tha same both in the VCF and reference FASTA files.
(opt) --start          Position of the reference genome from where to start FASTA sequences. If not specified, sequence will start on the beginning of the specified chromosome.
(opt) --en             Position of the reference genome on which to end FASTA sequences. If not specified, sequence will end on the end of the specified chromosome.
(opt) --des            Default description to be added to all multi FASTA entries. If not specified, no description will be added.
(opt) --out            Path to the output FASTA file. If not specified, the path will be the same as the VCF file. Output file will have the '.fa' extension.
(opt) --gaps, -g       Insert naive gaps on the output sequences. Gaps will be added such that SNPs in each position have the same size for every individual in the VCF file. This results in sequences of the same length in the multi FASTA file.
(opt) --verbose, -v    Verbose mode.
```
