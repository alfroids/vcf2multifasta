#!/bin/python3

import argparse
import vcf  # PyVCF
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def replace_char(string: str, char_at: int, replace_with: str) -> str:
    return string[:char_at] + replace_with + string[char_at + 1 :]


def initialize_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", help="Path to the VCF file.", required=True, type=str)
    p.add_argument(
        "--ref", help="Path to the reference FASTA file.", required=True, type=str
    )
    p.add_argument(
        "--chr",
        help="Chromosome to be converted to Multi-FASTA.",
        required=True,
        type=str,
    )
    p.add_argument(
        "--start",
        help="Position of the reference genome from where to start FASTA sequences.",
        required=False,
        type=int,
        nargs="?",
        default=1,
    )
    p.add_argument(
        "--end",
        help="Position of the reference genome on which to end FASTA sequences.",
        required=False,
        type=int,
        nargs="?",
        default=-1,
    )
    # p.add_argument(
    #     "--des",
    #     help="Default description for all FASTA entries.",
    #     required=False,
    #     type=str,
    #     nargs="+",
    #     default="",
    # )
    p.add_argument(
        "--out",
        help="Path to the output FASTA file.",
        required=False,
        type=str,
        nargs="?",
        default="",
    )
    return p


parser = initialize_parser()
args = parser.parse_args()


# Read VCF and create SNPs DataFrame ##########################################
vcf_reader = vcf.Reader(open(args.vcf, "rb"))


if args.end > -1:
    vcf_reader = vcf_reader.fetch(args.chr, start=args.start, end=args.end)
else:
    vcf_reader = vcf_reader.fetch(args.chr, start=args.start)

cols = ["POS"] + vcf_reader.samples
snps = {c: [] for c in cols}


for record in vcf_reader:
    snps["POS"].append(record.POS)
    ref = record.REF

    for sample in record.samples:
        if sample.gt_bases is None:
            snps[sample.sample].append(record.REF)
        else:
            snps[sample.sample].append(sample.gt_bases)

SNP = pd.DataFrame(snps)


# Read reference FASTA and extract the sequence at the specified range ########
fasta_dict = SeqIO.index(args.ref, "fasta")
chrom = fasta_dict[args.chr]
locus = chrom[args.start : args.end + 1]
seq = str(locus.seq)


# Substitute SNPs into reference sequence and output multi-FASTA file #########
if args.out != "":
    filepath = args.out
    if not filepath.endswith(".fa"):
        filepath += ".fa"
else:
    filepath = args.vcf
    if filepath.endswith(".vcf"):
        filepath = filepath[:-4] + ".fa"
    elif filepath.endswith(".vcf.gz"):
        filepath = filepath[:-7] + ".fa"
    else:
        print("Error: unexpected VCF file extension.")
        exit(1)

with open(filepath, "w") as output:
    for column in SNP:
        if column == "POS":
            for i in SNP.index:
                seq = replace_char(seq, SNP["POS"][i] - args.start, "_")
        else:
            strain_seq = seq
            snp_count = 0
            seq_count = 0

            while seq_count < len(strain_seq):
                if strain_seq[seq_count] == "_":
                    snp = SNP[column][snp_count]
                    if snp == "*":
                        strain_seq = replace_char(strain_seq, seq_count, "")
                        seq_count -= 1
                    else:
                        strain_seq = replace_char(strain_seq, seq_count, snp)
                    snp_count += 1
                seq_count += 1

            strain = SeqRecord(Seq(strain_seq), id=column, description="")
            output.write(strain.format("fasta"))
            output.write("\n")
