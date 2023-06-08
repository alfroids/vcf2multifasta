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
    p.add_argument(
        "--des",
        help="Default description for all mulsti FASTA entries.",
        required=False,
        type=str,
        nargs="+",
        default="",
    )
    p.add_argument(
        "--out",
        help="Path to the output multi FASTA file.",
        required=False,
        type=str,
        nargs="?",
        default="",
    )
    p.add_argument(
        "--gaps",
        "-g",
        help="Insert naive gaps on the output sequences.",
        required=False,
        action="store_true",
    )
    p.add_argument(
        "--verbose", "-v", help="Verbose mode.", required=False, action="store_true"
    )
    return p


parser = initialize_parser()
args = parser.parse_args()


# Read VCF and create SNPs DataFrame ##########################################
if args.verbose:
    print("Reading and filtering VCF file.")

vcf_reader = vcf.Reader(open(args.vcf, "rb"))


if args.end > -1:
    vcf_reader = vcf_reader.fetch(args.chr, start=args.start, end=args.end)
else:
    vcf_reader = vcf_reader.fetch(args.chr, start=args.start)

cols = ["POS"] + vcf_reader.samples
snps = {c: [] for c in cols}

if args.verbose:
    print("Creating and processing SNPs DataFrame.")

for record in vcf_reader:
    snps["POS"].append(record.POS)
    ref = record.REF

    for sample in record.samples:
        if sample.gt_bases is None:
            snps[sample.sample].append(record.REF)
        elif args.gaps and sample.gt_bases == "*":
            snps[sample.sample].append("")
        else:
            snps[sample.sample].append(sample.gt_bases)

SNP = pd.DataFrame(snps)

if args.gaps:
    for i in SNP.index:
        M = max(
            [len(s) for s in SNP.loc[i, SNP.columns != "POS"].values.flatten().tolist()]
        )
        m = min(
            [len(s) for s in SNP.loc[i, SNP.columns != "POS"].values.flatten().tolist()]
        )

        if M != m:
            for column in SNP:
                if column != "POS":
                    if len(SNP.at[i, column]) < M:
                        SNP.at[i, column] += (M - len(SNP.at[i, column])) * "-"


# Read reference FASTA and extract the sequence at the specified range ########
if args.verbose:
    print("Reading FASTA file and extracting sequence.")

fasta_dict = SeqIO.index(args.ref, "fasta")
chrom = fasta_dict[args.chr]
locus = chrom[args.start : args.end + 1]
seq = str(locus.seq)


# Substitute SNPs into reference sequence and output multi-FASTA file #########
if args.verbose:
    print("Computing sequences and outputing to FASTA file.")

if args.out != "":
    filepath = args.out
    if not filepath.endswith(".fasta"):
        filepath += ".fasta"
else:
    filepath = args.vcf
    if filepath.endswith(".vcf"):
        filepath = filepath[:-4] + ".fasta"
    elif filepath.endswith(".vcf.gz"):
        filepath = filepath[:-7] + ".fasta"
    else:
        print("Error: unexpected VCF file extension.")
        exit(1)

count = 0
desc = " ".join(args.des).strip('"')
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

            strain = SeqRecord(Seq(strain_seq), id=column, description=desc)
            output.write(strain.format("fasta"))
            output.write("\n")
            count += 1

if args.verbose:
    print("Successfully generated FASTA file with {} entries.".format(count))
