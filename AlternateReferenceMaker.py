#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Juliane Hartlich
# date : March 2018
# developed at Leibnz Institute DSMZ-German Collection of Migcroorganisms an Cell Cultures
# in cooperation with Technical University Braunschweig

# ---packages----------------------
import argparse, copy

# ---functions---------------------
def command():
    Parser = argparse.ArgumentParser(
        prog="AlternateReferenceMaker.py",
        usage="%(prog)s FASTA VCF",
        description="%(prog)s is used to alternate a FASTA sequence with variants called in a VCF file",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    Parser.add_argument(
        "INPUT",
        nargs=2,
        help="name of input FASTA file and name of VCF file, \nFASTA file must be first, VCF file must be second, \nplease add directory if file is not in current working directory",
    )
    Args = Parser.parse_args()
    return Args


# ---main program------------------
ARGs = command()

# generating new FASTA file
FANAME = ARGs.INPUT[0].split(".")[0]
if FANAME.find("/") != -1:
    NAME = FANAME.split("/")
    NAME = "_".join(NAME[1 : len(NAME) - 1])
else:
    NAME = FANAME
NewFAName = "%s_consensus.fasta" % NAME
NewFA = open(NewFAName, "w")

## accessing reference sequence
Rfile = open(ARGs.INPUT[0], "r")
Ref = Rfile.read()
Rfile.close()
Ref = Ref.replace("\r", "")
Ref = Ref.split(">")[1:]  # preparation of file if multi-FASTA

## accessing VCF file
Vfile = open(ARGs.INPUT[1], "r")
Vcf = Vfile.read()
Vfile.close()
Vcf = Vcf.replace("\r", "")
Vcf_rows = Vcf.split("\n")

# preparing reference sequence
for FaSe in Ref:  # if FASTA File contains multiple sequences
    # saving FASTA sequence label
    Label = FaSe.split("\n")[0]
    Label = Label.split(" ")[0]  # split at whitespaces
    # converting sequence into one string, adding '#' so sequence starts with 1 and not 0
    RefSeq = "#" + "".join(FaSe.split("\n")[1:])

    VCFinfo = []

    for row in Vcf_rows:  # here 'row' is a string
        if row.find("#") == -1 and len(row) != 0:  # ignoring the VCF file header
            row = row.split("\t")  # here 'row' is a list
            if row[0] == Label.replace(">", ""):  # needed to handle multi-FASTA files
                Pos = int(row[1])
                Ori = row[3]
                Alt = row[4]

                VCFinfo[len(VCFinfo) :] = [[Pos, Ori, Alt]]

    VCFinfo = sorted(VCFinfo, key=lambda VCFinfo: VCFinfo[0], reverse=True)

    for entry in VCFinfo:
        Pos = entry[0]
        Ori = entry[1]
        Alt = entry[2]

        lala = ""
        for i in range(len(Ori)):
            lala += RefSeq[Pos + i]

        lala = lala.upper()
        Ori = Ori.upper()

        if lala == Ori:
            if Alt.find("D") != -1 or Alt.find(".") != -1:
                EditSeq = RefSeq[:Pos] + RefSeq[Pos + len(Ori) :]
                RefSeq = EditSeq
            else:
                EditSeq = RefSeq[:Pos] + Alt + RefSeq[Pos + len(Ori) :]
                RefSeq = EditSeq
        elif Ori == ".":
            if Alt.find("I") != -1:  # handeling PacBio VCF output
                Alt = Alt.strip("I")
            EditSeq = RefSeq[: Pos + 1] + Alt + RefSeq[Pos + 1 :]
            RefSeq = EditSeq

    RefSeq = RefSeq[1:]

    # writing edited sequence into new file
    NewFA.write(">%s_consensus\n" % Label)
    line = ""
    SeqCop = copy.copy(RefSeq)
    runs = len(SeqCop) // 80
    rest = len(SeqCop) - (runs * 80)
    if rest != 0:
        runs += 1
    for blubb in range(0, runs):
        line = SeqCop[0:80]
        NewFA.write("%s\n" % line)
        SeqCop = SeqCop[80:]


NewFA.close()
