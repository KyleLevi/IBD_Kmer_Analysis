import sys
import os
import argparse
import re


def fastq_seqs(fname):
    i = 0
    with open(fname, 'r') as infile:
        for line in infile:
            if i%4 == 1:
                yield line
                i = 1
            i += 1


def regex_motifs(fname):
    """

    :param fname:
    :param kwargs:
    :return:
    """
    dict_motifs = {1: '(?=(CG[^CG]{1,5}CG))', 2: "(?=(CG[^CG]{6,10}CG))", 3: "(?=(CG[^CG]{11,15}CG))", 4: "(?=(CG[^CG]{16,}CG))"}

    motif_counts = {}

    for seq in fastq_seqs(fname):
        if "N" in seq:
            continue
        seq = seq.replace('\n', '')

        for k,v in dict_motifs.items():
            found_motifs = re.findall(v, seq)
            motif_counts[k] = motif_counts.get(k, 0) + len(found_motifs)

    return motif_counts.get(1, 0), motif_counts.get(2, 0), motif_counts.get(1, 3), motif_counts.get(1, 4),

def raw_motifs(fname):
    pairs = []
    for seq in fastq_seqs(fname):
        prevC = False
        in_cpg = False
        start = None
        for i, char in enumerate(seq):
            if in_cpg and char == 'G' and prevC:
                pairs.append(seq[start-1:i+1])
                start = i
            elif char == 'G' and prevC:
                in_cpg = True
                start = i
            if char == 'C':
                prevC = True
            else:
                prevC = False
            if start and i-start > 30:
                start = None
                in_cpg = False

    return pairs


# print(raw_motifs('testfiles/small.fastq'))




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input File', required=True)
    parser.add_argument('-o', '--output', help='output File')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(1)

    if not args.output:
        args.output = sys.stdout

    with open(args.output, 'w') as outfile:
        for seq in fastq_seqs(args.input):
            prevC = False
            in_cpg = False
            start = None
            for i, char in enumerate(seq):
                if in_cpg and char == 'G' and prevC:
                    outfile.write(seq[start - 1:i + 1] + '\n')
                    start = i
                elif char == 'G' and prevC:
                    in_cpg = True
                    start = i
                if char == 'C':
                    prevC = True
                else:
                    prevC = False
                if start and i - start > 30:
                    start = None
                    in_cpg = False
