"""
This file was used to generate those large strings from cpg_Stats.py
"""

import random
import itertools


def condense(my_dict):
    '''
    because we don't know which reads are on the same strand, we can combine kmer couns that are reverse complements of each other
    :param my_dict: kmer count dict
    :return: dict
    '''
    def rc(dna):
        MAPPING = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return ''.join([MAPPING[x] for x in dna][::-1])

    new_dict = set()
    for k in sorted(my_dict):
        if k not in new_dict and rc(k) not in new_dict:
            new_dict.add(k)
        elif k not in new_dict and rc(k) in new_dict:
                pass
        else:
            print(k)  # This should never happen
    return list(new_dict)


dna = ["A", "C", "G", "T"]
for i in [2, 3, 6]:
    kmers = list(itertools.product(dna, repeat=i))
    kmers = sorted([''.join(x) for x in kmers])
    print(', '.join(["sixmers.get('{}', 0)".format(x) for x in sorted(condense(kmers))]))
