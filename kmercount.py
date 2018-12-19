#!/lustre/apps/anaconda/3/bin/python
import sys
import argparse
import multiprocessing as mp
from Bio import SeqIO
from functools import partial
import re
from datetime import datetime
import os
import math
import numpy as np
from scipy import stats


__author__ = "Guillermo Luque, Alejandro Reyes"
__copyright__ = "Copyright 2017, kmercount"
__credits__ = ["BCEM"]
__license__ = "MIT"
__maintainer__ = "Guillermo Luque"
__email__ = "guillermo.luque.ds@gmail.com, a.reyes@uniandes.edu.co"
__version__ = '1.0.5'


def ck(args, x):
    if x in args.__dict__:
        return args.__dict__[x]
    else:
        return None


def _dna_iupac(s):
    return {
        'R': '[AG]',
        'Y': '[CT]',
        'S': '[GC]',
        'W': '[AT]',
        'K': '[GT]',
        'M': '[AC]',
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'V': '[ACG]',
        'N': '[ACGT]'
    }.get(s, s)


def regex_from(seq, disable_overlapping=False):
    pattern = ""
    for s in seq:
        pattern += _dna_iupac(s)
    if not disable_overlapping:
        return "(?=(%s))" % pattern
    else:
        return "(%s)" % pattern


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


class Counter:
    def __init__(self, args):
        self.fasta = ck(args, 'fasta')
        self.output = ck(args, 'output')
        self.kmer = ck(args, 'kmer')
        self.disable_overlapping = ck(args, 'disable_overlapping')
        self.reverse_complement = ck(args, 'reverse_complement')
        self.debug = ck(args, 'debug')
        self.split_size = ck(args, 'split_size')

        self.records = []
        self.step = 0

        if not self.fasta:
            raise Exception("No FASTA file set")

        if not self.output:
            raise Exception("No output file set")

        if not self.kmer:
            raise Exception("No K-mer given")

        # Processing regex
        kmer_lst = self.kmer.upper().split(",")
        self.kmer_regex_lst = []
        for kmer in kmer_lst:
            ele = (kmer, re.compile(regex_from(kmer, self.disable_overlapping)))
            self.kmer_regex_lst.append(ele)

        self.prefix = os.path.basename(self.fasta)
        self.prefix = os.path.splitext(self.prefix)[0]

        self.np = mp.cpu_count()
        self._debug("Available CPUs: %s\n" % self.np)

    def exec(self):
        self._debug("KMERCOUNT [run at %s]\n" % datetime.now())
        self._debug("Flags:\n")
        self._debug("--fasta: %s\n" % self.fasta)
        self._debug("--output: %s\n" % self.output)
        self._debug("--kmer: %s\n" % self.kmer)
        self._debug("--disable_overlapping: %s\n" % self.disable_overlapping)
        self._debug("--reverse_complement: %s\n" % self.reverse_complement)
        self._debug("--debug: %s\n" % self.debug)

        file_size = round(os.path.getsize(self.fasta) / 1024 / 1024 / 1024)
        self._debug("File size: %s GB\n" % file_size)

        if file_size == 0 and self.split_size > 0:
            self._debug("Your file size is less than 1 GB. The split_size flag (-s) will be ignored\n")
            self.split_size = 0

        total = 0
        with open(self.fasta) as handle:
            for _ in SeqIO.parse(handle, 'fasta'):
                total += 1

        if total == 0:
            raise Exception("There are no sequences in file")

        self._debug("Total of sequences: %s\n" % total)

        # Save header in output file
        with open(self.output, "w") as fh:
            fh.write("record.id\tkmer\tn.matches\tn.fw\tn.re\t"
                     "fw.min\tfw.max\tfw.mean\tfw.std\tre.min\tre.max\tre.mean\tre.std\n")

        cont = 0
        if self.split_size > 0:
            sequence_size = file_size / total
            chunk_size = math.ceil(self.split_size / sequence_size)
            self._debug("Sequences size (aprox.): %s GB\n" % sequence_size)
            self._debug("Sequences per chunk: %s\n" % chunk_size)

            record_iter = SeqIO.parse(open(self.fasta), "fasta")
            for i, batch in enumerate(batch_iterator(record_iter, chunk_size)):
                cont += 1
                filename = "%s_%i.fx" % (self.prefix, i + 1)
                with open(filename, "w") as handle:
                    x = SeqIO.write(batch, handle, "fasta")
                self._debug("Wrote %i records to %s\n" % (x, filename))

            self._debug("Total of files: %s\n" % cont)

            for i in range(cont):
                filename = "%s_%i.fx" % (self.prefix, i + 1)
                self._process_file(filename, i + 1)
                os.remove(filename)

            # consolidate output
            with open(self.output, "a") as fh:
                for i in range(cont):
                    filename = "%s.%i" % (self.output, i + 1)
                    with open(filename) as infile:
                        for line in infile:
                            fh.write(line)

            for i in range(cont):
                filename = "%s.%i" % (self.output, i + 1)
                os.remove(filename)
        else:
            self._process_file(self.fasta)

        self._debug('END OF PROCESSING\n')
        return None

    def _process_file(self, filename, _id=0):
        self._debug("Processing filename: %s\n" % filename)
        self.records = []
        with open(filename) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                self.records.append(record)

        total = len(self.records)
        self._debug("FASTA records: %s\n" % total)

        if total == 0:
            return None

        self.step = total // self.np
        offsets = [x for x in range(0, total, self.step)]

        self._debug("Steps: %s\n" % self.step)
        self._debug("Offsets: %s\n" % offsets)

        manager = mp.Manager()
        queue = manager.Queue()
        fn = partial(self._search_kmer, _id)

        pool = mp.Pool()
        pool.map(fn, offsets)
        pool.close()
        pool.join()

        self._debug('>>> %s processed\n' % filename)
        return None

    def _search_kmer(self, _id, offset):
        self._debug("Search Range: %s-%s\n" % (offset, offset + self.step))

        for record in [x for x in self.records[offset:offset+self.step]]:
            self._debug("SEQ: %s\t%s\n" % (record.id, record.seq[:10]))
            for kmer_str, kmer in self.kmer_regex_lst:
                self._debug("SEQ: %s\t%s\tkmer_str: %s\tPattern: %s\n" % (record.id, record.seq[:10], kmer_str, kmer))

                matches = []
                _fw_pos = []
                for _m in kmer.finditer(str(record.seq)):
                    matches.append(_m.group(1))
                    _fw_pos.append(_m.start(1))

                _fw_min = _fw_max = _fw_mean = _fw_std = 0
                _re_min = _re_max = _re_mean = _re_std = 0
                _n_re = 0

                _n_fw = len(matches)

                if _n_fw > 0:
                    try:
                        _fw_lens = np.array([len(x) for x in matches])
                        _fw_dist = np.diff(_fw_pos) - _fw_lens[1:]
                        _, (_fw_min, _fw_max), _fw_mean, _fw_std, _, _ = stats.describe(_fw_dist)
                        _fw_std = np.sqrt(_fw_std)
                    except Exception as _ex:
                        print(_ex)

                    self._debug("Forward matches (number): %s\n" % len(matches))
                    self._debug("Forward matches: %s...\n" % matches[:5])

                if self.reverse_complement:
                    m = []
                    _re_pos = []
                    for _m in kmer.finditer(str(record.seq.reverse_complement())):
                        m.append(_m.group(1))
                        _re_pos.append(_m.start(1))

                    _n_re = len(m)
                    if _n_re > 0:
                        try:
                            _re_lens = np.array([len(x) for x in m])
                            _re_dist = np.diff(_re_pos) - _re_lens[1:]
                            _, (_re_min, _re_max), _re_mean, _re_std, _, _ = stats.describe(_re_dist)
                            _re_std = np.sqrt(_re_std)
                        except Exception as _ex:
                            print(_ex)

                        matches.extend(m)

                        self._debug("Reverse matches (number): %s\n" % len(m))
                        self._debug("Reverse matches: %s...\n" % m[:5])

                outtxt = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%s\t%s\t%.2f\t%.2f\n" % (
                    record.id, kmer_str, len(matches), _n_fw, _n_re,
                    _fw_min, _fw_max, _fw_mean, _fw_std,
                    _re_min, _re_max, _re_mean, _re_std)

                self._write_results(outtxt, _id)

        self._debug("--> end offset %s\n" % offset)
        return None

    def _write_results(self, outtxt, _id=0):
        if _id == 0:
            filename = self.output
        else:
            filename = "%s.%i" % (self.output, _id)

        self._debug("Writing %s\n" % filename)

        with open(filename, "a") as fh:
            fh.write(outtxt)

    def _debug(self, logtxt):
        if self.debug:
            with open(self.debug, "a") as fh:
                fh.write(logtxt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Return the number of occurrences of a given k-mer per sequence")
    parser.add_argument("-f", "--fasta", default=None, metavar="FILE",
                        help="FASTA file containing DNA sequences.")
    parser.add_argument("-k", "--kmer", type=str, help="K-mer in IUPAC code")
    parser.add_argument("-o", "--output", default=None, metavar="FILE",
                        help="Output file.")
    parser.add_argument("-r", "--reverse_complement", help="Search in the reverse complement also",
                        action="store_true")
    parser.add_argument("-d", "--disable_overlapping", help="Search with no overlapping",
                        action="store_true")
    parser.add_argument("-g", "--debug", default=None, metavar="FILE",
                        help="Debugging output file.")
    parser.add_argument("-s", "--split_size", default=0, type=int,
                        help="Size in GB of each chunk in which the FASTA file will be splitted."
                             "Zero means no split at all.")
    args = parser.parse_args()

    try:
        counter = Counter(args)
        counter.exec()
    except Exception as e:
        print(e)
        sys.exit(-1)
