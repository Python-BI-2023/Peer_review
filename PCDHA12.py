import os
from abc import ABC, abstractmethod

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC


# 2.4
def filter_length(record: SeqRecord, length_bounds: tuple) -> bool:
    return length_bounds[0] <= len(record.seq) <= length_bounds[1]


def filter_gc(record: SeqRecord, gc_bounds: tuple) -> bool:
    return gc_bounds[0] <= GC(record.seq) <= gc_bounds[1]


def filter_quality(record: SeqRecord, quality_threshold: int) -> bool:
    avg_quality = np.mean(record.letter_annotations["phred_quality"])
    return avg_quality >= quality_threshold


def filter_fastq(input_path: str,
                 gc_bounds: tuple = (0, 100),
                 length_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> None:
    """ Reads a FASTQ file, filters sequences based on GC content, sequence
    length and quality threshold, and writes the filtered sequences to
    a new file """
    path, filename = os.path.split(input_path)
    name, ext = os.path.splitext(filename)
    output_path = path + "/" + f"{name}_filtered{ext}"
    input_seq_iterator = SeqIO.parse(input_path, 'fastq')
    filtered_seq_iterator = (record for record in input_seq_iterator
                             if filter_length(record, length_bounds)
                             and filter_gc(record, gc_bounds)
                             and filter_quality(record, quality_threshold))
    SeqIO.write(filtered_seq_iterator, output_path, "fastq")


# 2.5
class BiologicalSequence(ABC, str):
    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        raise NotImplementedError('An instance of this class cannot be created')

    def check_alphabet(self) -> bool:
        """Checks if the provided nucleic acid sequence contains
        only acceptable letters"""
        sequence_chars = set(self.sequence)
        return sequence_chars.issubset(type(self).ALPHABET)

    def complement(self):
        """Complements DNA or RNA sequence."""
        complement_seq = ''.join(type(self).COMPLEMENT_DICT.get(base)
                                 for base in self.sequence)
        return type(self)(complement_seq)

    def get_gc(self) -> float:
        """Calculates GC content for nucleic acid sequence"""
        return GC(self.sequence)


class RNASequence(NucleicAcidSequence):
    ALPHABET = ('A', 'U', 'G', 'C', 'N')
    COMPLEMENT_DICT = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                       'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}

    def __init__(self, sequence):
        self.sequence = sequence


class DNASequence(NucleicAcidSequence):
    ALPHABET = ('A', 'T', 'G', 'C', 'N')
    COMPLEMENT_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                       'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}

    def __init__(self, sequence):
        self.sequence = sequence

    def transcribe(self) -> RNASequence:
        """Transcribes DNA sequence into RNA sequence.
        Returns RMASequence object"""
        return RNASequence(self.sequence.replace('T', 'U').replace('t', 'u'))


class AminoAcidSequence(BiologicalSequence):
    ONE_LETTER_ALPHABET = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
    THREE_LETTER_ALPHABET = ('Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His',
                             'Ile', 'Lys', 'Leu', 'Met', 'Pro', 'Gln', 'Arg',
                             'Ser', 'Thr', 'Val', 'Trp', 'Tyr')

    def __init__(self, sequence: str):
        self.sequence = sequence

    def check_alphabet(self) -> bool:
        """Checks if the provided amino acid sequence contains only
        acceptable letters in one-letter or three-letter code"""
        sequence_chars = set(self.sequence)
        return (sequence_chars.issubset(self.ONE_LETTER_ALPHABET)
                or sequence_chars.issubset(self.THREE_LETTER_ALPHABET))

    def get_molecular_weight(self) -> float:
        """ calculate molecular weight for one-letter amino acid sequence"""
        WEIGHT_DICT = {
            'G': 57.051, 'A': 71.078, 'S': 87.077, 'P': 97.115, 'V': 99.131,
            'T': 101.104, 'C': 103.143, 'I': 113.158, 'L': 113.158, 'N': 114.103,
            'D': 115.087, 'Q': 128.129, 'K': 128.172, 'E': 129.114, 'M': 131.196,
            'H': 137.139, 'F': 147.174, 'R': 156.186, 'Y': 163.173, 'W': 186.210
        }
        terminal_h_oh_weight = 18.02
        weight = (terminal_h_oh_weight +
                  sum(WEIGHT_DICT[aa] for aa in self.sequence if aa in WEIGHT_DICT))
        return weight
