from __future__ import annotations
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from abc import ABC, abstractmethod


def check(args: list[str]):
    """
    Checks list of sequences for correct format.
    Arguments:
        -args (list) - Input format: at least 1 sequence and command, command on last position of list,
        only DNA/RNA seqeunces.
    Return:
        -ValueError - if something does not fit arguments input format.
    """
    if len(args) < 2:
        raise ValueError('Input sequences and command')
    elif args[-1] not in COMMANDS:
        raise ValueError('No such command')
    for seq in args[:-1]:
        if set(seq.lower()) == set('autgc') and set(seq.lower()) != set('atgc') and set(seq.lower()) != set('augc'):
            raise ValueError('Incorrect nucleic acid')


def filter_fastq(
        input_path: str,
        output_filename: str = None,
        gc_bounds: Union[int, tuple[int, int]] = (0, 100),
        length_bounds: Union[int, tuple[int, int]] = (0, 2 ** 32),
        quality_treshholds: int = 0
) -> dict[str, str]:
    """
    Filters fastq files by specifiable parameters.
    Arguments:
        -input_path(str) - path to input file.
        -output_filename (str) - name of output file, input file name will be used if none is given.
        -gc_bounds (int,tuple) - GC interval (in percent) for filtering, default is (0, 100). If input is single int - it will be ceiling (0, n).
        -length_bounds (int,tuple) - length interval for filtering. Works exactly as gc_bounds. Default is (0, 2**32)
        -quality_threshold (int) - The threshold value of average read quality for the filter is 0 by default (phred33 scale).
    Return:
        -output (dict[str,str]) - filtered dictionary, which consists of entities that fulfill entered parameters.
    """

    def average_quality(record):
        return np.mean(record.letter_annotations["phred_quality"])

    if output_filename is None:
        output_filename = input_path

    with open(input_path, "r") as input_handle, open(output_filename, "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "fastq")
        filtered = {}
        for seq in sequences:
            if len(seq) >= length_bounds[0] and len(seq) <= length_bounds[1] \
                    and average_quality(seq) >= quality_threshold \
                    and GC(seq.seq) >= gc_bounds[0] and GC(seq.seq) <= gc_bounds[1]:
                filtered[seq.id] = seq
        SeqIO.write(filtered.values(), output_handle, "fastq")
    return filtered


class BiologicalSequence(ABC):
    def __init__(self, sequence):
        self.sequence = sequence

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    sequence_type: str = None  # Это будет переопределено в подклассах

    def __init__(self, sequence):
        self.sequence = sequence
        if not self.check_alphabet():
            raise ValueError("Invalid sequence")

    def complement(self):
        COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'} if self.sequence_type == 'DNA' else {
            'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(COMPLEMENT_MAP[base] for base in self.sequence)

    def check_alphabet(self):
        valid_bases = set('ATGCU') if self.sequence_type == 'DNA' else set('AUGC')
        return set(self.sequence.upper()) <= valid_bases


class DNASequence(NucleicAcidSequence):
    sequence_type = 'DNA'

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def transcribe(self):
        return (self.complement()).replace('T', 'U')


class RNASequence(NucleicAcidSequence):
    sequence_type = 'RNA'

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence


class AminoAcidSequence(BiologicalSequence):
    sequence_type = 'AminoAcidSequence'

    def __len__(self):
        return len(self.sequence)

    def __init__(self, sequence):
        self.sequence = sequence
        if not self.check_alphabet():
            raise ValueError("Invalid sequence")

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def check_alphabet(self):
        return set(self.sequence.upper()) <= set('ACDEFGHIKLMNPQRSTVWY')

    def calculate_mm(self):
        MOLECULAR_MASS_MAP = {'A': 89.094, 'R': 174.203, 'N': 132.119, 'D': 133.104, 'C': 121.154,
                              'E': 147.131, 'Q': 146.146, 'G': 75.067, 'H': 155.156, 'I': 131.175,
                              'L': 131.175, 'K': 146.189, 'M': 149.208, 'F': 165.192, 'P': 115.132,
                              'S': 105.093, 'T': 119.119, 'W': 204.228, 'Y': 181.191, 'V': 117.148}
        return sum(MOLECULAR_MASS_MAP[aa] for aa in self.sequence)

