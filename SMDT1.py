from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Union
import numpy as np
import os
import string

COMP_BASES_DNA: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't',
                        'g': 'c', 'c': 'g', 't': 'a'}
COMP_BASES_RNA: dict = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'a': 'u',
                        'g': 'c', 'c': 'g', 'u': 'a'}
AA_MASSES: dict = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147, 'Z': 147,
                   'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105,
                   'T': 119, 'W': 204, 'Y': 181, 'V': 117}

class BiologicalSequence:
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return str(self.seq)

    def __repr__(self):
        return f"{self.seq}"

    def __getitem__(self, key):
        if isinstance(key, int):
            item = str(self).__getitem__(key)
            return item
        elif isinstance(key, slice):
            # Get the start, stop, and step from the slice
            piece = [self[i] for i in range(*key.indices(len(self)))]
            return ''.join(piece)
    def check_alphabet(self):
        return set(self.seq) <= set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')

class NucleicAcidSequence(BiologicalSequence):

    def complement(self):
        if 'T' in set(self.seq.upper()):
            comp_seq = [COMP_BASES_DNA[nuc] for nuc in list(self.seq)]
        else:
            comp_seq = [COMP_BASES_RNA[nuc] for nuc in list(self.seq)]
        return ''.join(comp_seq)

    def gc_content(self):
        gc_score = (100 * (self.seq.count('G') + self.seq.count('C')) / len(self.seq))
        return gc_score

class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        rna_seq = self.seq.replace('t', 'u').replace('T', 'U')
        return rna_seq

class RNASequence(NucleicAcidSequence):
    pass

class AminoAcidSequence(BiologicalSequence):
    """

    Calculates the mass (Da) of a protein based on its amino acids sequence.
    Takes a string of amino acids, returns the molecular weight in Da.
    Amino acids in the string should be indicated as one-letter symbols.

    """
    def calculate_protein_mass(self) -> float:
        aa_seq = self.seq.upper()
        mass = 0
        for amino_acid in aa_seq:
            mass += AA_MASSES[amino_acid]
        return mass


def sort_by_quality(input_path: str, quality_threshold: int = 0):
    """
        Returns a list of indexes of sequences that have a quality_score within the given range().
        Allows to sort out unsuitable sequences based on their length.

        Arguments:
        -input_path: str of the path to a fastq file
        -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                        0.2**32 by default

    """
    with open(input_path) as handle:
        quality_ids = [rec.id for rec in SeqIO.parse(handle, "fastq") if
                       (np.mean(rec.letter_annotations["phred_quality"]) >= quality_threshold)]

    return quality_ids


def sort_by_length(input_path: str, length_bounds: Union[int, tuple[int]] = (0.2 ** 32)):
    """

        Returns a list of indexes of sequences that have length within the given range().
        Allows to sort out unsuitable sequences based on their length.

        Arguments:
        -input_path: str of the path to a fastq file
        -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                        0.2**32 by default

    """
    with open(input_path) as handle:
        len_and_ids = [(len(rec), rec.id) for rec in SeqIO.parse(handle, "fastq")]
        record_index = SeqIO.index(input_path, 'fastq')
        if type(length_bounds) is tuple:
            length_ids = [i for (length, i) in len_and_ids if (length_bounds[0] < length < length_bounds[1])]
        else:
            length_ids = [i for (length, i) in len_and_ids if (length > length_bounds)]
        del len_and_ids
        return length_ids


def sort_by_gc(input_path: str, gc_bounds: Union[int, tuple[int]] = (0, 100)):
    """

    Returns a list of indexes of sequences that have gc percentage within of the given range (determined by gc_bounds argument).
    Allows to sort out unsuitable sequences based on the percentage of (G + C) in the nucleotide sequence.

    Arguments:
    -input_path: path to a fastq file
    -qc_bounds: int or tuple[int]: int = the maximum allowed gc percentage, tuple = allowed range of gc percentage
                (0,100) by default

    """
    with open(input_path) as handle:
        gc_values_and_ids = [(100 * gc_fraction(rec.seq), rec.id) for rec in SeqIO.parse(handle, "fastq")]
        if type(gc_bounds) is tuple:
            gc_ids = [i for (gc, i) in gc_values_and_ids if (gc_bounds[0] < gc < gc_bounds[1])]
        else:
            gc_ids = [i for (gc, i) in gc_values_and_ids if (gc > gc_bounds)]
        return gc_ids


def filter_fastq(input_path: str, output_path: str = 'input_filename', length_bounds: Union[int, tuple[int]] = (0.2 ** 32),
                 gc_bounds: Union[int, tuple[int]] = (0, 100), quality_threshold: int = 0):
    """

    Returns a file with sequences satisfying given criteria: length, gc proportion, quality threshold.

    Arguments:
    -input_path: path to a fastq file
    -output_filename: str with an output file name
                      by default is the same as an input file name
    -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                     0.2**32 by default
    -gc_bounds: int or tuple[int]: int indicates = maximum allowed gc percentage, tuple =  allowed range of percentage
                (0,100) by default
    -quality_threshold: int (phred33 quality scores)
                        0 by default


    """
    quality_ids = sort_by_quality(input_path, quality_threshold)
    length_ids = sort_by_length(input_path, length_bounds)
    gc_ids = sort_by_gc(input_path, gc_bounds)
    ids = list((set(quality_ids) & set(length_ids)) & set(gc_ids))
    records = SeqIO.index(input_path, 'fastq')
    if output_path == 'input_filename':
        output_path = os.path.basename(input_path)
    with open(output_path, 'w') as out_handle:
           SeqIO.write([records[i] for i in records if i in ids], out_handle, "fastq")


