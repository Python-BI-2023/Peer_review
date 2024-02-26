import os
from typing import Tuple, Union, Type
from abc import ABC, abstractmethod

import numpy as np
from Bio import SeqIO
from Bio import SeqUtils


class InvalidInput(ValueError):
    """
    Exception raised for invalid input values.

    Inherits from ValueError.

    Attributes:
        message (str): Explanation of the error.
    """
    pass

class NotImplementedError(ValueError):
    """
    Exception raised for call invalid class sample.

    Inherits from ValueError.

    Attributes:
        message (str): Explanation of the error.
    """
    pass


class BiologicalSequence(ABC, str):

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str, alphabet: dict = None, complement_alphabet: dict = None):
        self.sequence = None
        self.alphabet = None
        self.complement_alphabet = None


    def setter(self, sequence: str, alphabet: dict, complement_alphabet: dict = None):
        self.sequence = sequence
        self.alphabet = alphabet
        self.complement_alphabet = complement_alphabet

    def __str__(self):
        return (str(self.sequence))

    def check_alphabet(self) -> bool:
        """
        Check if a sequence is DNA or RNA.

        Returns:
            bool: True if the sequence is DNA or RNA, False otherwise.
        """

        unique_chars = set(self.sequence)
        if not (unique_chars <= self.alphabet):
            return False
        return True  

    def complement(self) -> str:
        """
        Find the complement of a DNA or RNA sequence.

        Returns:
            str: The complemented sequence.
        """
        if self.alphabet == None:
            raise NotImplementedError()
        new_seq = []
        for nucl in self.sequence:
            new_seq.append(self.complement_alphabet.get(nucl))
        return type(self)(''.join(new_seq))

    def gc_content(self):
        return SeqUtils.gc_fraction(self.sequence)


class DNASequence(NucleicAcidSequence):
    _alphabet = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
    _complement_alphabet = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a'}
    transcribe_alphabet = {
        "a": "a", "A": "A",
        "t": "u", "T": "U",
        "u": "t", "U": "T",
        "g": "g", "G": "G",
        "c": "c", "C": "C"
    }

    def __init__(self, sequence: str, alphabet: dict = _alphabet,
                 complement_alphabet: dict = _complement_alphabet):
        super().setter(sequence, alphabet, complement_alphabet)

    def transcribe(self):
        return RNASequence(''.join([self.transcribe_alphabet[i] for i in self.sequence]))


class RNASequence(NucleicAcidSequence):
    _alphabet = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
    _complement_alphabet = {'A': 'U', 'a': 'u', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'U': 'A', 'u': 'a'}

    def __init__(self, sequence: str, alphabet: dict = _alphabet,
                 complement_alphabet: dict = _complement_alphabet):
        super().setter(sequence, alphabet, complement_alphabet)


class AminoAcidSequence(BiologicalSequence):
    alphabet = {'V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P',
                'C', 'v', 'i', 'l', 'e', 'q', 'd', 'n', 'h', 'w', 'f', 'y', 'r', 'k', 's', 't', 'm', 'a', 'g', 'p', 'c'}
    average_weights = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
        'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
        'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
        'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }
    three_letter_alphabet = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
                             'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
                             'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr',
                             'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'}

    def __init__(self, sequence: str):
        self.sequence = sequence

    def __str__(self):
        return (str(self.sequence))
    
    
    def check_alphabet(self) -> bool:
        """
        Check if a sequence is DNA or RNA.

        Returns:
            bool: True if the sequence is DNA or RNA, False otherwise.
        """

        unique_chars = set(self.sequence)
        if not (unique_chars <= self.alphabet):
            return False
        return True  
    
    
    def aa_average_weight(self, weight: str = 'average') -> float:
        """
        Calculate the amino acids weight in a protein sequence.

        Args:
            seq (str): The amino acid sequence to calculate the weight for.
            weight (str, optional): The type of weight to use, either 'average' or 'monoisotopic'. Default is 'average'.

        Returns:
            float: The calculated weight of the amino acid sequence.
        """

        final_weight = 0
        for aa in self.sequence.upper():
            final_weight += self.average_weights[aa]
        return round(final_weight, 3)


    def one_to_three_letter_code(self) -> str:
        """
        This function converts a protein sequence from one-letter amino acid code to three-letter code.

        Args:
            sequence (str): The input protein sequence in one-letter code.

        Returns:
            str: The converted protein sequence in three-letter code.
        """

        three_letter_code = [self.three_letter_alphabet.get(aa.upper()) for aa in self.sequence]
        return '-'.join(three_letter_code)


def filter_dna(input_path: str, output_filename: str = '', gc_bounds: Union[Tuple[int, int], int] = (0, 100),
               length_bounds: Union[Tuple[int, int], int] = (0, 2 ** 32), quality_threshold: int = 0) -> None:
    """
    Filter and process a dictionary of FASTQ sequences based on specified criteria.

    Args:
        input_path (str): Path to the input FASTQ file. Please write your path with directory etc.
                          You can use os.path.join(dir_name, file_name)
        output_filename (str): Name of the output FASTQ file. By default, it is the same as the input name.
        gc_bounds (tuple or int): GC content filtering bounds.
            If a tuple, it represents the lower and upper bounds (inclusive) for GC content as percentages.
            If an int, it represents the upper bound for GC content as a percentage.
        length_bounds (tuple or int, optional): Length filtering bounds.
            If a tuple, it represents the lower and upper bounds (inclusive) for sequence length.
            If an int, it represents the upper bound for sequence length.
            Default is (0, 2**32).
        quality_threshold (int, optional): Quality threshold for filtering sequences based on average quality.
            Sequences with an average quality below this threshold will be discarded.
            Default is 0 (phred33 scale).

    Returns:
        filtered_seqs (dict): A file containing filtered FASTQ sequences.
            Key: Sequence name (string).
            Value: Tuple of two strings (sequence, comment, quality).

    Example:
        filtered_seqs = filtr_dna(seqs, gc_bounds=(20, 80), length_bounds=50, quality_threshold=30)
        - here gc_bounds interval will be [20, 80] and length_bounds will be [0, 50].
          quality_threshold will be 30 <= sequence score

    """

    good_reads = []
    records = SeqIO.parse(input_path, "fastq")
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    if not isinstance(gc_bounds, tuple) or not isinstance(length_bounds, tuple):
        raise InvalidInput()
    for record in records:
        if (np.mean(record.letter_annotations["phred_quality"]) >= quality_threshold) and \
                is_in_gc_bounds(gc_bounds, SeqUtils.gc_fraction(record) * 100) and \
                is_in_length_bounds(length_bounds, len(record.seq)):
            good_reads.append(record)

    if output_filename == '':
        output_filename = os.path.basename(input_path)
    if not os.path.exists('fastq_filtrator_results'):
        os.mkdir('fastq_filtrator_results')
    output_filename = os.path.join('fastq_filtrator_results', output_filename)
    SeqIO.write(good_reads, handle=output_filename, format="fastq")

    return None

def is_in_gc_bounds(bounds: tuple, gc_content: float) -> bool:
    """
    Check if the GC content of a DNA sequence falls within the specified bounds.

    Args:
        bounds (tuple): A tuple specifying the lower and upper bounds for GC content.
        gc_content (float): GC content of input DNA.

    Returns:
        bool: True if the GC content is within the bounds, False otherwise.
    """

    upper_bound = max(bounds[0], bounds[1])
    lower_bound = min(bounds[0], bounds[1])
    if upper_bound <= 0 or lower_bound < 0:
        raise InvalidInput("Invalid gc_bounds. Each value must be greater than zero.")
    if upper_bound == lower_bound:
        upper_bound += 1
    return lower_bound <= gc_content <= upper_bound


def is_in_length_bounds(bounds: tuple, dna_length: int) -> bool:
    """
    Check if the length of a DNA sequence falls within the specified bounds.

    Args:
        bounds (tuple): A tuple specifying the lower and upper bounds for sequence length.
        dna_length (int): The length of input DNA sequence.

    Returns:
        bool: True if the sequence length is within the bounds, False otherwise.
    """

    lower_bound = min(bounds[0], bounds[1])
    upper_bound = max(bounds[0], bounds[1])
    if upper_bound < 0 or lower_bound < 0:
        raise ValueError("Invalid length_bounds. Each value must be greater than zero.")
    if upper_bound == lower_bound:
        upper_bound += 1
    return lower_bound <= dna_length <= upper_bound
