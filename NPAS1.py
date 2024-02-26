from abc import ABC, abstractmethod
from random import choice

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def filter_fastq(input_path: str, output_filename: str = None, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0):
    """
    Filter a FASTQ file based on specified criteria and save the filtered sequences to a new FASTQ file.

    Args:
        input_path (str): The path to the input FASTQ file.
        output_filename (str, optional): The name for the output FASTQ file. If not provided, the input file name is used.
        gc_bounds (tuple, optional): GC content filter bounds (default is (0, 100)).
        length_bounds (tuple, optional): Sequence length filter bounds (default is (0, 2**32)).
        quality_threshold (int, optional): Quality score threshold (default is 0).

    Returns:
        None: This function does not return a value but saves the filtered sequences to a new FASTQ file.
    """
    filtered_seqs = []

    for record in SeqIO.parse(input_path, "fastq"):
        gc_percent = gc_fraction(str(record.seq)) * 100
        seq_len = len(record)
        mean_offset = sum(record.letter_annotations["phred_quality"]) / len(record)

        if gc_bounds[0] <= gc_percent <= gc_bounds[1] and \
                length_bounds[0] <= seq_len <= length_bounds[1] and \
                mean_offset >= quality_threshold:
            filtered_seqs.append(record)

    if output_filename is None:
        output_filename = input_path.split("/")[-1]
    else:
        output_filename = output_filename + ".fastq"

    output_path = output_filename
    SeqIO.write(filtered_seqs, output_path, "fastq")

class BiologicalSequence(ABC):
    """
        Abstract base class representing a biological sequence.

        Defines common operations for biological sequences, such as DNA, RNA, or protein sequences.

        Attributes:
            sequence (str): The sequence of characters representing the biological sequence.

        Methods:
            __len__(): Returns the length of the biological sequence.
            __getitem__(index): Gets the character at the specified index or returns a subsequence.
            __str__(): Returns a string representation of the biological sequence.
            is_valid_alphabet(): Checks if all characters in the sequence belong to a valid alphabet.
        """
    def __init__(self, sequence):
        self.sequence = sequence
        if not self.alphabet_checking():
            raise ValueError("Invalid characters in the sequence.")

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    def alphabet_checking(self):
        if not set(self.sequence) <= set(type(self).ALPHABET):
            return False
        return True


class NucleicAcidSequence(BiologicalSequence):
    """
        Represents a nucleic acid sequence, such as DNA or RNA.

        Inherits from BiologicalSequence and adds functionality specific to nucleic acid sequences.

        Attributes:
            sequence (str): The sequence of characters representing the nucleic acid sequence.

        Methods:
            complement(): Returns the complement of the nucleic acid sequence.
            gc_content(as_percentage=False): Calculates the GC content of the nucleic acid sequence.

        Overrides:
            is_valid_alphabet(): Checks if all characters in the sequence belong to a valid nucleic acid alphabet.
        """

    COMPLEMENT_DICT = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c'}

    def __init__(self, sequence):
        super().__init__(sequence)

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def complement(self):
        return ''.join(self.COMPLEMENT_DICT.get(base, base) for base in self.sequence)

    def gc_content(self, as_percentage=False):
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        total_count = len(self.sequence)
        gc_content = gc_count / total_count if total_count > 0 else 0
        return gc_content * 100 if as_percentage else gc_content


class DNASequence(NucleicAcidSequence):
    ALPHABET = set("ATGCatgc")

    TRANSCRIBE_DICT = {
        'T': 'U',
        't': 'u'
    }

    def __init__(self, sequence):
        super().__init__(sequence)

    def transcribe(self):
        transcribed_seq = ''.join(self.TRANSCRIBE_DICT.get(base, base) for base in self.sequence)
        return transcribed_seq


class RNASequence(NucleicAcidSequence):
    ALPHABET = set("AUGCaugc")


class AminoAcidSequence(BiologicalSequence):
    """
        This function takes  aminoacid sequence and translates in to the RNA.
        As most of the aminoacids are coded with several different codons,
        this function will take a random codon of the set for such aminoacids.

        Arguments:
            seq (str): A sequence of RNA molecule

        Output:
            returns sequence of aminoacids
        """

    ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def __init__(self, sequence):
        super().__init__(sequence)

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def translate_to_rna(self):
        AA_CODON_DICT = {
            "G": ["GGA", "GGU", "GGC", "GGG"],
            "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"],
            "S": ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"],
            "E": ["GAA", "GAG"],
            "P": ["CCA", "CCC", "CCG", "CCU"],
            "L": ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"],
            "V": ["GUA", "GUC", "GUG", "GUU"],
            "T": ["ACA", "ACC", "ACG", "ACU"],
            "A": ["GCA", "GCC", "GCG", "GCU"],
            "I": ["AUA", "AUC", "AUU"],
            "F": ["UUC", "UUU"],
            "H": ["CAC", "CAU"],
            "Y": ["UAC", "UAU"],
            "Q": ["CAA", "CAG"],
            "C": ["UGC", "UGU"],
            "N": ["AAC", "AAU"],
            "D": ["GAC", "GAU"],
            "K": ["AAA", "AAG"],
            "M": ["AUG"],
            "W": ["UGG"],
        }
        rna_seq = ""
        for aa in self.sequence:
            codon = choice(AA_CODON_DICT[aa])
            rna_seq += codon
        return rna_seq

