from __future__ import annotations
import os
import sys
from Bio import SeqIO, SeqUtils, SeqRecord
from Bio.SeqUtils import gc_fraction


class BiologicalSequence:
    """
    Represents a generic biological sequence.
    Attributes:
    - sequence (str): The biological sequence.
    Methods:
    - __len__(): Returns the length of the sequence.
    - __getitem__(index: int): Returns the character at the specified index.
    - __str__(): Returns a string representation of the sequence.
    - is_valid_alphabet(): Checks if the sequence uses a valid alphabet.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return f"{self.sequence}"

    def is_valid_alphabet(self) -> bool:
        """
        Checks if the sequence uses a valid alphabet.
        """
        try:
            unique_chars = set(self.sequence)
            return unique_chars <= self.ALPHABET
        except AttributeError as e:
            raise NotImplementedError(
                "Is valid alphabet method not implemented for this class."
            ) from e


class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence.
    Methods:
    - complement(): Returns the complemented sequence.
    - gc_content(percentage: bool = True) -> int | float: Calculates the GC content of the sequence.
    Raises:
    - NotImplementedError: If complement method is not implemented for this class.
    """

    def complement(self) -> DNASequence | RNASequence:
        try:
            complemented_sequence = "".join(
                [self.COMPLEMENT_MAP[base] for base in self.sequence]
            )
            return self.__class__(complemented_sequence)
        except AttributeError as e:
            raise NotImplementedError(
                "Complement method not implemented for this class."
            ) from e

    def gc_content(self, percentage: bool = True) -> int | float:
        """
        Calculates the GC content of the sequence.
        Args:
        - percentage (bool): If True, returns the result as a percentage.
        Returns:
        - int | float: GC content value.
        """
        gc_symbols = set("GCgc")
        gc_count = sum(1 for nucleotide in self.sequence if nucleotide in gc_symbols)
        if percentage:
            return (gc_count / self.sequence.__len__()) * 100
        else:
            return gc_count / self.sequence.__len__()


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.
    Methods:
    - transcribe(): Returns the transcribed RNA sequence.
    """

    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.ALPHABET = set("ATGCatgc")
        self.COMPLEMENT_MAP = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "a": "t",
            "c": "g",
            "g": "c",
            "t": "a",
        }

    def transcribe(self) -> RNASequence:
        """
        Returns the transcribed RNA sequence.
        Returns:
        - RNASequence: Transcribed RNA sequence.
        """
        transcribed_sequence = self.sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.
    """

    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.ALPHABET = set("AUGCaugc")
        self.COMPLEMENT_MAP = {
            "A": "U",
            "C": "G",
            "G": "C",
            "U": "A",
            "a": "u",
            "c": "g",
            "g": "c",
            "u": "a",
        }


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.
    """

    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.ALPHABET = set("GgLlYySsEeQqDdNnFfAaKkRrHhCcVvPpWwIiMmTt")

    def count_aa(self) -> dict:
        """
        Counts the number of given or all amino acids in a protein sequence.
        Arguments:
        - seq (str): sequence to count amino acids
        - aminoacids (str): which amino acids to count in sequence
        Return:
        - dict: a dictionary with amino acids and its count
        """

        aa_dict_count = {}
        for aa in set(self.sequence):
            aa_dict_count[aa] = self.sequence.count(aa)
        return aa_dict_count


def filter_fastq(
    path_to_seqs: str,
    output_file_name: str = None,
    gc_bounds: tuple | int = (0, 100),
    length_bounds: tuple | int = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filters FASTQ sequences based on the GC-content, length and quality parameters.

    Args:
    - path_to_seqs (str): the path to the FASTQ file to be filtered.
    - output_file_name (str): the name of the file where the filtered FASTQ sequences will be saved.
    - gc_bounds (tuple, int): GC content range (in percentages) for filtering. If you pass a single number
    to the argument, it is assumed to be an upper bound.
    - length_bounds (tuple, int): length range for filtering. If you pass a single number to the argument,
    it is assumed to be an upper bound.
    - quality_threshold (int): threshold value for average read quality filtering

    Returns:
    - None: the function doesn't return a value but writes the filtered FASTQ to the output file.
    """
    if output_file_name is None:
        output_file_name = os.path.basename(path_to_seqs)

    records = SeqIO.parse(path_to_seqs, "fastq")

    filtered_records = []

    for record in records:
        gc_content = gc_fraction(record.seq)
        if isinstance(gc_bounds, tuple):
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue
        elif gc_content >= gc_bounds:
            continue
        sequence_length = len(record.seq)
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= sequence_length <= length_bounds[1]):
                continue
        elif sequence_length >= length_bounds:
            continue
        mean_quality = sum(record.letter_annotations["phred_quality"]) / sequence_length
        if mean_quality <= quality_threshold:
            continue
        filtered_records.append(record)

    with open(output_file_name, "w") as output_handle:
        SeqIO.write(filtered_records, output_handle, "fastq")
