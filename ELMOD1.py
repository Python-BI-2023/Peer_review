from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from typing import List, Dict, Union, Tuple
from abc import ABC, abstractmethod
import os


class FastQFilter:
    def __init__(self, input_path: str, output_filename: str, gc_bounds: Union[int, Tuple[int, int]] = (0, 100),
                 length_bounds: Union[int, Tuple[int, int]] = (0, 2**32), quality_threshold: int = 0):
        """
        Initialize FastQFilter instance.

        Parameters:
        - input_path (str): Path to the input FastQ file.
        - output_filename (str): Name for the output FastQ file.
        - gc_bounds (Union[int, Tuple[int, int]]): GC content bounds for filtering.
        - length_bounds (Union[int, Tuple[int, int]]): Length bounds for filtering.
        - quality_threshold (int): Quality threshold for filtering.
        """
        self.input_path = input_path
        self.output_filename = output_filename
        self.gc_bounds = gc_bounds
        self.length_bounds = length_bounds
        self.quality_threshold = quality_threshold

    def filter_fastq(self) -> None:
        """
        Read, filter, and write FastQ records based on specified criteria.
        """
        records = self.read_fastq()
        filtered_records = self.apply_filters(records)
        self.write_fastq(filtered_records)

    def read_fastq(self) -> List[SeqRecord]:
        """
        Read FastQ file and return a list of SeqRecord objects.
        """
        records = list(SeqIO.parse(self.input_path, "fastq"))
        return records

    def apply_filters(self, records: List[SeqRecord]) -> List[SeqRecord]:
        """
        Filter SeqRecord objects based on specified criteria.

        Parameters:
        - records (List[SeqRecord]): List of SeqRecord objects to filter.

        Returns:
        - List[SeqRecord]: Filtered list of SeqRecord objects.
        """
        filtered_records = []

        for record in records:
            gc_content = gc_fraction(record.seq)
            length = len(record.seq)
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])

            gc_pass = self.check_bounds(gc_content, self.gc_bounds)
            length_pass = self.check_bounds(length, self.length_bounds)
            quality_pass = avg_quality >= self.quality_threshold

            if gc_pass and length_pass and quality_pass:
                filtered_records.append(record)

        return filtered_records

    def check_bounds(self, value: int, bounds: Union[int, Tuple[int, int]]) -> bool:
        """
        Check if a value is within the specified bounds.

        Parameters:
        - value (int): Value to check.
        - bounds (Union[int, Tuple[int, int]]): Bounds to check against.

        Returns:
        - bool: True if value is within bounds, False otherwise.
        """
        if isinstance(bounds, int):
            return value <= bounds
        else:
            return bounds[0] <= value <= bounds[1]

    def write_fastq(self, records: List[SeqRecord]) -> None:
        """
        Write SeqRecord objects to a new FastQ file.

        Parameters:
        - records (List[SeqRecord]): List of SeqRecord objects to write.
        """
        output_dir = "fastq_filtrator_results"
        output_path = f"{output_dir}/{self.output_filename}.fastq"

        os.makedirs(output_dir, exist_ok=True)

        SeqIO.write(records, output_path, "fastq")
        print(f"Filtered sequences written to {output_path}")

class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    """

    def __init__(self, sequence: str):
        """
        Initialize a BiologicalSequence instance.

        Parameters:
        - sequence (str): The biological sequence.

        The sequence is automatically converted to uppercase.
        """
        self.sequence = sequence.upper()
        self.validate_alphabet()

    def __len__(self) -> int:
        """
        Return the length of the sequence.

        Returns:
        - int: Length of the sequence.
        """
        return len(self.sequence)
    
    def __getitem__(self, index: int) -> str:
        """
        Get the character at the specified index in the sequence.

        Parameters:
        - index (int): The index to retrieve.

        Returns:
        - str: The character at the specified index.
        """
        return self.sequence[index]
    
    def __str__(self) -> str:
        """
        Return the string representation of the sequence.

        Returns:
        - str: String representation of the sequence.
        """
        return self.sequence
    
    @abstractmethod
    def is_valid_alphabet(self) -> bool:
        """
        Check if the sequence contains a valid alphabet.

        Returns:
        - bool: True if the alphabet is valid, False otherwise.
        """
        pass

    def validate_alphabet(self) -> None:
        """
        Validate the alphabet of the sequence.

        Raises:
        - ValueError: If the alphabet is not valid.
        """
        if not self.is_valid_alphabet():
            raise ValueError("Invalid alphabet for {}: {}".format(self.__class__.__name__, str(self)))

class NucleicAcidSequence(BiologicalSequence):
    """
    Class representing nucleic acid sequences.
    """
    def complement(self) -> 'NucleicAcidSequence':
        """
        Return the complement of the sequence.

        Returns:
        - NucleicAcidSequence: Complement of the sequence.
        """
        self.validate_alphabet()
        return self.__class__(''.join([self.complement_base(base) for base in self.sequence]))
    
    def gc_content(self) -> float:
        """
        Calculate the GC content of the sequence.

        Returns:
        - float: GC content of the sequence.
        """
        self.validate_alphabet()
        gc_count = sum(1 for base in self.sequence if base in {'G', 'C'})
        return gc_count / len(self)

class DNASequence(NucleicAcidSequence):
    """
    Class representing DNA sequences.
    """
    def transcribe(self) -> 'RNASequence':
        """
        Transcribe the DNA sequence into RNA.

        Returns:
        - RNASequence: Transcribed RNA sequence.
        """
        self.validate_alphabet()
        return RNASequence(''.join(['U' if base == 'T' else base for base in self.sequence]))
    
    def complement_base(self, base: str) -> str:
        """
        Return the complement of a DNA base.

        Parameters:
        - base (str): The DNA base.

        Returns:
        - str: Complement of the DNA base.
        """
        if base == 'A':
            return 'T'
        elif base == 'T':
            return 'A'
        elif base == 'C':
            return 'G'
        elif base == 'G':
            return 'C'
    
    def is_valid_alphabet(self) -> bool:
        """
        Check if the DNA sequence contains a valid alphabet.

        Returns:
        - bool: True if the alphabet is valid, False otherwise.
        """
        return set(self.sequence) <= {'A', 'T', 'C', 'G'}

class RNASequence(NucleicAcidSequence):
    """
    Class representing RNA sequences.
    """
    def is_valid_alphabet(self) -> bool:
        """
        Check if the RNA sequence contains a valid alphabet.

        Returns:
        - bool: True if the alphabet is valid, False otherwise.
        """
        return set(self.sequence) <= {'A', 'U', 'C', 'G'}
    
    def complement_base(self, base: str) -> str:
        """
        Return the complement of an RNA base.

        Parameters:
        - base (str): The RNA base.

        Returns:
        - str: Complement of the RNA base.
        """
        if base == 'A':
            return 'U'
        elif base == 'U':
            return 'A'
        elif base == 'C':
            return 'G'
        elif base == 'G':
            return 'C'

class AminoAcidSequence(BiologicalSequence):
    """
    Class representing amino acid sequences.
    """
    amino_brutto = {
        "A": (3, 7, 1, 2, 0),
        "R": (6, 14, 4, 2, 0),
        "N": (4, 8, 2, 3, 0),
        "D": (4, 7, 1, 4, 0),
        "V": (5, 11, 1, 2, 0),
        "H": (6, 9, 3, 2, 0),
        "G": (2, 5, 1, 2, 0),
        "Q": (5, 10, 2, 3, 0),
        "E": (5, 9, 1, 4, 0),
        "I": (6, 13, 1, 2, 0),
        "L": (6, 13, 1, 2, 0),
        "K": (6, 14, 2, 2, 0),
        "M": (5, 11, 1, 2, 1),
        "P": (5, 9, 1, 2, 0),
        "S": (3, 7, 1, 3, 0),
        "Y": (9, 11, 1, 3, 0),
        "T": (4, 9, 11, 1, 3, 0),
        "W": (11, 12, 2, 2, 0),
        "F": (9, 11, 1, 2, 0),
        "C": (3, 7, 1, 2, 1),
    }

    def is_valid_alphabet(self) -> bool:
        """
        Check if the amino acid sequence contains a valid alphabet.

        Returns:
        - bool: True if the alphabet is valid, False otherwise.
        """
        return set(self.sequence) <= {'A', 'R', 'N', 'D', 'V', 'H', 'G', 'Q', 'E', 'I', 'L', 'K', 'M',
                                      'P', 'S', 'Y', 'T', 'W', 'F', 'C'}

    def brutto_count(self) -> List[Dict[str, int]]:
        """
        Calculate the Brutto formula values for each amino acid in the sequence.

        Returns:
        - List[Dict[str, int]]: List of dictionaries containing Brutto formula values for each amino acid.
        """
        self.validate_alphabet()
        elements = ["C", "H", "N", "O", "S"]
        result = []
        for letter in self.sequence:
            brutto_values = self.amino_brutto.get(letter, (0, 0, 0, 0, 0))
            result.append(dict(zip(elements, brutto_values)))
        return result
    