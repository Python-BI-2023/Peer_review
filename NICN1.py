from typing import List, Union
import os
from Bio import SeqIO
from Bio import SeqUtils
from abc import ABC, abstractmethod


class IncorrectNucleotideError(ValueError):
    """Exception raised for invalid nucleotide in sequence."""

    pass


class IncorrectAminoacidEncodingError(ValueError):
    """Exception raised for invalid amino acid encoding."""

    pass


class IncorrectAminoacidError(ValueError):
    """Exception raised for invalid amino acid in sequence."""

    pass


class BiologicalSequence(ABC):
    """Abstract base class representing a biological sequence."""

    @abstractmethod
    def __len__(self) -> int:
        """Return the length of the sequence."""
        pass

    @abstractmethod
    def __getitem__(self, index: int) -> str:
        """Return the item at the given index in the sequence."""
        pass

    @abstractmethod
    def __str__(self) -> str:
        """Return the string representation of the sequence."""
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid elements."""
        pass


class NucleicAcidSequence(BiologicalSequence, ABC):
    """Class representing a nucleic acid sequence."""

    def __init__(self, sequence: str):
        self.sequence = sequence
        self._nucleotide_pairs = {}

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return str(self.sequence)

    def check_alphabet(self) -> bool:
        nucleotides = self._nucleotide_pairs.keys()

        for nucleotide in set(self.sequence):
            if nucleotide not in nucleotides:
                raise IncorrectNucleotideError(
                    f"Invalid nucleotide found: {nucleotide}"
                )
        return True

    def complement(self) -> "NucleicAcidSequence":
        """Return the complement of the sequence."""
        self.check_alphabet()

        complement_seq = "".join(
            self._nucleotide_pairs[nucleotide] for nucleotide in self.sequence
        )
        return type(self)(complement_seq)

    def get_gc_content(self) -> Union[float, int]:
        """Return the GC content of the sequence."""
        self.check_alphabet()

        gc_content = SeqUtils.GC(self.sequence.upper())
        return gc_content


class DNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self._nucleotide_pairs = {
            "A": "T",
            "a": "t",
            "G": "C",
            "g": "c",
            "T": "A",
            "t": "a",
            "C": "G",
            "c": "g",
        }

    def transcribe(self) -> "RNASequence":
        """Return the RNA transcript of the sequence."""
        transcribed_seq = self.sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    """Class representing an RNA sequence."""

    def __init__(self, sequence: str):
        super().__init__(sequence)
        self._nucleotide_pairs = {
            "A": "U",
            "a": "u",
            "G": "C",
            "g": "c",
            "U": "A",
            "u": "a",
            "C": "G",
            "c": "g",
        }


class AminoAcidSequence(BiologicalSequence, ABC):
    """Class representing an amino acid sequence."""

    def __init__(self, sequence, encoding):
        self.sequence = sequence.upper()
        self.encoding = encoding
        self._residue_names = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLN": "Q",
            "GLU": "E",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V",
        }
        self._residue_mass = {
            "A": 89,
            "R": 174,
            "N": 132,
            "D": 133,
            "C": 121,
            "Q": 146,
            "E": 147,
            "G": 75,
            "H": 155,
            "I": 131,
            "L": 131,
            "K": 146,
            "M": 149,
            "F": 165,
            "P": 115,
            "S": 105,
            "T": 119,
            "W": 204,
            "Y": 181,
            "V": 117,
        }

    def _reformat_based_on_encoding(self) -> Union[str, List[str]]:
        """Reformat the sequence based on the encoding."""
        if self.encoding == 1:
            return self.sequence
        elif self.encoding == 3:
            amino_acid_list = [
                self.sequence[letter: letter + 3]
                for letter in range(0, len(self.sequence), 3)
            ]
            return amino_acid_list
        else:
            raise IncorrectAminoacidEncodingError(
                f"{self.encoding}-letter encoding is unavailable. "
                f"Please, use 1 or 3 letter encoding "
            )

    def _make_one_letter(self) -> str:
        """Convert the sequence to one-letter code."""
        if self.encoding == 3:
            one_letter_seq = str()
            for aminoacid in self._reformat_based_on_encoding():
                one_letter_seq += self._residue_names[aminoacid]
            return one_letter_seq
        return self.sequence

    def __len__(self) -> int:
        return len(self._reformat_based_on_encoding())

    def __getitem__(self, index: int) -> Union[str, List[str]]:
        if self.encoding == 3:
            return "".join(self._reformat_based_on_encoding()[index])
        return self.sequence[index]

    def __str__(self) -> str:
        return str(self.sequence)

    def check_alphabet(self) -> bool:
        """Check if the sequence contains valid amino acids."""
        aminoacids = set(self._reformat_based_on_encoding())
        for aminoacid in aminoacids:
            if (
                aminoacid not in self._residue_names.keys()
                and aminoacid not in self._residue_names.values()
            ):
                raise IncorrectAminoacidError(
                    f"{aminoacid} is not a supported amino acid!"
                )
        return True

    def get_molecular_mass(self) -> Union[float, int]:
        """Calculate the molecular mass of the sequence."""
        sequence = self._make_one_letter()
        mass = 0
        for aminoacid in sequence:
            mass += self._residue_mass[aminoacid]

        return mass


def filter_fastq(
    input_path: str,
    output_filename=None,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    """
    Filter DNA sequences based on the GC-content, length and sequencing quality (phred33) from FASTQ file to a new
    FASTQ file and stores it in fastaq_filtered_results folder in the same directory.

    :param input_path: path to the sequences in fastq format
    :param output_filename: name for output fastq file (str);
    if not given the output file name will be filtered_*input file name*.fastq
    :param gc_bounds: given threshold for GC-content (tuple/int/float)
    :param length_bounds: given threshold for length (tuple/int)
    :param quality_threshold: given threshold for quality (int)
    """

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    output_dir = os.path.join(os.path.dirname(input_path), "fastq_filtrator_results")
    os.makedirs(output_dir, exist_ok=True)

    if output_filename:
        output_file_path = os.path.join(output_dir, output_filename)
    else:
        output_file_path = os.path.join(
            output_dir, f"filtered_{os.path.basename(input_path)}"
        )

    with open(output_file_path, "w") as f:
        for record in SeqIO.parse(input_path, "fastq"):
            gc_test = gc_bounds[0] < SeqUtils.GC(record.seq) < gc_bounds[1]
            len_test = length_bounds[0] < len(record.seq) < length_bounds[1]
            quality_test = (
                sum(record.letter_annotations["phred_quality"]) / len(record.seq)
                >= quality_threshold
            )
            if gc_test and len_test and quality_test:
                SeqIO.write(record, f, "fastq")
