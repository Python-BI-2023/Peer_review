import numpy as np
from abc import ABC, abstractmethod
from typing import (
    Union,
    Dict,
    Tuple,
    Set,
    Self
)  # для импорта Self нужен python >= 3.11 или используйте импорт ниже
# from typing_extensions import Self

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
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


class ExpressedSequences(BiologicalSequence):
    def __init__(self, seq: str, alphabet: Set[str]) -> None:
        self.letters = set(seq)
        self.seq = seq
        self._alphabet = alphabet

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, s):
        if not s:
            raise ValueError("Sequence can't be empty!")
        if not self.check_alphabet():
            raise ValueError("Sequence contains incorrect symbols!")
        self._seq = s

    def __len__(self) -> int:
        return len(self.seq)

    def __getitem__(self, index: int) -> "Self":
        return self.__class__(self.seq[index])

    def __str__(self) -> str:
        return self.seq

    def check_alphabet(self) -> bool:
        return all([x in self._alphabet for x in self.letters]) and not (
            "U" in self.letters and "T" in self.letters
        )


class NucleicAcidSequence(ExpressedSequences):
    def __init__(self, seq: str) -> None:
        self._alphabet = {"A", "T", "G", "C", "U"}
        super().__init__(seq, self._alphabet)
        self._complement_table = None

    def complement(self) -> "Self":
        if self._complement_table is None:
            raise NotImplementedError
        complemented = []
        for nucleotide in self.seq:
            complemented.append(self._complement_table[nucleotide])
        return self.__class__("".join(complemented))

    def gc_content(self) -> int:
        return (self.seq.count("G") + self.seq.count("C")) / len(self.seq) * 100


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str) -> None:
        super().__init__(seq)
        self._complement_table = {"T": "A", "A": "T", "C": "G", "G": "C"}
        self._transcribe_table = {
            "C": "C",
            "G": "G",
            "A": "A",
            "T": "U",
        }

    def transcribe(self) -> "RNASequence":
        transcribed = []
        for nucleotide in self.seq:
            transcribed.append(self._transcribe_table[nucleotide])
        return RNASequence("".join(transcribed))


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str):
        super().__init__(seq)
        self._complement_table = {"A": "U", "U": "A", "C": "G", "G": "C"}


class AminoAcidSequence(ExpressedSequences):
    def __init__(self, seq: str) -> None:
        self._alphabet = {
            "A",
            "G",
            "D",
            "L",
            "N",
            "P",
            "C",
            "Y",
            "S",
            "I",
            "H",
            "W",
            "E",
            "F",
            "R",
            "T",
            "V",
            "K",
            "M",
            "Q",
        }
        super().__init__(seq, self._alphabet)
        self._positive_aa = {"R", "K", "H"}
        self._negative_aa = {"D", "E"}

    def define_charge(self) -> Dict[str, int]:
        positive_count = 0
        negative_count = 0
        neutral_count = 0
        for aa in self.seq:
            if aa in self._positive_aa:
                positive_count += 1
            elif aa in self._negative_aa:
                negative_count += 1
            else:
                neutral_count += 1
        result = {
            "Positive": positive_count,
            "Negative": negative_count,
            "Neutral": neutral_count,
        }
        return result


def filter_fastq(
    input_path: str,
    output_filename: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (
        0,
        100,
    ),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filters appropriate sequences.

    Arguments
    ----------
    input_path: str
        Path to fastq file
    output_filename: str
        Name of filtered fastq file
    gc_bounds: Union[int, float, Tuple[Union[int,float], Union[int,float]]
        Filter parameter (in percents) of gc composition. Tuple associated with lowest and highest levels of gc content in  sequences, also it can get on input int or float associated with highest level of gc content, lowest level will be set to 0.
    length_bounds: Union[int, Tuple[int, int]]
        Filter parameter of sequences length. Logic the same as in gc_bounds.
    quality_threshold: int
        Filter parameter of mean sequence quality. Sets lowest mean quality of sequences.

    Saves filtered sequences in a fastq file named input_path/output_filename.
    """
    records = SeqIO.parse(input_path, "fastq")
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)
    records_filtered = []
    for rec in records:
        seq_length = len(rec)
        if seq_length < length_bounds[0] or seq_length > length_bounds[1]:
            continue
        gc_count = gc_fraction(rec.seq) * 100
        if gc_count < gc_bounds[0] or gc_count > gc_bounds[1]:
            continue
        mean_quality = np.mean(rec.letter_annotations["phred_quality"])
        if mean_quality < quality_threshold:
            continue
        records_filtered.append(rec)
    SeqIO.write(records_filtered, output_filename, "fasta")