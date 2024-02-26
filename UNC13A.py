from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Sequence


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
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class BiologicalSequenceInitializer(BiologicalSequence):
    def __init__(self):
        self._seq_set = None
        self._seq = None

    # инициализируем алфавит
    @property
    def seq_set(self) -> set:
        return self._seq_set

    @seq_set.setter
    def seq_set(self, seq_set: set):
        if not isinstance(seq_set, set):
            raise TypeError("Must be a set.")
        self._seq_set = seq_set

    # инициализируем последовательность
    @property
    def seq(self) -> str | Sequence[str]:
        return self._seq

    @seq.setter
    def seq(self, seq: str | Sequence[str]):
        if not (
            isinstance(seq, Sequence) and all(isinstance(elem, str) for elem in seq)
        ):
            raise TypeError("Must be a string or sequence of string objects.")
        if not set(seq).issubset(self._seq_set):
            raise TypeError(
                f"Sequence does not match alphabet for {repr(self)}. Please, use only: {self._seq_set}."
            )
        self._seq = seq

    # реализуем заданные в абстрактном классе методы
    def __len__(self):
        return self.seq.__len__()

    def __getitem__(self, index: int | slice):
        return self.seq.__getitem__(index)

    def __str__(self):
        seq = self.seq
        return "".join(seq) if not seq is None else ""

    def __repr__(self):
        return f"BiologicalSequenceInitializer({str(self.seq)})"

    def check_alphabet(self) -> bool:
        return set(self._seq).issubset(self._seq_set)


class NucleicAcidSequence(BiologicalSequenceInitializer):
    def __init__(self):
        super().__init__()  # https://stackoverflow.com/a/64504667
        self._complement_rule = None

    # инициализируем правило комплементарности
    @property
    def complement_rule(self) -> dict:
        return self._complement_rule

    @complement_rule.setter
    def complement_rule(self, complement_rule: dict):
        if not isinstance(complement_rule, dict):
            raise TypeError("Must be a dict.")
        self._complement_rule = complement_rule

    def __repr__(self):
        return f"NucleicAcidSequence({str(self.seq)})"

    def complement(self) -> RNASequence | DNASequence:
        self_class, self_seq_class = type(self), type(self.seq)
        complement_seq = str(self).translate(str.maketrans(self.complement_rule))
        return self_class(self_seq_class(complement_seq))

    def gc_content(self) -> float:
        return sum(map(str(self).count, ("G", "C", "g", "c"))) / len(self)


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("AUGCaugc")
        self.complement_rule = dict(zip("AUGCaugc", "UACGuacg"))
        self.seq = seq

    def __repr__(self):
        return f"RNASequence({str(self.seq)})"


class DNASequence(NucleicAcidSequence):
    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("ATGCatgc")
        self.complement_rule = dict(zip("ATGCatgc", "TACGtacg"))
        self.seq = seq

    def __repr__(self):
        return f"DNASequence({self.seq})"

    def transcribe(self) -> RNASequence:
        return RNASequence(str(self).translate(str.maketrans(("Tt"), ("Uu"))))


class AminoAcidSequence(BiologicalSequenceInitializer):

    aa_uniprot_content = {
        "A": 9.03,
        "R": 5.84,
        "N": 3.79,
        "D": 5.47,
        "C": 1.29,
        "Q": 3.80,
        "E": 6.24,
        "G": 7.27,
        "H": 2.22,
        "I": 5.53,
        "L": 9.85,
        "K": 4.93,
        "M": 2.33,
        "F": 3.88,
        "P": 4.99,
        "S": 6.82,
        "T": 5.55,
        "W": 1.30,
        "Y": 2.88,
        "V": 6.86,
    }

    def __init__(self, seq: str | Sequence[str]):
        super().__init__()
        self.seq_set = set("ARNDCQEGHILKMFPSTWYV")
        self.seq = seq

    def __repr__(self):
        return f"AminoAcidSequence({str(self.seq)})"

    def content_check(self) -> dict:
        "Returns aminoacids content of the protein"
        seq_content = dict.fromkeys(AminoAcidSequence.aa_uniprot_content.keys(), 0)
        for aacd in str(self.seq).upper():
            seq_content[aacd] = seq_content[aacd] + 1

        seq_length = len(self.seq)
        for aacd, occurence in seq_content.items():
            seq_content[aacd] = 100 * occurence / seq_length

        return seq_content


# def run_fastq_processor(*, input_path: str, output_filename: str = None, output_path: str = 'fastq_filtrator_results', gc_bounds: Union[Tuple[int, int], int] = (0, 100), length_bounds: Tuple[int, int] = (0, 2**32), quality_thershold: int = 0):
#     """Filters reads presented in input fasta file using three metrics:
#     - GC-content;
#     - sequence length;
#     - average phred quality.

#     Args:
#         input_file (str): path to input file.
#         output_filename (str): name of output file.
#         output_path (str): name of folder to store output file.
#         gc_bounds (tuple|int): desired interval for GC-content if argument
#             is tuple; upper bound of this interval if argument is integer
#             (lower will be 0).
#         length_bounds (tuple|int): desired interval for sequence length if
#             argument is tuple; upper bound of this interval if argument is
#             integer (lower will be 0).
#         quality_thershold (tuple|int): desired interval for average phred
#             quality if argument is tuple; upper bound of this interval if
#             argument is integer (lower will be 0).

#     Returns:
#         None.
#     """
#     output_filename = fp.process_paths(input_path, output_filename, output_path)
#     fastq_dict, result = fp.process_file(input_path), {}
#     for name, seq in fastq_dict.items():
#         is_seq_valid = fp.check_seq_and_bounds(seq, gc_bounds, length_bounds, quality_thershold)
#         if is_seq_valid:
#             result[name] = seq
#     fp.save_output(result, output_filename)
