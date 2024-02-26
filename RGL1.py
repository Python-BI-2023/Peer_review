from typing import Union
from Bio import SeqIO, SeqUtils
from abc import ABC, abstractmethod


def length_filter(record, lower_length_bound: int, upper_length_bound: int):
    if lower_length_bound <= len(record.seq) <= upper_length_bound:
        return record


def gc_filter(record, lower_gc_bound: float, upper_gc_bound: float):
    if lower_gc_bound <= SeqUtils.GC(record.seq) <= upper_gc_bound:
        return record


def quality_filter(record, quality_threshold: Union[int, float] = 0):
    if (quality_threshold < sum(record.letter_annotations['phred_quality']) /
            len(record.letter_annotations['phred_quality'])):
        return record


def filter_fastq(file_path: str, output_file: str = 'input_fasta_name_filter.fastq',
                 lower_gc_bound: int = 0, upper_gc_bound: int = 30,
                 lower_length_bound: int = 0, upper_length_bound: int = 2**32,
                 quality_threshold: Union[int, float] = 0) -> None:
    records = SeqIO.parse(file_path, 'fastq')
    if output_file == 'input_fasta_name_filter.fastq':
        output_filename = file_path.split('\\')[-1]
        output_file = output_filename.replace('.fastq', '_filter.fastq')
    with open(output_file, 'w') as output:
        for record in records:
            if length_filter(record, lower_length_bound, upper_length_bound) and \
                    gc_filter(record, lower_gc_bound, upper_gc_bound) and \
                    quality_filter(record, quality_threshold):
                SeqIO.write(record, output, 'fastq')


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, item: int):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class SequenceFunction(BiologicalSequence):
    alphabet = ()

    def __init__(self, seq: str) -> None:
        self.seq = seq

    def __len__(self) -> int:
        self.length = len(self.seq)
        return self.length

    def __getitem__(self, item: int) -> str:
        if 0 <= item < len(self.seq):
            return self.seq[item]
        else:
            raise IndexError('Your index is incorrect')

    def __str__(self) -> str:
        return str(self.seq)

    def check_alphabet(self) -> bool:
        return set(self.seq).issubset(set(self.alphabet))


class NucleicAcidSequence(SequenceFunction):
    complement_dict = {}
    def __init__(self, seq: str):
        super().__init__(seq)

    def complement(self) -> str:
        complement_seq = self.seq.translate(str.maketrans(self.complement_dict))
        return complement_seq

    def gc_content(self) -> Union[int, float]:
        gc_count = self.seq.count('C') + self.seq.count('G')
        return gc_count/len(self.seq)*100


class DNASequence(NucleicAcidSequence):
    alphabet = ('A', 'T', 'G', 'C')
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def __init__(self, seq: str):
        super().__init__(seq)

    def check_alphabet(self):
        return super().check_alphabet()

    def complement(self):
        return super().complement()

    def transcribe(self) -> str:
        transcribed_seq = self.seq.translate(str.maketrans('ATGC', 'UACG'))
        return transcribed_seq


class RNASequence(NucleicAcidSequence):
    alphabet = ('A', 'U', 'G', 'C')
    complement_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    def complement(self):
        return super().complement()

    def check_alphabet(self):
        return super().check_alphabet()


class AminoAcidSequence(SequenceFunction):
    alphabet = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
    amino_acid_weights = {
        'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121,
        'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131,
        'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115,
        'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117
    }

    def __init__(self, seq) -> None:
        super().__init__(seq)

    def check_alphabet(self):
        return super().check_alphabet()

    def count_molecular_weight(self, amino_acid_weights: dict) -> int:
        molecular_weight = sum(self.amino_acid_weights.get(aa, 0) for aa in self.seq)
        return molecular_weight
