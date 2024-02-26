from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC, abstractmethod


class FastQFilter:
    def __init__(
        self, input_file, output_file, min_length=0, min_quality=0, min_gc=0, max_gc=100
    ):
        self.input_file = input_file
        self.output_file = output_file
        self.min_length = min_length
        self.min_quality = min_quality
        self.min_gc = min_gc
        self.max_gc = max_gc

    def filter_fastq(self):
        with open(self.output_file, "w") as output_handle:
            for record in SeqIO.parse(self.input_file, "fastq"):
                if self._passes_filter(record):
                    SeqIO.write(record, output_handle, "fastq")

    def _passes_filter(self, record):
        if len(record.seq) < self.min_length:
            return False

        if min(record.letter_annotations["phred_quality"]) < self.min_quality:
            return False

        gc_content = gc_fraction(record.seq)
        if gc_content < self.min_gc or gc_content > self.max_gc:
            return False

        return True


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


class NucleicAcidSequence(BiologicalSequence):
    complement_dict = {"A": "", "T": "", "G": "", "C": ""}

    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def check_alphabet(self):
        valid_alphabet = set("ATGC")
        return set(self.sequence) <= valid_alphabet

    def complement(self):
        complemented_sequence = [
            self.complement_dict.get(base, base) for base in self.sequence
        ]
        return self._create_instance("".join(complemented_sequence))

    def _create_instance(self, sequence):
        return self.__class__(sequence)

    def gc_content(self):
        gc_count = sum(base in "GCgc" for base in self.sequence)
        return gc_count / len(self.sequence)


class DNASequence(NucleicAcidSequence):
    complement_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        transcription_rules = {"A": "U", "T": "A", "G": "C", "C": "G"}
        transcribed_sequence = "".join(
            transcription_rules.get(base, base) for base in self.sequence
        )
        return RNASequence(transcribed_sequence)


class RNASequence(NucleicAcidSequence):
    complement_dict = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def codons(self):
        return [self.sequence[i : i + 3] for i in range(0, len(self.sequence), 3)]


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def check_alphabet(self):
        valid_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        return set(self.sequence) <= valid_alphabet

    def get_molecular_weight(self):
        molecular_weights = {
            "A": 89.09,
            "C": 121.16,
            "D": 133.10,
            "E": 147.13,
            "F": 165.19,
            "G": 75.07,
            "H": 155.16,
            "I": 131.18,
            "K": 146.19,
            "L": 131.18,
            "M": 149.21,
            "N": 132.12,
            "P": 115.13,
            "Q": 146.15,
            "R": 174.20,
            "S": 105.09,
            "T": 119.12,
            "V": 117.15,
            "W": 204.23,
            "Y": 181.19,
        }
        return sum(molecular_weights[aa] for aa in self.sequence)
