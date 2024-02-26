import os
from typing import List, Union
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def calculate_gc_content(seq: Seq) -> float:
    """
    Calculate the GC content of a nucleotide sequence.

    Arguments:
        seq (Seq): Biological sequence.

    Returns:
        float: GC content as a percentage.
    """
    return GC(seq)


def calculate_average_quality(seq: Seq) -> float:
    """
    Calculate the average quality of a nucleotide sequence.

    Arguments:
        seq (Seq): Biological sequence.

    Return:
        float: Average sequence quality according to the phred33 scale.
    """
    return sum(seq.letter_annotations["phred_quality"]) / len(seq)   # обращаемся к методу letter_annotations из Bio.SeqRecord


def filter_fastq(
    input_path: str,
    output_filename: str = '',
    gc_bounds: Union[tuple, int, float] = (0, 100),
    length_bounds: Union[tuple, int] = (0, 2**32),
    quality_threshold: int = 0
) -> str:
    """
    Filter FASTQ records based on specified criteria and save them to a new file.

    Arguments:
        input_path (str): Path to the FASTQ file.
        output_filename (str): Optional, file name to save the result.
        gc_bounds (Union[tuple, int, float]): GC composition interval (percents) to filter, default is (0, 100).
        length_bounds (Union[tuple, int]): Sequence length interval to filter, default is (0, 2**32).
        quality_threshold (int): Threshold value of average read quality (phred33) to filter, default is 0.

    Return:
        str: Path for filtered file.
    """
    records = list(SeqIO.parse(input_path, "fastq"))  # 1 элемент списка = 1 последовательноть формата fastq

    filtered_records = []

    # Обрабатываем границы ГЦ-состава
    if isinstance(gc_bounds, (float, int)):
        gc_min = 0
        gc_max = gc_bounds
    else:
        gc_min = gc_bounds[0]
        gc_max = gc_bounds[1]

    # Обрабатываем границы по длине
    if isinstance(length_bounds, int):
        length_min = 0
        length_max = length_bounds
    else:
        length_min = length_bounds[0]
        length_max = length_bounds[1]

    for record in records:
        gc_content = calculate_gc_content(record.seq)
        length = len(record)
        avg_quality = calculate_average_quality(record)

        if (
            gc_min <= gc_content <= gc_max
            and length_min <= length <= length_max
            and avg_quality >= quality_threshold
        ):
            filtered_records.append(record)

    return save_records_to_fastq(filtered_records, output_filename, input_path)


def save_records_to_fastq(filtered_records: List[SeqIO.SeqRecord], output_filename: str, input_path: str) -> str:
    """
    Save filtered records to a new FASTQ file.

    Arguments:
        filtered_records (List[SeqIO.SeqRecord]): List of filtered SeqRecords.
        output_filename (str): Name of the file to save.
        input_path (str): Path to the original FASTQ file.

    Returns:
        str: Path to the saved file.
    """
    if output_filename == '':  # если имя не передано, используем имя входного файла
        output_filename = os.path.basename(input_path)
    if not output_filename.endswith(".fastq"):
        output_filename += ".fastq"
    output_dir = 'fastq_filtrator_results'  # задаём имя папки для сохранения результата
    os.makedirs(output_dir, exist_ok=True)  # создаем папку, если она не существует
    output_path = os.path.join(output_dir, output_filename)  # задаем путь, по которому сохранится файл

    SeqIO.write(filtered_records, output_path, "fastq")

    return output_path


class BiologicalSequence(ABC):

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index) -> str:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def alphabet_is_valid(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return self.sequence

    def alphabet_is_valid(self) -> bool:
        dna_alphabet = set("ATGCatgc")
        rna_alphabet = set("AGCUagcu")
        if isinstance(self, DNASequence):
            return set(self.sequence).issubset(dna_alphabet)
        if isinstance(self, RNASequence):
            return set(self.sequence).issubset(rna_alphabet)
        return False

    def complement(self) -> 'NucleicAcidSequence':
        complement_pairs = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a'
        }
        complement_sequence = ''
        for base in self.sequence:
            complement_sequence += complement_pairs.get(base, base)

        return type(self)(complement_sequence)

    def gc_content(self) -> float:
        gc_count = 0
        for base in self.sequence:
            if base in {'G', 'g', 'C', 'c'}:
                gc_count += 1
        total_bases = len(self.sequence)
        gc_percentage = gc_count / total_bases
        return gc_percentage


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.sequence = sequence

    def transcribe(self) -> 'RNASequence':
        transcribe_sequence = ''
        for base in self.sequence:
            if base == 'T':
                transcribe_sequence += 'U'
            if base == 't':
                transcribe_sequence += 'u'
            else:
                transcribe_sequence += base
        return RNASequence(transcribe_sequence)


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)
        self.sequence = sequence


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return self.sequence

    def alphabet_is_valid(self) -> bool:
        amino_acid_alphabet = set("ARNDCHGQEILKMPSYTWFV")
        return set(self.sequence).issubset(amino_acid_alphabet)

    def calculate_molecular_weight(self) -> float:
        amino_acid_weights = {
            'G': 57.051, 'A': 71.078, 'S': 87.077, 'P': 97.115, 'V': 99.131,
            'T': 101.104, 'C': 103.143, 'I': 113.158, 'L': 113.158, 'N': 114.103,
            'D': 115.087, 'Q': 128.129, 'K': 128.172, 'E': 129.114, 'M': 131.196,
            'H': 137.139, 'F': 147.174, 'R': 156.186, 'Y': 163.173, 'W': 186.210
        }
        molecular_weight = 0.0
        for aa in self.sequence:
            molecular_weight += amino_acid_weights[aa]

        return molecular_weight
