from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Any

class BiologicalSequence(ABC):
    """Абстрактный базовый класс для биологических последовательностей."""
    
    def __init__(self, sequence: str) -> None:
        if not self.is_valid(sequence):
            raise ValueError("Invalid sequence for the given type.")
        self._sequence = sequence
    
    @abstractmethod
    def is_valid(self, sequence: str) -> bool:
        """Проверяет, является ли последовательность валидной для данного типа."""
        pass

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index: Any) -> Any:
        # Поддержка получения элемента по индексу и срезов
        if isinstance(index, slice):
            return self.__class__(self._sequence[index])
        else:
            return self._sequence[index]

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._sequence})"

class NucleicAcidSequence(BiologicalSequence):
    """Класс для нуклеиновых кислот."""

    @abstractmethod
    def complement(self) -> 'NucleicAcidSequence':
        """Возвращает комплементарную последовательность."""
        pass

    def gc_content(self) -> float:
        """Вычисляет GC-содержание последовательности."""
        gc_count = sum(base in {'G', 'C'} for base in self._sequence)
        return gc_count / len(self._sequence)

class DNASequence(NucleicAcidSequence):
    """Класс для ДНК последовательности."""

    def is_valid(self, sequence: str) -> bool:
        return all(nucleotide in {'A', 'C', 'G', 'T'} for nucleotide in sequence)

    def complement(self) -> 'DNASequence':
        complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return DNASequence("".join(complement_dict.get(base, base) for base in self._sequence))

    def transcribe(self) -> 'RNASequence':
        """Транскрибирует ДНК последовательность в РНК."""
        return RNASequence(self._sequence.replace('T', 'U'))

class RNASequence(NucleicAcidSequence):
    """Класс для РНК последовательности."""

    def is_valid(self, sequence: str) -> bool:
        return all(nucleotide in {'A', 'C', 'G', 'U'} for nucleotide in sequence)

    def complement(self) -> 'RNASequence':
        complement_dict = {"A": "U", "U": "A", "C": "G", "G": "C"}
        return RNASequence("".join(complement_dict.get(base, base) for base in self._sequence))

class AminoAcidSequence(BiologicalSequence):
    """Класс для аминокислотных последовательностей."""

    def is_valid(self, sequence: str) -> bool:
        valid_amino_acids = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                             'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}
        return all(amino_acid in valid_amino_acids for amino_acid in sequence)

    def molecular_weight(self) -> float:
        """Вычисляет молекулярный вес аминокислотной последовательности."""
        amino_acid_weights = {
            'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
            'Q': 146.14, 'E': 147.13, 'G': 75.07, 'H': 155.15, 'I': 131.17,
            'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
            'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
        }
        return sum(amino_acid_weights.get(aa, 0) for aa in self._sequence)


def filter_fastq(
    input_path: str, 
    output_path: str, 
    gc_bounds: int | float | tuple = (0, 100), 
    length_bounds: int | tuple = (0, 2**32), 
    quality_threshold: float | int = 0
) -> dict:
    """Фильтрация FastQ-файлов по длине, качеству и GC-составу с использованием Biopython."""
    
    def filter_record(record: SeqRecord) -> bool:
        """Возвращает True, если запись удовлетворяет заданным критериям."""
        seq_len: int = len(record.seq)
        mean_quality: float = sum(record.letter_annotations["phred_quality"]) / seq_len
        gc_content: float = gc_fraction(record.seq)
        
        return (length_bounds[0] <= seq_len <= length_bounds[1] and
                mean_quality >= quality_threshold and
                gc_bounds[0] <= gc_content <= gc_bounds[1])
        
    # Читаем FastQ-файл и фильтруем записи
    with open(input_path, "r") as input_handle, open(output_path, "w") as output_handle:
        records = SeqIO.parse(input_handle, "fastq")
        filtered_records = (record for record in records if filter_record(record))

        # Записываем отфильтрованные записи в новый FastQ-файл
        SeqIO.write(filtered_records, output_handle, "fastq")
