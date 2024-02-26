import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Tuple


def filter_fastq(input_path: str, 
                 output_filename: str = None, 
                 gc_bounds: Tuple[int, int] = (0, 100), 
                 length_bounds: Tuple[int, int] = (0, 2**32), 
                 quality_threshold: int = 0) -> None:
    '''
    Filter FASTQ-sequences based on entered requirements.
    
    Arguments:
        - input_path (str): path to the file with FASTQ-sequences
        - output_filename (str): name of the output file with 
        filtered FASTQ-sequences
        - gc_bounds (tuple or int, default = (0, 100)): GC-content
        interval (percentage) for filtering. Tuple if contains 
        lower and upper bounds, int if only contains an upper bound.
        - length_bounds (tuple or int, default = (0, 2**32)): length 
        interval for filtering. Tuple if contains lower and upper 
        bounds, int if only contains an upper bound.
        - quality_threshold (int, default = 0): threshold value of average 
        read quality for filtering.

    Note: the output file is saved to the /fastq_filtrator_results 
    directory. The default output file name is the name of the input file.
    '''
    with open(input_path, "r") as handle:
        records = [record for record in SeqIO.parse(handle, "fastq")]

    filtered_records = []
    for record in records:
        gc_content = gc_fraction(record.seq) * 100
        seq_length = len(record.seq)
        avg_quality = sum(record.letter_annotations["phred_quality"]) / seq_length

        if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
            length_bounds[0] <= seq_length <= length_bounds[1] and
            avg_quality >= quality_threshold):
            filtered_records.append(record)

    if not os.path.exists("fastq_filtrator_results"):
        os.makedirs("fastq_filtrator_results")

    if output_filename is None:
        output_filename = input_path.split("/")[-1].split(".")[0] + "_filtered.fastq"

    output_path = "fastq_filtrator_results/" + output_filename
    SeqIO.write(filtered_records, output_path, "fastq")

    print(f"Filtered sequences saved to {output_path}")


class BiologicalSequence(ABC):
    def __init__(self, seq: str = None):
        self.seq = seq

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __len__(self):
        return len(self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __repr__(self):
        return self.seq

    def check_alphabet(self):
        return set(self.seq.upper()).issubset(self.ALPHABET)

    def complement(self):
        comp_seq = self.seq.translate(str.maketrans(self.COMPLEMENT_MAP))
        return type(self)(comp_seq)
    
    def reverse_complement(self):
        reverse_seq = self.complement().seq[::-1]
        return type(self)(reverse_seq)
    
    def gc_content(self):
        gc_count = self.seq.count('G') + self.seq.count('C')
        return gc_count / len(self.seq) * 100
    

class DNASequence(NucleicAcidSequence):
    ALPHABET = {'A', 'T', 'G', 'C'}
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                      'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    def transcribe(self):
        return RNASequence(self.seq.replace('T', 'U'))


class RNASequence(NucleicAcidSequence):
    ALPHABET = {'A', 'G', 'C', 'U'}
    COMPLEMENT_MAP = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                      'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


class AminoAcidSequence(BiologicalSequence):
    ALPHABET = {"A", "C", "D", "E", "F", "G", "H", "I","K", "L", 
                "M", "N","P", "Q", "R", "S", "T", "V", "W", "Y"}
    HYDROPHOBIC_AMINOACIDS = {"A", "V", "L", "I", "M", "F", "Y", "W"}
    
    def __len__(self):
        return len(self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __repr__(self):
        return self.seq

    def check_alphabet(self):
        return set(self.seq.upper()).issubset(self.ALPHABET)
    
    def compute_hydrophobicity(self):
        hydrophobic_count = sum(1 for aa in self.seq if aa.upper() in self.HYDROPHOBIC_AMINOACIDS)
        return (hydrophobic_count / len(self.seq)) * 100
