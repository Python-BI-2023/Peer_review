import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from typing import Union, Tuple, Dict
from abc import ABC, abstractmethod



def analyse_gc(records: str, min_gc: Union[int, float], max_gc: Union[int, float]) -> str:
    """
    Return filtered FASTQ-sequences by GC-content
    
    :param records: FASTQ-file
    :type records: str
    :param min_gc: Left boundary for GC-content filtration
    :type min_gc: Union[int, float]
    :tupe max_gc: Right boundary for GC-content filtration
    :type max_gc: Union[int, float]
    :rtype: str
    :return: filtered FASTQ-sequences
    """
    return {record.id: record for record in records if (min_gc/100) <= gc_fraction(record.seq) <= (max_gc/100)}


def filter_by_length(records: str, min_length: Union[int, float], max_length: Union[int, float]) -> str:
    """
    Return filtered FASTQ-sequences by GC-content
    
    :param records: FASTQ-file
    :type records: str
    :param min_gc: Left boundary for GC-content filtration
    :type min_gc: Union[int, float]
    :tupe max_gc: Right boundary for GC-content filtration
    :type max_gc: Union[int, float]
    :rtype: str
    :return: filtered FASTQ-sequences
    """

    return {record.id: record for record in records if min_length <= len(record.seq) <= max_length}


def filter_by_quality(records: str, quality_threshold: Union[int, float]) -> str:
    """
    Return filtered FASTQ-sequences by quality
    :param records: Sequnces filtered by GC-content and length
    :type records: str
    :param quality_threshold: boundary for quality filtration
    :type quality_threshold: Union[int, float]
    :rtype: str
    :return: filtered FASTQ-sequences   
    """
    
    return {record.id: record for record in records if min(record.letter_annotations["phred_quality"]) >= quality_threshold}



def write_filtered_sequences_to_fastq(filtered_sequences: Dict[str, str], unfiltered_sequences: Dict[str, str],
                                      output_file: str, folder_path: str = 'fastq_filtrator_results') -> str:
    '''
    The function writes filtered FASTQ reads into new file and save it into folder "fastq_filtrator_resuls"
    :param filtered_sequences: Dict of filtered sequences by GC-content, quality and length
    :filtered_sequences type: Dict[str]
    :param unfiltered_sequences: Dict of unfiltered sequences by GC-content, quality and length
    :filtered_sequences type: Dict[str]
    :output_file: name of output file
    :output_file type: str
    :folder_path: name of result folder, default = 'fastq_filtrator_resuls'
    :folder_path type: str
    :rtype: None
    :return: None
    '''

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    file_path_filtered = os.path.join(folder_path, output_file)
    file_path_not_filtered = os.path.join(folder_path, 'unfiltered_sequences.fasta')
    
    with open(file_path_filtered, "w") as output_handle:
        for key, value in filtered_sequences.items():
            records_filtered = []
            sequence = Seq(value)
            record = SeqRecord(sequence, id=key, description="")
            records_filtered.append(record) 
            SeqIO.write(records_filtered, output_handle, "fasta")

    with open(file_path_not_filtered, "w") as output_handle:
        for key, value in unfiltered_sequences.items():
            records_unfiltered = []
            sequence = Seq(value)
            record = SeqRecord(sequence, id=key, description="")
            records_unfiltered.append(record) 
            SeqIO.write(records_unfiltered, output_handle, "fasta")


def filter_fastq(input_path: str, 
                  gc_bounds: Union[int, float, Tuple [int], Tuple [float]] = (0, 100), 
                  length_bounds: Union[int, Tuple [int]] = (0, 2**32),
                  quality_threshold: float = 0.0, filtered_file_name: Union[None, str] = None) -> Dict[str,str]:
    """
    This function help analyze a set of reads obtained from next-generation sequencing. 
    
    The function allow to filter the desired reads according to three parameters:
    GC-content, length and reading quality.
    
    :param seqs: 
    Path to the file with FASTQ-sequences in the format. 
    :type seqs: str
    
    :param gc_bounds: 
    Boundary parameters for filtering sequences by GC-content. Save only reads with a GC-content between boundaries 
    or lower than one boundary. Lower boundary cannot be less than 0 and upper boundary cannot be greater than 100. 
    gc_bounds default value is (0,100).
    :type param gc_bounds: Union[int, float, Tuple [int], Tuple [float]]
    
    :param length_bounds: 
    Boundary parameters for filtering sequences by length. Works the same as gc_bounds. Lower boundary cannot be less 
    than 0 and upper boundary cannot be greater than 2^32. length_bounds default value is (0,2^32)
    :type param length_bounds: Union[int, Tuple [int]
    
    :param quality_threshold: 
    Threshold for quality of each nucleotide in read. Quality incodes by ASCII codes. The threshold cannot be more 
    than 40. quality_threshold default value is 0 
    :type param quality_threshold: float
    
    :return: 
    New dictionaries with fastq sequence.The first one consisting of filtered fastq sequences and the other one with 
    sequences that did not pass filters.
    :rtype: Dict[str]
    
    :raises ValueError: if sequence not RNA or DNA, also if the argument values are outside the allowed ones
    """

    if type(gc_bounds) == float or type(gc_bounds) == int:
        gc_bounds = (0, gc_bounds)
        if gc_bounds[0] < 0 or gc_bounds[1] > 100:
            raise ValueError(f'Wrong boundaries!')
    min_gc, max_gc = gc_bounds
    if type(length_bounds) == int:
        length_bounds = (0,length_bounds)
        if length_bounds[0] < 0 or length_bounds[1] > 2**32:
            raise ValueError(f'Wrong boundaries!')    
    min_length, max_length = length_bounds    
    if quality_threshold > 40:
        raise ValueError(f'Wrong quality threshold!')

    records = list(SeqIO.parse(input_path, "fastq"))

    filtered_by_gc = analyse_gc(records, min_gc, max_gc)

    filtered_by_length = filter_by_length(filtered_by_gc.values(), min_length, max_length)

    filtered_by_quality = filter_by_quality(filtered_by_length.values(), quality_threshold)

    filtered_seq = {record_id: str(record.seq) for record_id, record in filtered_by_quality.items()}
    unfiltered_seq = {record.id: str(record.seq) for record in records if record.id not in filtered_seq}

    if filtered_file_name == None:
        new_file_name = "filtered_sequences.fasta"
    else:
        new_file_name = filtered_file_name

    write_filtered_sequences_to_fastq(filtered_seq, unfiltered_seq, new_file_name)

    return ('Sequences are filtered!')


class BiologicalSequence(ABC):
    '''
    Abstract class for different biological sequences
    '''

    @abstractmethod
    def __len__(self):
        '''
        Method for working with the Python len function. !Needs to be overridden in child class!
        '''
        
        pass

    @abstractmethod
    def __getitem__(self):
        '''
        Method for get elements by index and slice the sequence. !Needs to be overridden in child class!
        '''
        
        pass

    @abstractmethod
    def __str__(self):
        '''
        Method for convertion sequence to a string. !Needs to be overridden in child class!
        '''
        
        pass

    @abstractmethod
    def is_alphabet_correct(self):
        '''
        Method for checking that a sequence is written correctly. 
        '''
        pass


class NucleicAcidSequnce(BiologicalSequence):
    '''
    Class for DNA or RNA molecules
    '''

    def __init__(self, seq) -> None:
        self.seq = seq
        self.dna_alphabet = set('AaTtGgCc')
        self.rna_alphabet = set('AaUuGgCc')
        self.complement_dict = {'A': 'T', 'C': 'G', 
                    'G': 'C', 'T': 'A', 'U': 'A', 'a': 't',
                    'c': 'g', 'g': 'c', 't': 'a', 'u': 'a'}        


    def __len__(self) -> int:
        return len(self.seq)
    

    def __getitem__(self, item) -> int:
        return self.seq[item]
    

    def __str__(self) -> str:
        return self.seq
    


    def is_alphabet_correct(self) -> bool:

        '''
        Method for checking of standard nucleotide content in sequence

        :param self: DNA or RNA sequence
        '''
        
        if (set(self.seq).issubset(self.dna_alphabet) and isinstance(self, DNASequence)) or (set(self.seq).issubset(self.rna_alphabet) and isinstance(self, RNASequence)):
            return True
        raise TypeError(f'{self.seq} is not correct nucleic acid')          


    def complement(self):
        """
        Function return complement sequence.
        
        :param self: DNA or RNA sequence
        :rtype: str
        :return: complement sequence   
        """
        if self.is_alphabet_correct():
            complement_seq = str()
            length = len(self.seq)
            for i in range (length):
                if self.seq[i] in self.complement_dict:
                    complement_seq += (self.complement_dict[self[i]])
            if isinstance(self, DNASequence):
                return DNASequence(complement_seq)
            if isinstance(self, RNASequence):
                return RNASequence(complement_seq)
        
    
    def gc_calculate(self) -> float:
        """
        Function return sequence GC-content in percent.
        
        :param seq: DNA or RNA sequence
        :type seq: str
        :rtype: float
        :return: GC-contentn percent 
        """
        length = len(self.seq)
        gc_content = 0.0
        seq_up = self.seq.upper()
        c = seq_up.count("C")
        g = seq_up.count("G")
        gc_content = round(((c+g)/length*100),2)
        return gc_content


class DNASequence(NucleicAcidSequnce):
    '''
    Class for DNA sequence
    '''

    def __init__(self, seq) -> None:
        super().__init__(seq)

    def transcribe(self):
        '''
        Method return transcribed sequence.
        '''

        if super().is_alphabet_correct():
            return RNASequence(self.seq.replace('T', 'U').replace('t', 'u'))
    

class RNASequence(NucleicAcidSequnce):
    '''
    Class for RNA sequence
    '''
    
    def __init__(self, seq) -> None:
        super().__init__(seq)


class AminoAcidSequence (BiologicalSequence):
    '''
    Class for protein sequence
    '''
    
    def __init__(self, seq):
        self.seq = seq
        self.protein_alphabet = set('ACDEFGHIKLMNPQRSTVWY')


    def __len__(self):
        return len(self.seq)


    def __getitem__(self, item):
        return self.seq[item]


    def __str__(self):
        return self.seq


    def is_alphabet_correct(self):
        if set(self.seq).issubset(self.protein_alphabet):
            return True
        raise TypeError(f'{self.seq} is not a protein')    

    def calculate_protein_mass(self):
        '''
        Method return mass of residues in seq in Da.
        '''
        
        if self.is_alphabet_correct:
            weights = {'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.15,
                   'E': 147.13, 'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17,
                   'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
                   'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15}
            return sum(weights.get(aa, 0) for aa in self.seq)
