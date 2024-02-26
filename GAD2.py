import os
import re
from typing import TextIO, Optional, Union
from abc import ABC, abstractmethod
from Bio import SeqIO, SeqUtils

class BiologicalSequence(ABC):
    @abstractmethod
    def __len__():
        pass

    def __getitem__():
        pass

    def __str__():
        pass

    def check_sequence():
        pass

class NucleicAcidSequence(BiologicalSequence):
    """Proccising nucleic acid sequences.
    This class is parental for DNASequence and RNASequence classes.
    """
    def __init__(self, seq) -> None:
        self.seq = seq

    def check_sequence(self):
        """
        Checkes the correctnessis of input string (DNA or RNA).
        """
        for nucleotide in self.seq:
            if not nucleotide in type(self).complement_rule:
                raise ValueError (f'Incorrect nucleotide in {self.seq}')
        return True

    def complement(self):
        """
        Create complement sequences of RNA or DNA acoording to the complementary base pairing rule. 
        """
        if self.check_sequence():
            complement_sequence = ''
            for nucleotide in self.seq:
                if nucleotide in type(self).complement_rule:
                    complement_sequence += type(self).complement_rule[nucleotide]
            return type(self)(complement_sequence)

    def gc_content(self):
        return (sum(1 for _ in re.finditer(r'[GCgc]', self.seq)))/self.__len__()
    
    def __len__(self):
        return len(self.seq)

    def __getitem__(self, key):
        return self.seq[key]

    def __str__(self):
        return self.seq

class DNASequence(NucleicAcidSequence):
    """
    Procissing DNA sequences.
    Argumemt: seq (str) - DNA sequence, letters case do not matter.
    -Example nucl_seq = DNASequence('ATGCGC' OR 'ATgCgc' OR 'atgcgc).

    Valid operations:
    - transcribe - return transcibed sequence of DNA coding strand;
    - complement - return complementary sequence according to complementary rule
    - check sequence - is input sequence correct
    - gc_content - return percent of GC in sequence
    """
    complement_rule = {'a': 't', 'A': 'T', 't': 'a', 'T': 'A',
                       'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}
    
    def __init__(self, seq: str) -> None:
        super().__init__(seq = seq)

    def transcribe(self):
        """
        Return transcribed sequence of DNA acoording to the complementary base pairing rule.  
        Function can procced only DNA seqences.
        """
        transcribed_sequence = ''
        transcribed_sequence = self.seq.replace('t', 'u').replace('T', 'U')
        return type(self)(transcribed_sequence)

class RNASequence(NucleicAcidSequence):
    """
    Procissing DNA sequences.
    Argumemt: seq (str) - RNA sequence, letters case do not matter.
    -Example nucl_seq = RNASequence('AUGCGC' OR 'AUgCgc' OR 'augcgc).

    Valid operations:
    - complement - return complementary sequence according to complementary rule
    - check sequence - is input sequence correct
    - gc_content - return percent of GC in sequence
    """
    complement_rule = {'a': 'u', 'A': 'U', 'u': 'a', 'U': 'A',
                       'g': 'c', 'G': 'C', 'c': 'g', 'C': 'G'}
    
    def __init__(self, seq) -> None:
        super().__init__(seq = seq)


def check_gc(fastq_read: str, gc_params: tuple) -> bool:
    """
    Filters sequences in FASTQ file by GC percentage. 
    
    Input:
    - fastq_read (str): FASTQ sequence read.
    - gc_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """
    gc_result = SeqUtils.GC123(fastq_read)[0]
    if gc_params[0] < gc_result < gc_params[1]:
        return True


def check_length(fastq_read: str, len_bound_params: tuple) -> bool:
    """
    Filters sequences in FASTQ file by sequence length.
    
    Input:
    - fastq_read (str): FASTQ sequence read.
    - len_bound_params (tuple): range for filtration, accepts upper or upper/lower border. Both borders are included.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """
    len_of_seq = len(fastq_read)
    if len_bound_params[0] < len_of_seq < len_bound_params[1]:
        return True


def check_quality(fastq_quality: str, quality_params: int) -> bool:
    """
    Filters sequences in FASTQ file by mean quality score.
    
    Input:
    - fastq_quality (str): FASTQ read quality
    - quality_params (int): threshold value of reads quality in phred33 scale.
    
    Output:
    Returns boolean value is this read satisfied the criteria.
    """ 
    mean_quality = sum(fastq_quality.letter_annotations["phred_quality"])/len(fastq_quality.seq)
    if mean_quality >= quality_params:
        return True


def int_to_tuple(input_parameters) -> tuple:
    """
    Converts input parameters to tuple format.
    If input is already a tuple, it will be return without changes.
    If input parameter is int (only upper threshold is entered), function will return a tuple like (0, 'input').
    
    Input:
    - input_parameters (tuple or int).
    
    Output:
    Lower and upper threshold for functions in tuple format.
    """
    if isinstance(input_parameters, tuple):
        return input_parameters
    return (0, input_parameters)


def run_fastq_filter(input_path: Optional[str] = None, output_filename: Optional[str] = None, gc_bounds: Union[int, tuple] = (0, 100), length_bounds: Union[int, tuple] = (0, 2**32), quality_threshold: int = 0) -> TextIO:
    """
    FASTQ Filtraror with BIOPYTHON utilities
    Performs filter of input FASTQ file according to input parameters. 
    Input will be filtered by: 
        - GC content (gc_bounds);
        - reads length (length_bounds);
        - reads quality score (quality_threshold).
        
    Input:
    - input_path (str): path to .fastq file; include 4 strings: 1 - read ID, 2 - sequence, 3 - comment, 4 - quality. Default - None.
    - output_filename (str): name of output file, by default, it will be saved in the directory 'fastq_filtrator_resuls'. Default name will be name of input file.
    - gc_bounds (tuple or int): GC content filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 100).
    - length_bounds (tuple or int): read length filter parameters, it accepts lower and upper (tuple), or only upper threshold value (int). Default value (0, 2**32).
    - quality_threshold (int): upper quality threshold in phred33 scale. Reads with average quality below the threshold are discarded. Default value - 0. 
    
    Output:
    Returns FASTQ only with filtered reads which satisfied all input/default conditions.
    """

    "Specify input path"
    if not input_path.endswith('.fastq'):
        raise ValueError('Incorrect input file extension, should be .fastq')   

    "Specify output path"
    output_path = 'fastq_filtrator_resuls'
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if output_filename is None:
        output_filename = os.path.basename(input_path)
    "Passed the parameters"
    gc_params = int_to_tuple(gc_bounds)
    len_bound_params = int_to_tuple(length_bounds)    
    "Filter and record results"
    filtererd_fastq = open(os.path.join(output_path, output_filename), mode='w')
    for seq_record in SeqIO.parse(input_path, "fastq"):
        if check_gc(seq_record.seq, gc_params) and check_length(seq_record.seq, len_bound_params) and check_quality(seq_record, quality_threshold):
                SeqIO.write(seq_record, filtererd_fastq, "fastq")  
    filtererd_fastq.close()   
    return filtererd_fastq
