# Importing modules
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Union


def fastq_filter(input_path: str = None, output_filename: str = None, *,
                 gc_bound: Union[tuple, int, float] = (0, 100),
                 length_bound: Union[tuple, int, float] = (0, 2**32),
                 quality_threshold: Union[int, float] = 0) -> None:
    """
    This function work with FASTQ files and filters them by
    GC content, length and Q-score.

    Arguments (positional):
    - input_path (str): full path to the file that you want to work with
    - output_filename (str): enter just a name of the file, don't add extention

    Arguments (keyword):
    - gc_bound (tuple, int, float): tuple of required range of GC percentage (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - length_bound (tuple, int, float): tuple of required range of sequences length (inclusive),
    num or float if only higher border of the range is needed (exclusive).
    - quality_threshold (int): int of lowest level of Q-score (inclusive).

    Output:
    - list of BioSeq records. Write file to .fastq
    """
    # Chech if PATH for input file is given
    if input_path is None:
        raise ValueError("You didn't enter any PATH to file")
    # Chech if folder exist and create outout_path if not given
    input_folder = input_path.rsplit('/', 1)[0]
    input_name = input_path.rsplit('/', 1)[1]
    is_exist = os.path.exists(f'{input_folder}/fastq_filtrator_resuls/')
    if not is_exist:
        os.makedirs(f'{input_folder}/fastq_filtrator_resuls/')
    if output_filename is None:
        output_path = f'{input_folder}/fastq_filtrator_resuls/{input_name}'
    else:
        output_path = f'{input_folder}/fastq_filtrator_resuls/{output_filename}.fastq'
    # Create dict from FASTQ
    seqs = list(SeqIO.parse(input_path, "fastq"))
    # Check that this dict is not empty
    if len(seqs) <= 0:
        raise ValueError('There are no fastq sequences')
    # Check if all given argumets have relevant type
    gc_bound_type = isinstance(gc_bound, (tuple, int, float))
    length_bound_type = isinstance(length_bound, (tuple, int, float))
    quality_thr_type = isinstance(quality_threshold, (int, float))
    if not (gc_bound_type and length_bound_type and quality_thr_type):
        raise ValueError('Your arguments are not suitable!')
    # Create filtered list
    filtered_fastq = []

    for line in seqs:
        if isinstance(gc_bound, (int, float)):
            gc_check = gc_fraction(line)*100 >= gc_bound
        else:
            gc_check = gc_fraction(line)*100 < gc_bound[0] or gc_fraction(line)*100 > gc_bound[1]
        if isinstance(length_bound, (int, float)):
            len_check = len(line) >= length_bound
        else:
            len_check = len(line) < length_bound[0] or len(line) > length_bound[1]
        quality_check = sum(line.letter_annotations['phred_quality'])/len(line.letter_annotations['phred_quality']) < quality_threshold
        if not (gc_check or len_check or quality_check):
            filtered_fastq.append(line)

    # Write  filtered data into new .fastq file
    SeqIO.write(filtered_fastq, output_path, 'fastq')

# Custom error
class WrongSequence(ValueError):
    pass


class BiologicalSequence(str):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)
    
    def slice(self, start_index, stop_index):
        return self.sequence[start_index:stop_index]
    
    def alphabet_checking(self):
        if not set(self.sequence) <= set(type(self).dictionary.keys()):
            raise WrongSequence('Wrong sequence')
        return True


class NucleicAcidSequence(BiologicalSequence):
    
    def __init__(self, sequence):
        super().__init__(sequence)
        if not self.alphabet_checking():
            del self.sequence
            raise WrongSequence('You have entered a wrong sequence')
        self.gc_cont = None

    def complement(self):
        return type(self)(''.join([type(self).dictionary[i] for i in self.sequence]))

    def gc_content(self):
        gc = 0
        for nucleotide in self.sequence:
            if nucleotide in set('CGcg'):
                gc += 1
        self.gc_cont = round(gc / len(self.sequence), 3)
        return self.gc_cont


class DNAsequence(NucleicAcidSequence):
    dictionary = {
        'A': 'T',
        'G': 'C',
        'T': 'A',
        'C': 'G',
        'a': 't',
        'g': 'c',
        't': 'a',
        'c': 'g',
}
    trans_dict = {
        'A': 'U',
        'G': 'C',
        'T': 'A',
        'C': 'G',
        'a': 'u',
        'g': 'c',
        't': 'a',
        'c': 'g',
}
    def transcribe(self):
        return RNAsequence(''.join([self.trans_dict[i] for i in self.sequence]))


class RNAsequence(NucleicAcidSequence):
    dictionary = {
        'A': 'U',
        'G': 'C',
        'U': 'A',
        'C': 'G',
        'a': 'u',
        'g': 'c',
        'u': 'a',
        'c': 'g',
}


class AminoAcidSequence(BiologicalSequence):
    dictionary = {
        "A": 71.03711,
        "C": 103.00919,
        "D": 115.02694,
        "E": 129.04259,
        "F": 147.06841,
        "G": 57.02146,
        "H": 137.05891,
        "I": 113.08406,
        "K": 128.09496,
        "L": 113.08406,
        "M": 131.04049,
        "N": 114.04293,
        "P": 97.05276,
        "Q": 128.05858,
        "R": 156.10111,
        "S": 87.03203,
        "T": 101.04768,
        "V": 99.06841,
        "W": 186.07931,
        "Y": 163.06333,
    }
    def __init__(self, sequence):
        super().__init__(sequence)
        if not self.alphabet_checking():
            del self.sequence
            raise WrongSequence('You have entered a wrong sequence')

    def protein_mass(self):
        mass = sum(self.dictionary.get(aa) for aa in self.sequence)
        return mass
