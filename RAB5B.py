import numpy as np
from abc import ABC, abstractmethod
from Bio import SeqIO, SeqRecord, SeqUtils


def length_check(sequence, length_bounds):
    """
    Checking if length of a sequence is inside bounds
    """
    if type(length_bounds) is tuple:
        min, max = length_bounds
    else:
        min, max = 0, length_bounds
    return min <= len(sequence) <= max


def gc_check(sequence, gc_bounds):
    """
    Checking if amount of GC in a sequence is inside bounds
    """
    if type(gc_bounds) is tuple:
        min, max = gc_bounds
    else:
        min, max = 0, gc_bounds
    G_amount = sequence.count('G')
    C_amount = sequence.count('C')
    GC_amount = 100 * (G_amount + C_amount) / len(sequence)
    if GC_amount >= min and GC_amount <= max:
        return True


def score_check(quality, quality_threshold):
    """
    Checking if mean q-score of a sequence passes the threshold
    """
    q_scores = []
    for char in quality:
        q_scores.append(ord(char) - 33)
    mean_q_score = sum(q_scores) / len(q_scores)
    if mean_q_score >= quality_threshold:
        return True


def write_output(name, sequence,
                 commentary, quality,
                 gc_bounds, length_bounds,
                 quality_threshold,
                 passed, failed,
                 save_filtered):
    """
    Filtering sequences based on their parameters
    and writing filtered sequences into different output files
    """
    length_result = length_check(sequence, length_bounds)
    GC_result = GC_check(sequence, gc_bounds)
    score_result = score_check(quality, quality_threshold)
    if length_result and GC_result and score_result:
        print(name, sequence, commentary, quality, file=passed, sep='\n')
    elif save_filtered == True:
        print(name, sequence, commentary, quality, file=failed, sep='\n')


def main(input_fastq, output_file_prefix,
         gc_bounds=(0, 100), length_bounds=(0, 2 ** 32),
         quality_threshold=0, save_filtered=False):
    """
    Filters fastq file based on several parameters

    :param input_fastq: path to file as 'file'
    :param output_file_prefix: prefix of filtered output files
    :param gc_bounds: allowed GC-percent given as a tuple or as one upper bound,
    :param length_bounds: allowed length of a sequence, given a as tuple or as
    one upper bound
    :param quality_threshold: reads with quality below this threshold will be
    filtered
    :param save_filtered: save reads that did not pass filters

    """
    with open(input_fastq) as fastq:
        with open(f'{output_file_prefix}_passed.fastq', 'w') as passed:
            with open(f'{output_file_prefix}_failed.fastq', 'w') as failed:
                while True:
                    name = fastq.readline().strip()
                    if name == '':
                        break
                    else:
                        sequence = fastq.readline().strip()
                        commentary = fastq.readline().strip()
                        quality = fastq.readline().strip()
                    write_output(name, sequence,
                                 commentary, quality,
                                 gc_bounds, length_bounds,
                                 quality_threshold,
                                 passed, failed,
                                 save_filtered)


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def alphabet_check(self):
        pass


class Sequence(BiologicalSequence):
    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        return str(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def alphabet_check(self):
        alphabet = set(self.sequence)
        if not alphabet.issubset(self.allowed):
            raise ValueError('Wrong alphabet')


class NucleicAcidSequence(Sequence):
    def __init__(self, sequence):
        self.sequence = sequence
        self.allowed = {'A', 'T', 'G', 'C', 'U'}

    def complement(self):
        self.complemented = []
        for letter in self.sequence:
            self.complemented.append(self.complement_dict[letter])
        return type(self)(''.join(self.complemented))

    def gc_content(self):
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        content = (g_count + c_count) / self.sequence.__len__()
        return content


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence)
        self.allowed = {'A', 'T', 'G', 'C'}
        self.complement_dict = {'A': 'T', 'T': 'A',
                                'G': 'C', 'C': 'G'}

    def transcribe(self):
        return self.sequence.replace('T', 'U')


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super.__init__(sequence)
        self.allowed = {'A', 'C', 'G', 'U'}
        self.complement_dict = {'A': 'U', 'U': 'A',
                                'G': 'C', 'C': 'G'}


class AminoAcidSequence(Sequence):
    def __init__(self, sequence):
        self.sequence = sequence
        self.allowed = {'A', 'R', 'N', 'D',
                        'C', 'Q', 'E', 'G',
                        'H', 'I', 'L', 'K',
                        'M', 'F', 'P', 'S',
                        'T', 'W', 'Y', 'V'}

        self.composition = {'A': 0, 'R': 0, 'N': 0, 'D': 0,
                            'C': 0, 'Q': 0, 'E': 0, 'G': 0,
                            'H': 0, 'I': 0, 'L': 0, 'K': 0,
                            'M': 0, 'F': 0, 'P': 0, 'S': 0,
                            'T': 0, 'W': 0, 'Y': 0, 'V': 0}

    def count_composition(self):
        for am_acid in self.sequence:
            self.composition[am_acid] += 1
        for am_acid in self.composition:
            self.composition[am_acid] = round(self.composition[am_acid] * 100 / len(self.sequence), 2)

        return self.composition


def length_check(sequence, length_bounds):
    """
    Checking if length of a sequence is inside bounds
    """
    if type(length_bounds) is tuple:
        min, max = length_bounds
    else:
        min, max = 0, length_bounds
    if len(sequence) >= min and len(sequence) <= max:
        return True

def GC_check(sequence, gc_bounds):
    """
    Checking if amount of GC in a sequence is inside bounds
    """
    if type(gc_bounds) is tuple:
        min, max = gc_bounds
    else:
        min, max = 0, gc_bounds
    G_amount = sequence.count('G')
    C_amount = sequence.count('C')
    GC_amount = 100*(G_amount + C_amount)/len(sequence)
    if GC_amount >= min and GC_amount <= max:
        return True

def score_check(quality, quality_threshold):
    """
    Checking if mean q-score of a sequence passes the threshold
    """
    q_scores = []
    for char in quality:
        q_scores.append(ord(char)-33)
    mean_q_score = sum(q_scores)/len(q_scores)
    if mean_q_score >= quality_threshold:
        return True


def write_output(name, sequence,
                 commentary, quality,
                 gc_bounds, length_bounds,
                 quality_threshold,
                 passed, failed,
                 save_filtered):
    """
    Filtering sequences based on their parameters 
    and writing filtered sequences into different output files
    """
    length_result = length_check(sequence, length_bounds)
    GC_result = GC_check(sequence, gc_bounds)
    score_result = score_check(quality, quality_threshold)
    if length_result and GC_result and score_result:
        print(name, sequence, commentary, quality, file = passed, sep = '\n')
    elif save_filtered == True:
        print(name, sequence, commentary, quality, file = failed, sep = '\n')


def main(input_fastq, output_file_prefix,
         gc_bounds = (0, 100), length_bounds = (0, 2**32),
         quality_threshold = 0, save_filtered = False):
    """
    Filters fastq file based on several parameters
    
    :param input_fastq: path to file as 'file'
    :param output_file_prefix: prefix of filtered output files
    :param gc_bounds: allowed GC-percent given as a tuple or as one upper bound,
    :param length_bounds: allowed length of a sequence, given a as tuple or as
    one upper bound
    :param quality_threshold: reads with quality below this threshold will be 
    filtered
    :param save_filtered: save reads that did not pass filters

    """
    with open(input_fastq) as fastq:
        with open(f'{output_file_prefix}_passed.fastq', 'w') as passed:
            with open(f'{output_file_prefix}_failed.fastq','w') as failed:
                while True:
                    name = fastq.readline().strip()
                    if name == '':
                        break
                    else:
                        sequence = fastq.readline().strip()
                        commentary = fastq.readline().strip() 
                        quality = fastq.readline().strip()
                    write_output(name, sequence,
                        commentary, quality,
                        gc_bounds, length_bounds,
                        quality_threshold,
                        passed, failed,
                        save_filtered)
