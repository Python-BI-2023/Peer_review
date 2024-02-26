import os
from abc import ABC, abstractmethod
from typing import Optional

import numpy as np
from Bio import SeqIO, SeqRecord, SeqUtils


def is_in_gc_bounds(seq_record: SeqRecord, gc_bounds: tuple) -> bool:
    """
    Check if the sequence is in the range of GC-content bounds

    Arguments:
    - seq_record(SeqRecord): the sequence object to check
    - gc_bounds(tuple): contain min and max GC-content bounds

    Return:
    - bool: the result of the check
    """
    gc_min, gc_max = gc_bounds[0], gc_bounds[1]
    gc_content = SeqUtils.GC123(seq_record.seq)[0]
    return gc_min <= gc_content <= gc_max


def is_in_length_bounds(seq_record: SeqRecord, length_bounds: tuple) -> bool:
    """
    Check if the sequence is in the range of length bounds

    Arguments:
    - seq(SeqRecord): the sequence object to check
    - length_bounds(tuple): contain min and max length bounds

    Return:
    - bool: the result of the check
    """
    length_min, length_max = length_bounds[0], length_bounds[1]
    return length_min <= len(seq_record) <= length_max


def is_above_quality_threshold(seq_record: SeqRecord,
                               quality_threshold: float) -> bool:
    """
    Check if the mean of the sequence quality values in FASTQ exceeds
    the quality threshold of the interest

    To convert quality values, the function uses phred+33

    Arguments:
    - seq_record(Bio.SeqRecord): the sequence object to check
    - quality_threshold(int or float): the quality threshold of the interest

    Return:
    - bool: the result of the check
    """
    quality = np.mean(seq_record.letter_annotations["phred_quality"])
    return quality >= quality_threshold


def filter_fastq(input_path: str,
                 output_filename: Optional[str] = None,
                 gc_bounds=(0, 100),
                 length_bounds=(0, 2**32),
                 quality_threshold=0) -> None:
    """
    Main function to select reads in FASTQ format according
    to three main requirements:

    -correspondence to the range of GC-content bounds
    The range of GC-content bounds is determined with gc_bounds argument

    -correspondence to the range of length bounds
    The range of length bounds is determined with length_bounds argument

    -exceeds the quality threshold of the interest
    The quality threshold is determined with quality_threshold argument

    Function takes the path to the file in input_path argument.
    Use files with fastq extension only.

    Function output the result of checking in a file that is named according
    to output_filename(Optional).

    The output file also has fastq extension.

    The output file is saved in the 'fastq_filtrator_results' directory.

    If 'fastq_filtrator_results' directory doesn't exist the program creates it
    in a current directory.

    Without full names use arguments in a certain order
    Example: filter_fastq(input_path, output_filename, (0,100), (0,200), 0)
           # filter_fastq(input_path, output_filename, gc_bounds=(0,100),
             length_bounds=(0,200), quality_threshold=0)

    In case of changing only one argument, provide its full name!
    Example: filter_fastq(input_path, output_filename,
                          length_bounds=(50, 100))

    Arguments:

    - input_path(str): the path to the file

    - output_filename(str): the name for output file with obtained result
    By default output_filename=None
    Without output_filename argument the output file is named as input file
    Name without fastq extention is acceptible.
    Example: output_filename='result'  # 'result.fastq'
             output_filename='result.fastq'

    - gc_bounds(tuple or int or float): contain min and max GC-content bounds
    By default gc_bounds=(0,100)
    If input contains one number the function accepts it as a maximum bound
    Examples: gc_bounds=(20,40)
              gc_bounds=40  # (0,40)

    - length_bounds(tuple or int or float): contain mini and max length bounds
    By default length_bounds=(0,4294967296)
    If input contains one number the function accepts it as a maximum bound
    Examples: length_bounds=(10,90)
              length_bounds=90  # (0,90)

    - quality_threshold(int or float): the quality threshold of the interest
    By default quality_threshold=0
    Examples: quality_threshold=10

    There are three functions that are used in the main function:

    - is_in_gc_bounds(seq_record, gc_bounds):
    Check if the sequence falls in the range of GC-content bounds

    - is_in_length_bounds(seq_record, length_bounds):
    Check if the sequence falls in the range of length bounds

    - is_above_quality_threshold(seq_record, quality_threshold):
    Check if the mean of quality values exceeds the quality threshold

    Return:
    - file: file with fastq extension containing selected fragments.

    For more information please see README

    """
    if type(length_bounds) is int:
        length_bounds = 0, length_bounds
    if type(gc_bounds) is int or type(gc_bounds) is float:
        gc_bounds = 0, gc_bounds
    sequences = SeqIO.parse(input_path, "fastq")
    good_reads = []
    for seq_record in sequences:
        if (is_in_gc_bounds(seq_record, gc_bounds) and
            is_in_length_bounds(seq_record, length_bounds) and
            is_above_quality_threshold(seq_record, quality_threshold)):
            good_reads += [seq_record]
    if len(good_reads) == 0:
        raise ValueError('There are no sequences suited to requirements')
    if output_filename is None:
        input_filename = os.path.split(input_path)[-1]
        output_filename = input_filename
    if not (output_filename.endswith('.fastq')):
        output_filename = output_filename + '.fastq'
    current_directory = os.getcwd()
    path = os.path.join(current_directory, 'fastq_filtrator_results')
    if not (os.path.exists(path)):
        os.mkdir(path)
    output_path = os.path.join(path, output_filename)
    if os.path.exists(output_path):
        error = 'File with such name exists! Change output_filename arg!'
        raise ValueError(error)
    SeqIO.write(good_reads, output_path, "fastq")


class BiologicalSequence(ABC):
    """
    The abstract class for biological sequences
    """
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, item):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def __repr__(self):
        pass

    @abstractmethod
    def is_alphabet_correct(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleic acids
    """
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, item):
        return self.seq[item]

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.seq

    def is_alphabet_correct(self) -> bool:
        """
        The function check does the sequence contain
        standard A, G, C, T, U nucleotides

        Returns: bool - the result of check
        """
        return type(self).is_correct(self)

    def complement(self) -> BiologicalSequence:
        """
        Output the complementary sequence
        The complementarity rule could be found here:
        https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)

        Arguments:
        - self: the sequence obj to change

        Return: DNAsequence or RNAsequence - the result sequence object
        """
        complement = self.seq.translate(type(self).rule_complement)
        complement_seq = type(self)(complement)
        return complement_seq

    def gc_content(self) -> float:
        """
        Check the gc content in the sequence object
        Returns: int - the percentage of GC content
        """
        gc_counter = 0
        for nucl in self.seq:
            if nucl in ('G', 'C', 'g', 'c'):
                gc_counter += 1
        gc_share = gc_counter / len(self.seq) * 100
        return gc_share


class RNASequence(NucleicAcidSequence):
    """
    The class for RNA sequences
    """
    rule_complement = 'AUCG'.maketrans('AaUuCcGg', 'UuAaGgCc')

    def __init__(self, seq):
        super().__init__(seq)
        if not (super().is_alphabet_correct()):
            raise ValueError('The sequence does not correspond to RNA')

    def is_correct(self) -> bool:
        """
        Check if the sequence is RNA

        Arguments:
        - self: the sequence object to check

        Return:
        - bool: the result of the check
        """
        return set(self.seq) <= set('AaGgCcUu')


class DNASequence(NucleicAcidSequence):
    """
    The class for DNA sequences
    """
    rule_complement = 'ATCG'.maketrans('AaTtCcGg', 'TtAaGgCc')
    rule_transcription = 'AUCG'.maketrans('Tt', 'Uu')

    def __init__(self, seq):
        super().__init__(seq)
        if not (super().is_alphabet_correct()):
            raise ValueError('The sequence does not correspond to DNA')

    def is_correct(self) -> bool:
        """
        Check if the sequence is DNA

        Arguments:
        - self: the sequence object to check

        Return:
        - bool: the result of the check
        """
        return set(self.seq) <= set('AaGgCcTt')

    def transcribe(self) -> RNASequence:
        """
        Transcribe DNA sequence to RNA

        Arguments:
        - self: the sequence object to change

        Return: str: RNASequence object
        """
        transcribe_seq = self.seq.translate(self.rule_transcription)
        rna_seq = RNASequence(transcribe_seq)
        return rna_seq


class AminoAcidSequence(BiologicalSequence):
    """
    The class for amino acid sequences
    """
    aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY'

    def __init__(self, seq):
        self.seq = seq
        if not (self.is_alphabet_correct()):
            error = 'The sequence does not match standard protein code'
            raise ValueError(error)

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, item):
        return self.seq[item]

    def __str__(self):
        return self.seq

    def __repr__(self):
        return self.seq

    def is_alphabet_correct(self) -> bool:
        """
        The function checks if the sequence object contains
        standard amino acid code
        Returns: list of alternative frames (AminoAcidSequence)
        """
        return set(self.seq.upper()) <= set(self.aa_alphabet)

    def search_for_alt_frames(self, alt_start_aa='M') -> list:
        """
        Search for alternative frames in a protein sequences
        Use only one-letter code

        Search is not sensitive for letter case
        Without an alt_start_aa argument search for
        frames that start with methionine ('M')
        To search frames with alternative start codon
        add alt_start_aa argument

        The function ignores the last three amino acids in sequences

        Arguments:
        - self: sequence object to check
        - alt_start_aa (str): the name of an amino acid that
        is encoded by alternative start AA (Optional)
        Default: alt_start_aa='I'

        Return: the list of alternative frames
        """
        alternative_frames = []
        num_position = 0
        for amino_acid in self.seq[1:-3]:
            alt_frame = ''
            num_position += 1
            if amino_acid.upper() == alt_start_aa:
                alt_frame += self.seq[num_position:]
                alt_frame = AminoAcidSequence(alt_frame)
                alternative_frames.append(alt_frame)
        return alternative_frames
