import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio import SeqUtils


class BiologicalSequence(ABC, str):

    def __init__(self, seq):
        self.seq = seq

    #да, Никита уже сказал, что так лучше не делать, а оставить сдесь только абстрактные методы,
    # но я зацепилась за идею не повторять проверки, раз суть для всех одна, только алфавиты разные,
    # а переделать уже не успеваю
    def check_alphabet(self):
        if set(self.seq) <= self.alphabet:
            raise UnexpectedSymbolInSeqError

    # не поняла, не доделала
    @abstractmethod
    def to_print_seq(self):
        return self.seq
        pass


class UnexpectedSymbolInSeqError(ValueError):
    pass


class AminoAcidSequence(BiologicalSequence):
    alphabet = {'M', 'O', 'v', 'D', 'f', 'N', 'c', 'A', 'R', 'W', 'I', 'm', 'L', 's', 'H', 'q', 'w', 'V', 'n', 'i',
                   'g', 'F', 'S', 'e', 'l', 'U', 'P', 'Q', 'K', 'Y', 'u', 'y', 'd', 'h', 'k', 'r', 't', 'G', 'o', 'E',
                   'p', 'T', 'C', 'a'}
    masses = {'A': 71.08, 'R': 156.2, 'N': 114.1, 'D': 115.1, 'C': 103.1, 'E': 129.1, 'Q': 128.1, 'G': 57.05, 'H': 137.1,
              'I': 113.2, 'L': 113.2, 'K': 128.2, 'M': 131.2, 'F': 147.2, 'P': 97.12, 'S': 87.08, 'T': 101.1,
              'W': 186.2, 'Y': 163.2, 'V': 99.13, 'U': 168.05, 'O': 255.3, 'a': 71.08, 'r': 156.2, 'n': 114.1, 'd': 115.1,
              'c': 103.1, 'e': 129.1, 'q': 128.1, 'g': 57.05, 'h': 137.1, 'i': 113.2, 'l': 113.2, 'k': 128.2, 'm': 131.2,
              'f': 147.2, 'p': 97.12, 's': 87.08, 't': 101.1, 'w': 186.2, 'y': 163.2, 'v': 99.13, 'u': 168.05, 'o': 255.3}

    def __init__(self, seq):
        super().__init__(seq)
        self.seq = seq
        self.check_alphabet()
    def molecular_weight(self) -> float:
        """
        Function calculates molecular weight of the amino acid chain
        Input Parameters:
               each letter refers to one-letter coded proteinogenic amino acids
        Returns:
            (float) Molecular weight of tge given amino acid chain in Da
        """
        m = 0
        for acid in str(self.seq):
            m += self.masses[acid]
        return m

class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        super().__init__(seq)
        self.check_alphabet()

    def gc_content(self) -> float:
        """
        Function counts GC-content in sequence, and returns result in %
        """
        n = 0
        for nucl in self.seq:
            if nucl == 'c' or nucl == 'g' or nucl == 'C' or nucl == 'G':
                n += 1
        return 100 * n / len(self.seq)

    def complement(self):
        """
        Function return complement sequence
        """
        complementary_dna = self.seq.maketrans(self.dict_comp)
        res = self.seq.translate(complementary_dna)
        return NucleicAcidSequence(res)

class DNASequence(BiologicalSequence):
    alphabet = {'A', 'g', 't', 'G', 'T', 'a', 'c', 'C'}
    dict_trans = {'A': 'A', 'C': 'C', 'T': 'U', 'G': 'G', 'a': 'a', 'c': 'c', 't': 'u', 'g': 'g'}
    dict_comp = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
    def __init__(self, seq):
        super().__init__(seq)


    def transcribe(self):
        """
        Function return transcript of DNA sequence
        """
        transcribe = self.seq.maketrans(self.dict_trans)
        res = self.seq.translate(transcribe)

        return RNASequence(res)


class RNASequence(BiologicalSequence):
    alphabet = {'U', 'A', 'g', 'G', 'a', 'c', 'C', 'u'}
    dict_comp = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C', 'a': 'u', 'c': 'g', 'u': 'a', 'g': 'c'}
    def __init__(self, seq):
        super().__init__(seq)


 #фильтратор
def filter_gc(records, gc_bounds_both_side=(0, 100)) -> list:
    """
    This function selects sequences with the GC content of your interest
    :parameters:
        records: records from fastq parced by SeqIO
        gc_bound: interval for the of acceptable GC content, in %
    :return:(dict) new dictionary consists of selected sequences
    """
    new_records = []
    for record in records:
        if (gc_bounds_both_side[1]/100  >= SeqUtils.gc_fraction(record.seq) >= gc_bounds_both_side[0]/100):
            new_records.append(record)

    return new_records


def filter_length(records, length_bounds_both_side=(0, 2 ** 32)) -> list:
    """
    This function selects sequences with the length of your interest
    :parameters:
        records: records from fastq parced by SeqIO
        length_bound: interval for the of acceptable sequense length in number of nucleotide
    :return:(dict) new dictionary consists of selected sequences
    """
    new_records = []
    for record in records:
        if (length_bounds_both_side[1] >= len(record.seq) >= length_bounds_both_side[0]):
            new_records.append(record)
    print(new_records)
    return new_records



def filter_quality(records, quality_threshold=0) -> list:
    """
    This function selects  FASTQ sequences with appropriate average nucleotide read quality
    parameters:
        seqs: dictionary of FASTQ sequences {name: (sequence, quality)}
        quality_treshold: threshold value for average quality per nucleotide (phred33 scale)
    :return:(dict) recordes for selected sequences
    """
    new_records = []
    for record in records:
        if (sum(record.letter_annotations["phred_quality"])/len(record.seq) >= quality_threshold):
            new_records.append(record)
    print(new_records)
    return new_records


def fastq_filtration(input_fastq, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_treshold=0, output_fastq=''):
    """
    This function provides you the opportunity to filter the FASTQ file to select sequences
    that fit  requirements on 5 parameters: input and output(optional) files, length, GC composition,
    and quality of the reed
    :parameters
        input_fastq: path to fastq file
        gc_bounds: (tuple) interval for the of acceptable GC content, in %
        length_bounds: (tuple) interval for the of acceptable sequense length in number of nucleotide
        quality_treshold: (float) threshold value for average quality per nucleotide (phred33 scale)
        output_fastq = name of output file, ./fastq_filtrator_resuls/output_fastq, if it is not defined,
        it will be the same of the input file

    """
    if not os.path.isdir('fastq_filtrator_resuls'):
        os.mkdir('fastq_filtrator_resuls')
    if output_fastq == '':
        output_fastq = os.path.join('fastq_filtrator_resuls', os.path.basename(input_fastq))
    else:
        output_fastq = os.path.join('fastq_filtrator_resuls', output_fastq + ".fasta")
    if type(gc_bounds) == float or type(gc_bounds) == int:
        gc_bounds_both_side = (0, gc_bounds)
    else:
        gc_bounds_both_side = gc_bounds
    if type(length_bounds) == int:
        length_bounds_both_side = (0, length_bounds)
    else:
        length_bounds_both_side = length_bounds
    records = list(SeqIO.parse(input_fastq, "fastq"))
    filtered_records_l = filter_length(records, length_bounds_both_side)
    filtered_records_gc = filter_gc(filtered_records_l, gc_bounds_both_side)
    filtered_records_q = filter_quality(filtered_records_gc, quality_treshold)
    SeqIO.write((record for record in filtered_records_q), output_fastq, "fastq")
    return

#fastq_filtration('example_fastq.fastq', length_bounds= (10, 30), gc_bounds=(20, 30), output_fastq='')
