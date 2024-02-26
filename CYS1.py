import os
import re
from typing import Tuple, List
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from abc import ABC


def run_dna_rna_tools(*args: List[str]) -> List[str]:
    """
    run_dna_rna_tools(*args: List[str])
    args: sequences of DNA or RNA
    last element in args: name of procedure
    available procedures:
        - transcribe
        - reverse
        - complement
        - reverse_complement
        - get_nucl_acid_type

    Returns:
        - list of sequences modified according to procedure
        - or string if only one sequence was given
    """
    procedure = args[-1]
    seqs = args[:-1]

    for seq in seqs:
        is_valid_dna_rna(seq)

    processed_result = []
    for seq in seqs:
        match procedure:
            case 'transcribe':
                processed_result.append(transcribe(seq))
            case 'reverse':
                processed_result.append(reverse(seq))
            case 'complement':
                processed_result.append(complement(seq))
            case 'reverse_complement':
                processed_result.append(reverse_complement(seq))
            case 'get_nucl_acid_type':
                processed_result.append(get_nucl_acid_type(seq))
            case _:
                return 'Procedure is not defined in the function'

    if len(processed_result) == 1:
        return processed_result[0]
    return processed_result


def protein_analysis(*args: str, procedure: str, cell_type: str = None,
                    letter_format: int = 1) -> list:
    """
    Function protein_analysis:
    - calculates predicted molecular weight of amino acid sequences in kDa
    (procedure name: molecular_weight)
    - translate aa sequences from one-letter to three-letter code
    (procedure name: one_letter_to_three)
    - calculates total amount of each amino acid in the sequences
    (procedure name: get_amino_acid_sum)
    - makes DNA based codon optimization for
    the introduced amino acid sequences, support 3 types of cells:
      Esherichia coli, Pichia pastoris, Mouse
      (procedure name: codon_optimization)
    - calculates length of amino acid sequences
    (procedure name: length)
    - counts the number of atoms of each type in a sequence
    (procedure name: brutto_count)

    Arguments:
    - one or multiple string of protein sequences written one letter
    or three letter code (not mixed)
    - name of procedure as string
    - cell type (required only for codon_optimization procedure)
    - letter_format of code for the protein sequences as int: 1 for one letter,
     3 for three letter code

    Return:
    - molecular_weight procedure returns list of floats
    - one_letter_to_three procedure returns list of strings
    - get_amino_acid_sum procedure returns list of dictionaries
    - codon_optimization procedure returns list of strings
    - length procedure returns list of int values
    - brutto_count procedure returns list of dictionaries with counts
    of atoms in the sequence
    """
    amino_acid_seqs = name_transform(args, letter_format)
    procedures = {
        "molecular_weight": molecular_weight,
        "one_letter_to_three": one_letter_to_three,
        "get_amino_acid_sum": get_amino_acid_sum,
        "codon_optimization": codon_optimization,
        "length": length,
        "brutto_count": brutto_count,
    }
    if procedure not in procedures.keys():
        raise ValueError("Requested procedure is not defined")
    elif procedure == "codon_optimization":
        return procedures.get(procedure)(amino_acid_seqs, cell_type)
    else:
        return procedures.get(procedure)(amino_acid_seqs)


def filter_fastq(input_path: str, output_filename: str = '',
                 gc_bounds: Tuple[int, int] = (0, 100),
                 length_bounds: Tuple[int, int] = (0, 2 ** 32),
                 quality_threshold: int = 0):
    """
        filter_fastq(seqs: Dict[str, Tuple[str, str]],
                           gc_bounds: Tuple[int, int] = (0, 100),
                           length_bounds: Tuple[int, int] = (0, 2 ** 32),
                           quality_threshold: int = 0)
        :param input_path: path to fastq file
        :param output_filename: desired name for the filered fasta file
        (if not given _output is added to input file name)
        :param gc_bounds: borders of GC-content that will be used
        to filter sequence
        :param length_bounds: borders of sequence length that will
        be used to filter sequence
        :param quality_threshold: borders of quality that will
        be used to filter sequence (mean quality of the sequence
        is considered to pass the treshold)
        All the borders include upper and lower values.

        :return: dictionary of the same structure as input,
        that only contains sequences that passed all the
        thresholdings
        """
    dir_name = os.path.dirname(input_path)
    if output_filename == '':
        output_filename = re.sub(r'\.[^.]*$', '',
                                os.path.basename(input_path)) + "_output"

    initial_sequences = []

    with open(input_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fastq"):
            initial_sequences.append(record)

    if isinstance(gc_bounds, int):
        gc_lower_bound = 0
        gc_upper_bound = gc_bounds
    else:
        gc_lower_bound = gc_bounds[0]
        gc_upper_bound = gc_bounds[1]

    if isinstance(length_bounds, int):
        length_lower_bound = 0
        length_upper_bound = length_bounds
    else:
        length_lower_bound = length_bounds[0]
        length_upper_bound = length_bounds[1]

    if not os.path.exists(os.path.join(dir_name, 'fastq_filtrator_resuls')):
        os.makedirs(os.path.join(dir_name, 'fastq_filtrator_resuls'))

    try:
        with open(os.path.join(dir_name, 'fastq_filtrator_resuls',
                               f'{output_filename}.fasta'), mode='x') as file:
            for seq in initial_sequences:
                if is_in_gc_bounds(seq, gc_lower_bound, gc_upper_bound) and \
                        is_in_length_bounds(seq, length_lower_bound,
                                            length_upper_bound) and \
                        is_above_quality_threshold(seq, quality_threshold):
                    file.write(seq.id+'\n')
                    file.write(str(seq.seq) + '\n')
    except FileExistsError:
        print('File with the provided name exist. Please use another name.')


def is_in_gc_bounds(seq, gc_lower_bound: int, gc_upper_bound: int) -> bool:
    if gc_lower_bound <= gc_fraction(seq.seq) * 100 <= gc_upper_bound:
        return True
    return False


def is_in_length_bounds(seq, length_lower_bound: int, length_upper_bound: int):
    if length_lower_bound <= len(seq) <= length_upper_bound:
        return True
    return False


def is_above_quality_threshold(seq, quality_threshold):
    total_q_score = sum(seq.letter_annotations["phred_quality"])
    if total_q_score/len(seq) >= quality_threshold:
        return True
    return False


AMINO_SHORT_NAMES_DIC = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "V": "Val",
    "H": "His",
    "G": "Gly",
    "Q": "Gln",
    "E": "Glu",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "P": "Pro",
    "S": "Ser",
    "Y": "Tyr",
    "T": "Thr",
    "W": "Trp",
    "F": "Phe",
    "C": "Cys",
}

AMINO_NAMES_DIC = {
    "ala": "A",
    "arg": "R",
    "asn": "N",
    "asp": "D",
    "val": "V",
    "his": "H",
    "gly": "G",
    "gln": "Q",
    "glu": "E",
    "ile": "I",
    "leu": "L",
    "lys": "K",
    "met": "M",
    "pro": "P",
    "ser": "S",
    "tyr": "Y",
    "thr": "T",
    "trp": "W",
    "phe": "F",
    "cys": "C",
}

AMINO_NAMES_DIC_REVERSE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Val": "V",
    "His": "H",
    "Gly": "G",
    "Gln": "Q",
    "Glu": "E",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Pro": "P",
    "Ser": "S",
    "Tyr": "Y",
    "Thr": "T",
    "Trp": "W",
    "Phe": "F",
    "Cys": "C",
}

amino_weights = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.16,
    "E": 147.13,
    "Q": 146.15,
    "G": 75.07,
    "H": 155.16,
    "I": 131.18,
    "L": 131.18,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

amino_brutto = {
    "A": (3, 7, 1, 2, 0),
    "R": (6, 14, 4, 2, 0),
    "N": (4, 8, 2, 3, 0),
    "D": (4, 7, 1, 4, 0),
    "V": (5, 11, 1, 2, 0),
    "H": (6, 9, 3, 2, 0),
    "G": (2, 5, 1, 2, 0),
    "Q": (5, 10, 2, 3, 0),
    "E": (5, 9, 1, 4, 0),
    "I": (6, 13, 1, 2, 0),
    "L": (6, 13, 1, 2, 0),
    "K": (6, 14, 2, 2, 0),
    "M": (5, 11, 1, 2, 1),
    "P": (5, 9, 1, 2, 0),
    "S": (3, 7, 1, 3, 0),
    "Y": (9, 11, 1, 3, 0),
    "T": (4, 9, 11, 1, 3, 0),
    "W": (11, 12, 2, 2, 0),
    "F": (9, 11, 1, 2, 0),
    "C": (3, 7, 1, 2, 1),
}

ecoli_triplets = {
    "A": "GCG",
    "C": "TGC",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGC",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCG",
    "Q": "CAG",
    "R": "CGT",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAT",
}

ppastoris_triplets = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAG",
    "L": "TTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCA",
    "Q": "CAA",
    "R": "AGA",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAC",
}

mouse_triplets = {
    "A": "GCC",
    "C": "TGC",
    "D": "GAC",
    "E": "GAG",
    "F": "TTC",
    "G": "GGC",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCC",
    "Q": "CAG",
    "R": "CGG",
    "S": "AGC",
    "T": "ACC",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAC",
}


def molecular_weight(amino_acid_seqs: list) -> float:
    """
    Calculates predicated molecular weight of aa sequences.

     Arguments:
    - amino_acid_seqs (list): list of string with the protein sequences

    Return:
    - List of floats corresponding to the molecular weight in kDa
    """
    molecular_weights = []
    for seq in amino_acid_seqs:
        total_weight = 0
        for aa in seq:
            aa = aa.upper()
            total_weight += amino_weights[aa]
        molecular_weights.append(round(total_weight / 1000, 2))
    return sum(molecular_weights)


def one_letter_to_three(amino_acid_seqs: list) -> list:
    """
    Translates one-letter coded amino acid sequences to three-letter coded
    Arguments:
    - amino_acid_seqs (list): list of string with the protein sequences

    Return:
    - List of of strings with three-letter coded sequences
    """
    three_letters_seqs = []
    for seq in amino_acid_seqs:
        three_letters_seq = []
        for aa in seq:
            aa = aa.upper()
            three_letters_seq.append(AMINO_SHORT_NAMES_DIC[aa])
        three_letters_seqs.append("".join(three_letters_seq))
    return three_letters_seqs


def get_amino_acid_sum(protein_sequences: list) -> list:
    """
    Counts the amount of each amino acid in the injected protein sequences

    Arguments:
    - protein_sequences (list): list of injected protein sequence

    Return:
    - List of dictionary with amino acid amount"""
    result = []
    for protein_sequence in range(len(protein_sequences)):
        d_amino_counts = {}
        for key in AMINO_SHORT_NAMES_DIC.keys():
            d_amino_counts[key] = 0
        for amino_acid in protein_sequences[protein_sequence]:
            d_amino_counts[amino_acid] += 1
        result.append(d_amino_counts)
    return result


def codon_optimization(protein_sequences: list, cell_type: str) -> list:
    """
    Makes codon-optimized DNA based on the introduced amino acid sequences
    for 3 types of cells:
    Esherichia coli, Pichia pastoris, Mouse

    Arguments:
    - protein_sequences (list): list of injected protein sequence
    - cell_type (str): user-entered cell type for codon optimization

    Return:
    - List of codon-optimized DNA"""

    if cell_type == "Esherichia coli" or cell_type == "E.coli":
        codon_optimization_ecoli = []
        replacer_ecoli = ecoli_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_ecoli += [
                "".join([replacer_ecoli(n, n) \
                         for n in protein_sequences[amino_acid]])
            ]
        return codon_optimization_ecoli

    if cell_type == "Pichia pastoris" or cell_type == "P.pastoris":
        codon_optimization_ppastoris = []
        replacer_ppastoris = ppastoris_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_ppastoris += [
                "".join(
                    [replacer_ppastoris(n, n)  \
                     for n in protein_sequences[amino_acid]]
                )
            ]
        return codon_optimization_ppastoris

    if cell_type == "Mouse" or cell_type == "mouse":
        codon_optimization_mouse = []
        replacer_mouse = mouse_triplets.get
        for amino_acid in range(len(protein_sequences)):
            codon_optimization_mouse += [
                "".join([replacer_mouse(n, n) \
                         for n in protein_sequences[amino_acid]])
            ]
        return codon_optimization_mouse
    else:
        raise ValueError(
            f'Type {cell_type} is not supported. \
            The following types of organisms are \
            available for codon optimization: \
            Esherichia coli, Pichia pastoris, Mouse'
        )


def length(seqs: list) -> list:
    """
    Counts total length of amino acid sequence.

    Arguments:
    - seqs (list): list of string with the protein sequences

    Return:
    - list of int values corresponding to the length of sequences"""
    result = [len(seq) for seq in seqs]
    return result


def name_transform(seqs: tuple, letter_format: int) -> list:
    """
    Transforms the amino acid sequences given to protein_analysis function
    from three-letter code to one-letter code,
    makes sequences unified (for one-letter letter_format all letters
    to upper and for three-letter letter_format to lower).

    Arguments:
      - seqs (tuple): tuple of string with the protein sequences

      Return:
      - list of strings with the transformed sequences"""
    result = []
    multiple_of_three = []
    test_three_letters = []
    if letter_format == 1:
        for seq in seqs:
            multiple_of_three.append(is_length_divisible_by_3(seq))
            test_three_letters.append(is_amino_acid_three_letter(seq))
            seq = seq.upper()
            for letter in seq:
                if is_amino_acid(letter):
                    pass
            result.append(seq)
        if all(multiple_of_three) and all(test_three_letters):
            print(
                "Warning: all your sequences are similar to three-letter ones.\
                Check the letter_format value"
            )
        return result
    elif letter_format == 3:
        for seq in seqs:
            seq = seq.lower()
            seq3 = [seq[i: i + 3] for i in range(0, len(seq), 3)]
            for triplet in seq3:
                if is_amino_acid(triplet):
                    pass
            seq_transformed = "".join([AMINO_NAMES_DIC.get(seq) \
                                       for seq in seq3])
            result.append(seq_transformed)
        return result
    else:
        raise ValueError(
            "Error unsupported letter_format. Only letter_formats 1 \
             and 3 are supported"
        )


def is_amino_acid(input_amino: str) -> bool:
    """
    Checks whether the entered string is an amino acid (either three-letter
    encoding or one-letter encoded).

    Arguments:
      - input_amino (str): string corresponding to one amino acid
      (in three-letter code or one-letter code)

      Return:
      - bool: True if amino acid is a valid amino acid, otherwise
      ValueError is amino acid is not correct
    """
    if len(input_amino) == 1:
        letter = input_amino
        if letter not in AMINO_SHORT_NAMES_DIC.keys():
            raise ValueError(f"Error {letter} \
            is not an amino acid. Correct your input")
        return True
    elif len(input_amino) == 3:
        triplet = input_amino
        if triplet not in AMINO_NAMES_DIC.keys():
            raise ValueError(
                f"Error {triplet} is not an amino acid. Correct your input"
            )
        return True
    else:
        raise ValueError(
            f"Error {input_amino} is incorrect form \
            of amino acid notation. Correct your input"
        )


def brutto_count(seqs: list) -> list:
    """
    Calculates the brutto formula of the amino acid sequences.

    Arguments:
      - seqs (list): list of string with the protein sequences

      Return:
      - list of dictionaries with counts of each elemet included
      (elements C,H,N,O,S)"""
    elements = ["C", "H", "N", "O", "S"]
    result = []
    for seq in seqs:
        brutto_list = [amino_brutto.get(letter) for letter in seq]
        brutto_pair = list(zip(*brutto_list))
        brutto = [sum(i) for i in brutto_pair]
        brutto_dict = dict(zip(elements, brutto))
        result.append(brutto_dict)
    return result


def is_length_divisible_by_3(seq: str) -> bool:
    """
    Checks if the sequence is divisible by three.

    Arguments:
      - seq (str): string of protein sequence

      Return:
      - bool: True if sequence is divisible by three, otherwise False"""
    seq_len = len(seq)
    if seq_len % 3 == 0:
        return True
    else:
        return False


def is_amino_acid_three_letter(seq: str) -> bool:
    """
    Checks whether all elements of a sequence are three-letter
    amino acid symbols.

    Arguments:
      - seq (str): string of protein sequence

    Return:
      - bool: True if sequence is corresponding to the valid
      three-letter amino acid, otherwise False
    """
    seq = seq.lower()
    seq3 = [seq[i: i + 3] for i in range(0, len(seq), 3)]
    for triplet in seq3:
        if triplet not in AMINO_NAMES_DIC.keys():
            return False
        else:
            return True


def is_valid_dna_rna(seq: str) -> bool:
    set_of_seq = set(seq.upper())
    if set_of_seq > {'A', 'T', 'G', 'C', 'U'}:
        raise ValueError('Invalid alphabet')
    elif ('T' in set_of_seq) and ('U' in set_of_seq):
        raise ValueError('Invalid alphabet')
    else:
        return True


def transcribe(seq: str) -> str:
    transcribed_seq = ''
    for nucleotide in seq:
        if nucleotide == 'T':
            transcribed_seq += 'U'
        elif nucleotide == 't':
            transcribed_seq += 'u'
        else:
            transcribed_seq += nucleotide
    return transcribed_seq


def reverse(seq: str) -> str:
    reversed_seq = seq[::-1]
    return reversed_seq


def complement(seq: str) -> str:
    complement_seq = ''
    dna_complement_map = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'a': 't',
        't': 'a',
        'g': 'c',
        'c': 'g'
    }
    for nucleotide in seq:
        complement_seq += dna_complement_map[nucleotide]
    return complement_seq


def reverse_complement(seq: str) -> str:
    return complement(reverse(seq))


def get_nucl_acid_type(seq: str) -> str:
    set_of_seq = set(seq)
    if ('U' in set_of_seq) or ('u' in set_of_seq):
        return 'RNA'
    elif ('T' in set_of_seq) or ('t' in set_of_seq):
        return 'DNA'
    return 'ND'


class BiologicalSequence(ABC):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self) -> int:
        """Returns the length of the sequence."""
        return len(self.sequence)

    def __getitem__(self, index) -> str:
        """Returns the item pointed by index or a slice of the sequence."""
        return self.sequence[index]

    def __str__(self) -> str:
        """Returns a string representation of the sequence."""
        return str(self.sequence)

    def check_alphabet(self, alphabet: str) -> bool:
        """Checks if the sequence consists of allowed symbols."""
        flag = True
        for letter in self.sequence:
            if letter not in alphabet:
                flag = False
        return flag


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        super().__init__(sequence)

    def check_alphabet(self) -> bool:
        return super().check_alphabet("ATCGUN")

    def complement(self) -> str:
        return complement(self.sequence)

    def gc_content(self) -> float:
        return gc_fraction(self.sequence)


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        return RNASequence(transcribe(self.sequence))


class RNASequence(NucleicAcidSequence):
    pass


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        super().__init__(sequence)

    def check_alphabet(self) -> bool:
        return super().check_alphabet("ACDEFGHIKLMNPQRSTVWY")

    def molecular_weight(self):
        return molecular_weight(self.sequence)
