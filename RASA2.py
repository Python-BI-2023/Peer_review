import os.path
from typing import Tuple, Optional, Union, List, Any
from Bio import SeqIO
from statistics import mean
from Bio import SeqUtils
from abc import ABC, abstractmethod


def filter_fastq(input_path: str,
                 gc_bounds: int | float | Tuple = (20, 80),
                 length_bounds: int | float | Tuple = (0, 2 ** 32),
                 quality_threshold: int = 0, output_filename='') -> None:
    seqs = SeqIO.parse(input_path, 'fastq')
    if seqs is None:
        raise ValueError('Your fastq_files are None')

    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)

    if output_filename == '':
        output_filename = os.path.basename(input_path)

    if '.fastq' in output_filename:
        output_filename = output_filename.replace('.fastq', '')

    output_dir = os.path.join(os.path.dirname(input_path),
                              'filter_fastq_results')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(os.path.join(output_dir, f'{output_filename}.fastq'), 'w') as fq:
        for record in seqs:
            quality = mean(record.letter_annotations["phred_quality"])
            if quality >= quality_threshold and \
                    length_bounds[0] <= len(record.seq) <= \
                    length_bounds[1]:
                gc_content = SeqUtils.GC123(record.seq)
                if gc_bounds[1] >= gc_content[0] >= gc_bounds[0]:
                    SeqIO.write(record, fq, 'fastq')


class BiologicalSequence:
    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, item) -> 'BiologicalSequence':
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def print_seq(self) -> 'BiologicalSequence':
        pass

    @abstractmethod
    def is_valid_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def print_seq(self):
        print(str(self.sequence))

    def is_valid_alphabet(self):
        valid_nucleotides = {'A', 'C', 'G', 'T', 'U'}
        return all(
            nucleotide in valid_nucleotides for nucleotide in self.sequence)

    def complement(self) -> 'NucleicAcidSequence':
        return self.exact_complement()

    def exact_complement(self) -> 'NucleicAcidSequence':
        raise NotImplemented

    def count_gc_content(self):
        c_count = self.sequence.count('C')
        g_count = self.sequence.count('G')
        gc_content = (c_count + g_count) / len(self.sequence)
        return gc_content

    def reverse(self) -> 'NucleicAcidSequence':
        raise NotImplemented


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence)

    COMPLEMENT_DNA = {
        'A': 'T',
        'a': 't',
        'G': 'C',
        'g': 'c',

        'T': 'A',
        't': 'a',
        'C': 'G',
        'c': 'g',
    }

    def transcribe(self) -> 'DNASequence':
        transcribe_sequence = ''
        for nucl in self.sequence:
            if nucl.upper() == 'A':
                transcribe_sequence += nucl
            elif nucl.upper() == 'G':
                transcribe_sequence += nucl
            elif nucl.upper() == 'C':
                transcribe_sequence += nucl
            elif nucl.upper() == 'T':
                if nucl == 't':
                    transcribe_sequence += 'u'
                elif nucl == 'T':
                    transcribe_sequence += 'U'

        return DNASequence(transcribe_sequence)

    def exact_complement(self) -> 'NucleicAcidSequence':
        complemented = ''
        for nucl in self.sequence:
            complemented += self.COMPLEMENT_DNA[nucl]
        return DNASequence(complemented)

    def reverse(self) -> 'NucleicAcidSequence':
        return DNASequence(self.sequence[::-1])


class RNASequence(NucleicAcidSequence):
    COMPLEMENT_RNA = {
        'A': 'U',
        'a': 'u',
        'G': 'C',
        'g': 'c',

        'U': 'A',
        'u': 'a',
        'C': 'G',
        'c': 'g',
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def exact_complement(self) -> 'NucleicAcidSequence':
        complemented = ''
        for nucl in self.sequence:
            complemented += self.COMPLEMENT_RNA[nucl]
        return RNASequence(complemented)

    def reverse(self) -> 'NucleicAcidSequence':
        return RNASequence(self.sequence[::-1])


class AminoAcidSequence(BiologicalSequence):
    AMINOACIDS_DICT = {
        'Ala': {'TO_1': 'A',
                'PROTEIN_TO_RNA_COMBINATION': {'GCU', 'GCC', 'GCA', 'GCG'},
                'PKA_AMINOACIDS': [2.34, 9.69],
                'MOLECULAR_WEIGHTS': 89},
        'Arg': {'TO_1': 'R',
                'PROTEIN_TO_RNA_COMBINATION': {'CGU', 'CGC', 'CGA', 'CGG',
                                               'AGA',
                                               'AGG'},
                'PKA_AMINOACIDS': [2.17, 9.04, 12.68],
                'MOLECULAR_WEIGHTS': 174},
        'Asn': {'TO_1': 'N',
                'PROTEIN_TO_RNA_COMBINATION': {'AAU', 'AAC'},
                'PKA_AMINOACIDS': [1.88, 9.60, 3.65],
                'MOLECULAR_WEIGHTS': 132},
        'Asp': {'TO_1': 'D',
                'PROTEIN_TO_RNA_COMBINATION': {'GAU', 'GAC'},
                'PKA_AMINOACIDS': [1.88, 9.60, 3.65],
                'MOLECULAR_WEIGHTS': 133},
        'Cys': {'TO_1': 'C',
                'PROTEIN_TO_RNA_COMBINATION': {'UGU', 'UGC'},
                'PKA_AMINOACIDS': [1.96, 10.28, 8.18],
                'MOLECULAR_WEIGHTS': 121},
        'Glu': {'TO_1': 'Q',
                'PROTEIN_TO_RNA_COMBINATION': {'GAA', 'GAG'},
                'PKA_AMINOACIDS': [2.19, 9.67, 4.25],
                'MOLECULAR_WEIGHTS': 147},
        'Gln': {'TO_1': 'E',
                'PROTEIN_TO_RNA_COMBINATION': {'CAA', 'CAG'},
                'PKA_AMINOACIDS': [2.17, 9.13],
                'MOLECULAR_WEIGHTS': 146},
        'Gly': {'TO_1': 'G',
                'PROTEIN_TO_RNA_COMBINATION': {'GGU', 'GGC', 'GGA', 'GGG'},
                'PKA_AMINOACIDS': [2.34, 9.60],
                'MOLECULAR_WEIGHTS': 75},
        'His': {'TO_1': 'E',
                'PROTEIN_TO_RNA_COMBINATION': {'CAU', 'CAC'},
                'PKA_AMINOACIDS': [1.82, 9.17],
                'MOLECULAR_WEIGHTS': 155},
        'Ile': {'TO_1': 'I',
                'PROTEIN_TO_RNA_COMBINATION': {'AUU', 'AUC', 'AUA'},
                'PKA_AMINOACIDS': [2.36, 9.68],
                'MOLECULAR_WEIGHTS': 131},
        'Leu': {'TO_1': 'L',
                'PROTEIN_TO_RNA_COMBINATION': {'CUU', 'CUC', 'CUA', 'CUG'},
                'PKA_AMINOACIDS': [2.36, 9.60],
                'MOLECULAR_WEIGHTS': 131},
        'Lys': {'TO_1': 'K',
                'PROTEIN_TO_RNA_COMBINATION': {'AAA', 'AAG'},
                'PKA_AMINOACIDS': [2.18, 8.95, 10.53],
                'MOLECULAR_WEIGHTS': 146},
        'Met': {'TO_1': 'M',
                'PROTEIN_TO_RNA_COMBINATION': {'AUG'},
                'PKA_AMINOACIDS': [2.28, 9.21],
                'MOLECULAR_WEIGHTS': 149},
        'Phe': {'TO_1': 'F',
                'PROTEIN_TO_RNA_COMBINATION': {'UUU', 'UUC'},
                'PKA_AMINOACIDS': [2.20, 9.13],
                'MOLECULAR_WEIGHTS': 165},
        'Pro': {'TO_1': 'P',
                'PROTEIN_TO_RNA_COMBINATION': {'CCU', 'CCC', 'CCA', 'CCG'},
                'PKA_AMINOACIDS': [1.99, 10.96],
                'MOLECULAR_WEIGHTS': 115},
        'Ser': {'TO_1': 'S',
                'PROTEIN_TO_RNA_COMBINATION': {'UCU', 'UCC', 'UCA', 'UCG'},
                'PKA_AMINOACIDS': [2.21, 9.15],
                'MOLECULAR_WEIGHTS': 105},
        'Thr': {'TO_1': 'T',
                'PROTEIN_TO_RNA_COMBINATION': {'ACU', 'ACC', 'ACA', 'ACG'},
                'PKA_AMINOACIDS': [2.11, 9.62],
                'MOLECULAR_WEIGHTS': 119},
        'Tyr': {'TO_1': 'W',
                'PROTEIN_TO_RNA_COMBINATION': {'UAU', 'UAC'},
                'PKA_AMINOACIDS': [2.20, 9.11, 10.07],
                'MOLECULAR_WEIGHTS': 181},
        'Trp': {'TO_1': 'Y',
                'PROTEIN_TO_RNA_COMBINATION': {'UGG'},
                'PKA_AMINOACIDS': [2.38, 9.39],
                'MOLECULAR_WEIGHTS': 204},
        'Val': {'TO_1': 'V',
                'PROTEIN_TO_RNA_COMBINATION': {'GUU', 'GUC', 'GUA', 'GUG'},
                'PKA_AMINOACIDS': [2.32, 9.62],
                'MOLECULAR_WEIGHTS': 117},
    }

    TRANSCRIBE_DICT: dict = {'A': 'A',
                             'U': 'T',
                             'G': 'G',
                             'C': 'C',
                             'a': 'a',
                             'u': 't',
                             'g': 'g',
                             'c': 'c'}

    TO_3_DICT = {nested_dict['TO_1']: key for key,
    nested_dict in AMINOACIDS_DICT.items()}

    def is_valid_alphabet(self) -> bool:
        for letter in self.sequence:
            if letter not in self.TRANSCRIBE_DICT.keys():
                return False
        return True

    def __init__(self, sequence: str):
        self.sequence = sequence
        self.sequence_converted = ''
        self.check_input()

    def __len__(self) -> int:
        return len(self.sequence_converted)

    def __getitem__(self, index) -> 'AminoAcidSequence':
        return AminoAcidSequence(self.sequence_converted[index])

    def __str__(self) -> str:
        return self.sequence_converted

    def print_seq(self):
        print(str(self.sequence_converted))

    def is_one_letter(self, seq: str) -> bool:
        """
        Defines whether the sequence is 1 coded.

        Args:
        - seq - sequence to check

        Returns:
        - bool
        """
        return all(aa.isalpha() and aa.isupper() for aa in seq)

    def convert_aa_coding(self, seq: str) -> str:
        """
        Translate 1-letter to 3-letter encoding if 1-letter
        encoded sequence is given and vice versa.

        Args:
        - seq - sequence or list of sequences to recode

        Returns:
        - function_result - a dictionary containing recoded sequences as values
        for original sequences keys
        """

        if self.is_one_letter(seq):
            three_letter_sequence = ""
            for aa in seq:
                three_letter_code = self.TO_3_DICT.get(aa, aa)
                three_letter_sequence += three_letter_code
            return three_letter_sequence
        one_letter_sequence = ""
        for aa in range(0, len(seq), 3):
            amino_acid = seq[aa:aa + 3]
            one_letter_sequence += self.AMINOACIDS_DICT[amino_acid]['TO_1']
        return one_letter_sequence

    def check_input(self):
        if self.is_one_letter(self.sequence):
            self.sequence_converted = self.convert_aa_coding(self.sequence)
        else:
            self.sequence_converted = self.sequence

    def isoelectric_point_calculating(self) -> float:
        divided_acids = [self.sequence_converted[i:i + 3] for i in range(0,
                                                            len(self),
                                                            3)]
        for divided_acid in divided_acids:
            if divided_acid not in self.AMINOACIDS_DICT.keys():
                raise ValueError('Non-protein aminoacids in sequence')

        isoelectric_point = 0
        count_groups = 0
        for acid_index, aminoacid in enumerate(divided_acids):
            if acid_index == 0:
                isoelectric_point \
                    += (self.AMINOACIDS_DICT[aminoacid]['PKA_AMINOACIDS'][0])
                count_groups += 1
            elif acid_index == len(divided_acids) - 1:
                isoelectric_point = (isoelectric_point
                                     + (self.AMINOACIDS_DICT[aminoacid]
                        ['PKA_AMINOACIDS'][-1]))
                count_groups += 1
            else:
                if len(self.AMINOACIDS_DICT[aminoacid][
                           'PKA_AMINOACIDS']) > 2:
                    isoelectric_point = (isoelectric_point
                                         + (self.AMINOACIDS_DICT[aminoacid]
                            ['PKA_AMINOACIDS'][1]))
                    count_groups += 1
        return isoelectric_point / count_groups


def run_dna_rna_tools(*parameters: str) -> (list[NucleicAcidSequence]
                                            | NucleicAcidSequence):
    """
        Run DNA and RNA sequence manipulation tools.

        This function accepts a variable number of parameters, with the last parameter
        specifying the name of the tool to be used. The preceding parameters should
        contain one or more DNA or RNA sequences as strings, represented by a combination
        of characters from the set ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u'].

        Parameters:
        *parameters (str): Variable number of DNA or RNA sequences and the tool name.

        Returns:
        [list[str], str]: Depending on the tool used, it returns a list of
        sequences or a single sequence as a string. If only one sequence is processed,
        it returns the sequence as a string. If the tool name is invalid or the answer
        is None, it raises a ValueError.

        Raises:
        - ValueError: If the parameters are None, there are not enough parameters,
          a given parameter is None, the sequences contain characters other than
          ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u'], or the tool name is unknown.
        - ValueError: If an RNA sequence contains 'T' or a DNA sequence contains 'U'.
        - ValueError: If the answer is None or does not contain valid sequences.

        Example:
        run_dna_rna_tools("ATGC", "AUG", "transcribe")
        ['AUGC', 'AUG', 'transcribe']

        run_dna_rna_tools("ATGC", "UAGC", "transcribe")
        Traceback (most recent call last):
          ...
        ValueError: RNA sequence cannot contain T
        """
    available_rna_dna_symbols = \
        ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u']

    if parameters is None:
        raise ValueError('Parameters are None!')
    if len(parameters) < 2:
        raise ValueError('Parameters are not enough!')

    tool_name = parameters[-1]
    sequences = parameters[:-1]

    for sequence in sequences:
        if sequence is None:
            raise ValueError('Given parameter was None')
        else:
            for nucl in sequence:
                if not (available_rna_dna_symbols
                        .__contains__(nucl)):
                    raise ValueError('Parameters are not nucleotide sequences')
    sequences_objects = []
    is_dna = False
    is_rna = False
    for sequence in sequences:
        for nucl in sequence:
            if nucl.upper() == 'T':
                if not is_rna:
                    is_dna = True
                    sequences_objects.append(DNASequence(sequence))
                    break
                else:
                    raise ValueError('RNA sequence cannot contain T')
            if nucl.upper() == 'U':
                if not is_dna:
                    is_rna = True
                    sequences_objects.append(RNASequence(sequence))
                    break
                else:
                    raise ValueError('DNA sequence cannot contain U')
    answer = []

    for nucleic_acid in sequences_objects:
        if tool_name == 'transcribe':
            answer.append(nucleic_acid.transcribe())
        elif tool_name == 'reverse':
            answer.append(nucleic_acid.reverse())
        elif tool_name == 'complement':
            answer.append(nucleic_acid.complement())
        elif tool_name == 'reverse_complement':
            complemented_seq_obj = None
            if is_dna:
                complemented_seq_obj = DNASequence(nucleic_acid.complement())
            else:
                complemented_seq_obj = RNASequence(nucleic_acid.complement())

            answer.append(complemented_seq_obj.reverse())
        else:
            raise ValueError('Unknown tool')

    if len(answer) < 2:
        return answer[0]
    elif answer is None:
        raise ValueError('Answer is None')
    else:
        return answer