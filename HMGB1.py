from Bio import SeqIO
from abc import ABC, abstractclassmethod
from Bio.SeqUtils import GC123


RETRANSLATION_DICT = {
        'F': 'TTC', 'f': 'ttc',
        'L': 'TTA', 'l': 'tta',
        'S': 'TCG', 's': 'tcg',
        'Y': 'TAC', 'y': 'tac',
        'C': 'TGC', 'c': 'tgc',
        'W': 'TGG', 'w': 'tgg',
        'P': 'CCC', 'p': 'ccc',
        'H': 'CAT', 'h': 'cat',
        'Q': 'GAA', 'q': 'gaa',
        'R': 'CGA', 'r': 'cga',
        'I': 'ATT', 'i': 'att',
        'M': 'ATG', 'm': 'atg',
        'T': 'ACC', 't': 'acc',
        'N': 'AAT', 'n': 'aat',
        'K': 'AAA', 'k': 'aaa',
        'V': 'GTT', 'v': 'gtt',
        'A': 'GCA', 'a': 'gca',
        'D': 'GAT', 'd': 'gca',
        'E': 'GAG', 'e': 'gag',
        'G': 'GGG', 'g': 'ggg'
    }


THREE_LETTER_ALPHABET = {
                'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': "ASP", 'V': 'VAL',
                'H': 'HIS', 'G': "GLY", 'Q': "GLN", 'E': 'GLU', 'I': 'ILE',
                'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'P': 'PRO', 'S': 'SER',
                'Y': 'TYR', 'T': 'THR', 'W': 'TRP', 'F': 'PHE', 'C': 'CYS',
                'a': 'ala', 'r': 'arg', 'n': 'asn', 'd': "asp", 'v': 'val',
                'h': 'his', 'g': "gly", 'q': "gln", 'e': 'glu', 'i': 'ile',
                'l': 'leu', 'k': 'lys', 'm': 'met', 'p': 'pro', 's': 'ser',
                'y': 'tyr', 't': 'thr', 'w': 'trp', 'f': 'phe', 'c': 'cys'
    }


class BiologicalSequence(ABC):
    @abstractclassmethod
    def __len__():
        pass

    @abstractclassmethod
    def __getitem__():
        pass

    @abstractclassmethod
    def __str__():
        pass

    @abstractclassmethod
    def __is_correct_alphabet__():
        pass


class InvalidInputError(ValueError):
    pass


class NucleicAcidSequence(str, BiologicalSequence):
    def __str__(self):
        sequence_to_str_list = []
        for nucleotide in self:
            sequence_to_str_list.append(nucleotide)
        return ''.join(sequence_to_str_list)

    def is_correct_alphabet(self):
        for nucleotide in self:
            if nucleotide not in self.alphabet:
                return False
        return True

    def complement(self):
        if not self.is_correct_alphabet():
            raise InvalidInputError('Cannot complement: incorrect input sequence. Only nucleotides (in both cases) are supported!')
        sequence_complement_list = []
        for nucleotide in self:
            sequence_complement_list.append(self.complement_dictionary[nucleotide])
        return type(self)(''.join(sequence_complement_list))

    def gc_content(self):
        gc_count = 0
        for nucleotide in self:
            if nucleotide == 'G' or nucleotide == 'C' or nucleotide == 'g' or nucleotide == 'c':
                gc_count += 1
        return gc_count / len(self)


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        self.alphabet = 'AUGCaugc'
        self.complement_dictionary = {
            'A': 'U',
            'U': 'A',
            'G': 'C',
            'C': 'G',
            'a': 'u',
            'u': 'a',
            'g': 'c',
            'c': 'g'
        }


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence: str):
        self.alphabet = 'ATGCatgc'
        self.complement_dictionary = {
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'a': 't',
            't': 'a',
            'g': 'c',
            'c': 'g'
        }

    def transcribe(self) -> RNASequence:
        transcription_dictionary = {
            'A': 'U',
            'T': 'A',
            'G': 'C',
            'C': 'G',
            'a': 'u',
            't': 'a',
            'g': 'c',
            'c': 'g'
        }
        dna_transcript_list = []
        for nucleotide in self:
            dna_transcript_list.append(transcription_dictionary[nucleotide])
        return RNASequence(''.join(dna_transcript_list))


class AminoAcidSequence(str, BiologicalSequence):
    def __init__(self, sequence: str):
        self.alphabet = 'ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv'
        if not self.is_correct_alphabet():
            raise InvalidInputError('Only amino acid names (one-letter) are supported. Please, check your input!')

    def __str__(self) -> str:
        sequence_to_str_list = []
        for amino_acid in self:
            sequence_to_str_list.append(amino_acid)
        return ''.join(sequence_to_str_list)

    def is_correct_alphabet(self) -> bool:
        for amino_acid in self:
            if amino_acid not in self.alphabet:
                return False
        return True

    def info_amino_acid_percentage(self) -> dict:
        info_dict = {}
        for amino_acid in self:
            if amino_acid not in info_dict:
                info_dict[amino_acid] = 1
            else:
                info_dict[amino_acid] += 1
        info_dict.update((key, round(value / len(self) * 100, 2)) for key, value in info_dict.items())
        info_dict = {key: value for key, value in sorted(info_dict.items(), key=lambda item: item[1], reverse=True)}
        return info_dict


def reverse(seqs: list) -> list:
    """
    Produce a list of reverse sequences
    arguments:
    - seqs (list): a list of sequences
    return
    - reverse_seqs (list): a list of reverse sequences
    """
    reverse_seqs = []
    for seq in seqs:
        reverse_seqs.append(seq[::-1])
    return reverse_seqs


def reverse_complement(seqs: list) -> list:
    """
    Produce a list of reverse complementary sequences
    arguments:
    - seqs (list): a list of sequences
    return
    - (list): a list of reverse complementary sequences
    """
    reverse_complement_seqs = []
    for seq in seqs:
        if 'U' in seq or 'u' in seq:
            reverse_complement_seqs.append(str(RNASequence(seq).complement())[::-1])
        else:
            reverse_complement_seqs.append(str(DNASequence(seq).complement())[::-1])
    return reverse_complement_seqs


def calculate_similarity(sequences: list, precision: int = 3, percentages: bool = False) -> dict:
    """
    Calculate similarity in AMINOACIDS between reference sequence and other sequences
    arguments:
    - sequences (list): reference sequence and other sequences for comparison
    - precision (int): a number of decimals to round the number to
    - percentages (bool): whether percentages are returned instead of fractions
    return:
    - similarities (dict): dictionary with compared sequences as keys and percentages/fractions as their values
    """
    similarities = {}
    for i in range(1, len(sequences)):
        similarity = []
        for j in range(0, len(sequences[i])):
            similarity.append(sequences[0][j] == sequences[i][j])
        if percentages:
            similarities[sequences[i]] = round(sum(similarity) * 100 / len(sequences[i]), precision)
        else:
            similarities[sequences[i]] = round(sum(similarity) / len(sequences[i]), precision)
    return similarities


def find_pattern(sequences: list, pattern: str) -> dict:
    """
    Find all non-overlaping instances of a given pattern in sequences
    arguments:
    - sequences (list): sequences to find the pattern in
    - pattern (str): pattern in question
    return
    - finds(dict): dictionary with sequences as keys and lists of indexes of patterns and the number of patterns as values
    """
    finds = {}
    for j in range(0, len(sequences)):
        find = []
        i = 0
        while i < len(sequences[j]):
            pattern_index = sequences[j].find(pattern, i)
            if pattern_index != -1:
                find.append(pattern_index)
                i = pattern_index + len(pattern)
            else:
                break
        finds[sequences[j]] = [len(find)] + find
    return finds


def convert_to_gene(protein: str) -> str:
    """
    Transforming of an amino acid sequence/protein to DNA sequence
    :param protein: amino acid sequence of protein
    :return: sequence of protein in the DNA sequence form
    """
    return ''.join([RETRANSLATION_DICT[aa] for aa in protein])


def recode_3letter_to_1letter(seqs: list, sep='') -> list:
    """
    Transform into a three-letter amino acids entry.
    arguments:
        - seqs (list): list of sequences for transforming to three-letter entire
        - sep (str): separator between AMINOACIDS, default = ''
    return:
        - three_letter_result (list): transformed sequences with separators
    """
    three_letter_result = []
    for seq in seqs:
        threel_form = ''
        for aa in seq:
            threel_form += THREE_LETTER_ALPHABET[aa] + sep
        if sep:
            threel_form = threel_form[:-1]
        three_letter_result.append(threel_form)
    return three_letter_result


def check_fastq_file(path: str):
    """
    Check whether the path and file provided are valid for further processing
    arguments:
    - path (str): path to .fastq file
    return:
    - no return
    """
    import os
    if '.fastq' in os.path.basename(path):
        if os.path.isfile(path):
            with open(path) as fastq:
                line_count = 0
                total_line_count = 0
                for line in fastq:
                    line_count += 1
                    total_line_count += 1
                    if line_count == 1:
                        if not line.startswith('@'):
                            raise ValueError('Incorrect file content! Suggestion: check the lines near line number ' + str(total_line_count) + ' in the provided file.')
                    if line_count != 1 and line_count != 4 and line.startswith('@'):
                        raise ValueError('Incorrect file content! Suggestion: check the lines near line number ' + str(total_line_count) + ' in the provided file.')
                    if line_count == 4:
                        line_count = 0
                if line_count != 0:
                    raise ValueError('Incorrect file content! Suggestion: check the end of the file')
        else:
            raise ValueError(path + ' is not a correct path to a file!')
    else:
        raise ValueError('Invalid input: path to a .fastq file was expected!')


def run_dna_rna_tools(inputs: tuple) -> list | str:
    """
    Produce a list of either transcripts, reverse sequences, complementary sequences or reverse complementary sequences
    arguments:
    - inputs (tuple): an arbitrary amount of strings where the last one is the name of desired operation, and other strings are sequences
    return
    - complement_seqs (list): a list of complementary sequences
    """
    if len(inputs) < 2:
        raise ValueError('Invalid input: the function requires at least one sequence and an operation name!')
    *seqs, operation = inputs
    if operation == 'transcribe':
        transcribed_seqs = []
        for seq in seqs:
            transcribed_seqs.append(str(DNASequence(seq).transcribe()))
        return transcribed_seqs
    elif operation == 'reverse':
        result = reverse(seqs)
    elif operation == 'complement':
        complement_seqs = []
        for seq in seqs:
            if 'U' in seq or 'u' in seq:
                complement_seqs.append(str(RNASequence(seq).complement()))
            else:
                complement_seqs.append(str(DNASequence(seq).complement()))
        return complement_seqs
    elif operation == 'reverse_complement':
        result = reverse_complement(seqs)
    else:
        raise ValueError('Invalid input: unknown operation! Check the last argument.')
    if len(result) == 1:
        result = ''.join(result)
    return result


def run_protein_tools(inputs: tuple, options: str = None) -> list | dict:
    """
    Produce a list or dictionary according to the option specified
    arguments:
    - inputs (tuple): a tuple of inputs
    - options (str): option name
    """
    operations = {
        'similarity': calculate_similarity,
        'pattern': find_pattern,
        '3letter_name': recode_3letter_to_1letter,
        'dna_code': convert_to_gene
    }

    if options == 'similarity':
        result = operations[options](inputs[:-2], inputs[-2], inputs[-1])
        return result
    elif options == 'pattern':
        result = operations[options](inputs[1:len(inputs)], inputs[0])
        return result
    elif options == '3letter_name':
        result = operations[options](inputs[:-1], inputs[-1])
        return result
    elif options == 'dna_code':
        result = []
        for inpt in inputs:
            res = operations[options](inpt)
            result.append(res)
        return result
    elif options == 'length':
        lengths = []
        for inpt in inputs:
            lengths.append(len(AminoAcidSequence(inpt)))
        return lengths
    elif options == 'percentage':
        percentages = []
        for inpt in inputs:
            percentages.append(AminoAcidSequence(inpt).info_amino_acid_percentage())
        return percentages
    else:
        raise ValueError('Incorrect options input, please try again')


def create_results_dir_if_doesnt_exist():
    import os
    if not os.path.exists('fastq_filtrator_results'):
        os.makedirs('fastq_filtrator_results')


def run_beginner_bioinf_tools(*input_data: str, toolbox: str = None, **kwargs: str | tuple | int) -> str | list | dict:
    """
    Performs various operations on nucleic acid, protein and fastq sequences
    arguments:
    - input_data (str or dict): data to be processed
    - toolbox (str): determines which of the three toolkits is used
    return:
    - (str or list or dict): processed data
    """
    if toolbox == 'dna_rna':
        return run_dna_rna_tools(input_data)
    elif toolbox == 'proteins':
        return run_protein_tools(input_data, kwargs['options'])
    else:
        raise ValueError('Invalid input: there is no toolbox corresponding to value ' + str(toolbox) + '!')


def filter_fastq(path_to_fastq: str, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0):
    check_fastq_file(path_to_fastq)
    seqs_filtered = []
    if isinstance(gc_bounds, int):
        gc_bounds = tuple([0, gc_bounds])
    if isinstance(length_bounds, int):
        length_bounds = tuple([0, length_bounds])
    for i, record in enumerate(SeqIO.parse(path_to_fastq, "fastq")):
        gc_check = gc_bounds[0] < GC123(record.seq)[0] < gc_bounds[1]
        length_check = length_bounds[0] < len(record.seq) < length_bounds[1]
        quality_per_letter = record.letter_annotations["phred_quality"]
        quality_check = sum(quality_per_letter) / len(quality_per_letter) >= quality_threshold
        if gc_check and length_check and quality_check:
            seqs_filtered.append(record)
    return seqs_filtered
