def transcribe(sequns: str) -> str:
    """
    That function transcribe DNA to RNA
    :param sequns: DNA sequence
    :return: RNA sequence
    """
    return "".join([transcription_dict[i] for i in sequns])


def reverse(sequns: str) -> str:
    """
    That function reverses DNA sequence
    :param sequns:DNA sequense
    :return: reversed DNA sequence
    """
    return sequns[::-1]


def complement(sequns: str) -> str:
    """
    That function create a complementarity chain fot the DNA sequence
    :param sequns: DNA sequense
    :return: complementarity chain to the DNA sequence
    """
    return "".join([complement_dict[i] for i in sequns])


def reverse_complement(sequns: str) -> str:
    """
    That function create reverse complementarity chain for the DNA sequence
    :param sequns: sequense
    :return: reversed complimentarity chain to the DNA sequence
    """
    rev_sequns = sequns[::-1]
    return "".join([complement_dict[i] for i in rev_sequns])


operations = {
    "transcribe": transcribe,
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
}


transcription_dict = {
    "A": "U",
    "a": "u",
    "T": "A",
    "t": "a",
    "G": "C",
    "g": "c",
    "C": "g",
    "c": "g",
}


complement_dict = {
    "A": "T",
    "a": "t",
    "T": "A",
    "t": "a",
    "G": "C",
    "g": "c",
    "C": "g",
    "c": "g",
}


def dna_rna_tools(*args):
    *sequnses, def_name = args
    for sequns in sequnses:
        if RNA_checking(sequns) == True:
            continue
        else:
            raise ValueError("Incorrect options input, please try again")
    result = []
    for sequns in sequnses:
        res = operations[def_name](sequns)
        result.append(res)
    return result


from os import makedirs
import os.path


def GC_cont_check(sequence_for_filtering: str, gc_bounds: list) -> bool:
    """
    Check the GC content of the sequence for filtering
    :param sequence_for_filtering: analyzed sequence
    :param gc_bounds: list with limitation for GC content
    :return: bool argument for filtering in fasta_filtering function
    """
    GC_counter = 0
    for i in sequence_for_filtering:
        if i == "C" or i == "G" or i == "c" or i == "g":
            GC_counter += 1
    GC_calculc = GC_counter / len(sequence_for_filtering) * 100
    if len(gc_bounds) == 1:
        if GC_calculc <= gc_bounds[0]:
            return True
    if len(gc_bounds) == 2:
        if GC_calculc <= gc_bounds[1] and GC_calculc >= gc_bounds[0]:
            return True


def lenght_chech(sequence_for_filtering: str, length_bounds: list) -> bool:
    """
    Check the lenght of the sequence for filtering
    :param sequence_for_filtering: analyzed sequence
    :param length_bounds: list with limitations for lenght
    :return: bool argument for filtering in fasta_filtering function
    """
    if len(length_bounds) == 1:
        if len(sequence_for_filtering) <= length_bounds[0]:
            return True
    if len(length_bounds) == 2:
        if (
            len(sequence_for_filtering) <= length_bounds[1]
            and len(sequence_for_filtering) >= length_bounds[0]
        ):
            return True


def quality_chech(
    quality_of_sequence_for_filterring: str, quality_threshold: int
) -> bool:
    """
    Check the quality of the sequence for filtering
    :param quality_of_sequence_for_filterring: 2nd key for the sequense with quality of each nucleotide reading
    :param quality_threshold: limitation for the quality of nucleotides reading
    :return: bool argument for filtering in fasta_filtering function
    """
    quality_threshold = str(quality_threshold)
    quality_counter = 0
    for i in quality_of_sequence_for_filterring:
        quality_counter += ord(i)
    quality = quality_counter / len(quality_of_sequence_for_filterring)
    if quality >= ord(quality_threshold):
        return True


def seqs_creation(input_path_input: str) -> dict:
    """
    Reading of the input file and creating the dictinary of the seqs
    with structure
     {name' : ('sequence', 'comment' 'quality')
     }
    :param input_path_input:
    :return: dictinary of esquenses
    """

    path_input = str(input_path_input)
    inline_new_dit_fasta = {}
    outline_new_dict_fasta = {}
    py_file = open(path_input)
    lines = py_file.readlines()
    i = 0
    seq = {}
    seqs = {}
    while i != (len(lines)):
        key = lines[i]
        key = key[:-1]
        value_1 = lines[i + 1]
        value_1 = value_1[:-1]
        value_3 = lines[i + 2]
        value_3 = value_3[:-1]
        value_2 = lines[i + 3]
        if value_2[-1] == "\n":
            value_2 = value_2[:-1]
        value = [value_1, value_3, value_2]
        seq[key] = value
        seqs = {**seq}
        i += 4
    return seqs


def output_creating(input_path: str, output_filename: str, outline_new_dict_fasta: str):
    """
    Create the output file and folder for the file
    :param input_path: the name of input file for extracting sequenses
    :param output_filename:  the name of output file for writing sequenses
    :param outline_new_dict_fasta: the name of the folder fot the keeping output filename
    :return: nothing
    """
    if output_filename != None:
        os.makedirs("fastq_filtrator_resuls", exist_ok=True)
        file_for_output_filename = (
            "fastq_filtrator_resuls/" + output_filename + ".fastq"
        )
        with open(file_for_output_filename, mode="w") as f:
            for key, value in outline_new_dict_fasta.items():
                f.write(key + "\n")
                f.write(value[0] + "\n")
                f.write(value[1] + "\n")
                f.write(value[2] + "\n")
    else:
        os.makedirs("fastq_filtrator_resuls", exist_ok=True)
        file_for_output_filename = "fastq_filtrator_resuls/" + input_path
        with open(file_for_output_filename, mode="w") as f:
            for key, value in outline_new_dict_fasta.items():
                f.write(key + "\n")
                f.write(value[0] + "\n")
                f.write(value[1] + "\n")
                f.write(value[2] + "\n")


def fasta_filtering(
    seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0
):
    if type(gc_bounds) != tuple:
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) != tuple:
        length_bounds = (0, length_bounds)
    inline_new_dit_fasta = {}
    for key, value in seqs.items():
        sequence_for_filtering = value[0]
        quality_of_sequence_for_filterring = value[1]
        if (
            GC_cont_check(sequence_for_filtering, gc_bounds)
            and lenght_chech(sequence_for_filtering, length_bounds)
            and quality_chech(quality_of_sequence_for_filterring, quality_threshold)
        ):
            inline_new_dit_fasta[key] = value
    return inline_new_dict_fasta




def length_info(protein: str) -> int:
    """
    Сounting the length of an amino acid sequence/protein in the number of amino acids
    :param protein:  sequence of protein
    :return: number of amino acids in an amino acid sequence/protein
    """
    return len(protein)


def count_percentage_aa(seq: str) -> dict:
    """
    Count percentage of each amino acid in sequence
    arguments:
        - seq (str): sequence for counting
    return:
        - dict: dictionary with counted percentage
    """
    l = count_length(seq)
    result = {}
    for aa in seq:
        if aa not in result:
            result[aa] = 1
        else:
            result[aa] += 1
    result.update((key, round(value / l * 100, 2)) for key, value in result.items())
    res = {
        key: value
        for key, value in sorted(result.items(), key=lambda item: item[1], reverse=True)
    }
    return res


def compare_pattern(sequence: str, pattern: str) -> bool:
    """
    Compare a given pattern to a fragment of sequence of the same length
    arguments:
    - sequence (str): sequence fragment to compare with the pattern
    - pattern (str): pattern for comparison
    return:
    - (bool): whether pattern and fragment match
    """
    for i in range(0, len(sequence)):
        if not sequence[i] == pattern[i]:
            return False
    return True


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
        for i in range(0, len(sequences[j])):
            if compare_pattern(sequences[j][i : i + len(pattern)], pattern):
                find.append(i)
                i += len(pattern)
        finds[sequences[j]] = [len(find)] + find
    return finds


def get_protein_gene(protein):
    """
    Transforming of an amino acid sequence/protein to DNA sequence
    :param protein: amino acid sequence of protein
    :return: sequence of protein in the DNA sequence form
    """
    return "".join([retrnaslation_dict[aa] for i in protein])


def rename_three_letter_name(seqs: list, sep="") -> list:
    """
    Transform into a three-letter amino acids entry.
    arguments:
        - seqs (list): list of sequences for transforming to three-letter entire
        - sep (str): separator between aminoacids, default = ''
    return:
        - list: transformed sequences with separators
    """
    res = []
    for seq in seqs:
        threel_form = ""
        for aa in seq:
            threel_form = threel_form + threel[aa] + sep
        if sep:
            threel_form = threel_form[:-1]
        res.append(threel_form)
    return res


coperations = {
    "length": length_info,
    "percentage": count_percentage_aa,
    "pattern": find_pattern,
    "3Letter_name": rename_three_letter_name,
    "DNA_code": get_protein_gene,
}


retrnaslation_dict = {
    "F": "TTC",
    "f": "ttc",
    "L": "TTA",
    "l": "tta",
    "S": "TCG",
    "s": "tcg",
    "Y": "TAC",
    "y": "tac",
    "C": "TGC",
    "c": "tgc",
    "W": "TGG",
    "w": "tgg",
    "P": "CCC",
    "p": "ccc",
    "H": "CAT",
    "h": "cat",
    "Q": "GAA",
    "q": "gaa",
    "R": "CGA",
    "r": "cga",
    "I": "ATT",
    "i": "att",
    "M": "ATG",
    "m": "atg",
    "T": "ACC",
    "t": "acc",
    "N": "AAT",
    "n": "aat",
    "K": "AAA",
    "k": "aaa",
    "V": "GTT",
    "v": "gtt",
    "A": "GCA",
    "a": "gca",
    "D": "GAT",
    "d": "gca",
    "E": "GAG",
    "e": "gag",
    "G": "GGG",
    "g": "ggg",
}


threel = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "V": "VAL",
    "H": "HIS",
    "G": "GLY",
    "Q": "GLN",
    "E": "GLU",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "P": "PRO",
    "S": "SER",
    "Y": "TYR",
    "T": "THR",
    "W": "TRP",
    "F": "PHE",
    "C": "CYS",
    "a": "ala",
    "r": "arg",
    "n": "asn",
    "d": "asp",
    "v": "val",
    "h": "his",
    "g": "gly",
    "q": "gln",
    "e": "glu",
    "i": "ile",
    "l": "leu",
    "k": "lys",
    "m": "met",
    "p": "pro",
    "s": "ser",
    "y": "tyr",
    "t": "thr",
    "w": "trp",
    "f": "phe",
    "c": "cys",
}


def protein_tool(*proteins, options=None):
    proteins = list(proteins)

    operations = {
        "compare": compare,
        "length": length_info,
        "percentage": count_percentage_aa,
        "pattern": find_pattern,
        "3Letter_name": rename_three_letter_name,
        "DNA_code": get_protein_gene,
    }

    if options == "compare":
        result = operations[options](proteins[:-2], proteins[-2], proteins[-1])
        return result
    elif options == "3Letter_name":
        result = operations[options](proteins[:-1], proteins[-1])
        return result
    elif options == "length" or options == "percentage" or options == "DNA_code":
        result = []
        for protein in proteins:
            res = operations[options](protein)
            result.append(res)
        return result
    else:
        raise ValueError("Incorrect options input, please try again")


from Bio import SeqIO
#from Bio.SeqUtils import GC

class FastQFilter:
    def __init__(self, input_file, output_file, min_length, min_quality, min_gc):
        self.input_file = input_file
        self.output_file = output_file
        self.min_length = min_length
        self.min_quality = min_quality
        self.min_gc = min_gc

    def filter_fastq(self):
        with open(self.output_file, 'w') as output_handle:
            for record in SeqIO.parse(self.input_file, 'fastq'):
                if (
                    len(record.seq) >= self.min_length
                    and min(record.letter_annotations["phred_quality"]) >= self.min_quality
                    and Bio.SeqUtils.GC(record.seq) >= self.min_gc
                ):
                    SeqIO.write(record, output_handle, 'fastq')