import os
from Bio import SeqIO
from Bio.SeqUtils import GC


def fastq_filtrator(input_path: str, output_filename=None,
                    gc_bounds=None, length_bounds=(0, 2 ** 32), quality_threshold=0):
    """
    Filters fastq file by GC content, length, and quality.
    :param input_path: str, path to fastq_file.
    :param output_filename: str, name of output fastq file with filtered data.
    Optional, takes input file name if not mentioned otherwise.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input records.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: fastq file in fastq_filtrator_results folder.
    Raises ValueError("Too strict conditions") if no record passed the conditions.
    """

    input_path = os.path.abspath(input_path)

    if not os.path.exists('fastq_filtrator_results'):
        os.mkdir('fastq_filtrator_results')

    if output_filename is None:
        output_filename = os.path.basename(input_path)

    output_path = os.path.join('fastq_filtrator_results', output_filename)

    passed_filters = []

    for record in SeqIO.parse(input_path, "fastq"):
        seq_len = len(record.seq)
        gc_content = GC(record.seq)
        quality = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])

        # Working with GC-bounds
        if gc_bounds is not None:
            if isinstance(gc_bounds, tuple):
                if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                    continue
            else:
                if not gc_content <= gc_bounds:
                    continue

        # Working with lenght-bounds
        if isinstance(length_bounds, tuple):
            if not (length_bounds[0] <= seq_len <= length_bounds[1]):
                continue
        else:
            if not seq_len <= length_bounds:
                continue

        # Working with quality threshold-bounds
        if quality < quality_threshold:
            continue

        passed_filters.append(record)

    if not passed_filters:
        raise ValueError("Too strict conditions")

    with open(output_path, "w") as output_fastq:
        SeqIO.write(passed_filters, output_fastq, "fastq")


DNA_NUCLEOTIDES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}


RNA_NUCLEOTIDES = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


def is_dna(seq: str) -> bool:
    """
    Checks whether sequence is DNA.
    :param seq: str, nucleotide sequence
    :return: bool
    """
    unique_chars = set(seq)
    nucleotides = set('ATGCatgc')
    return unique_chars.issubset(nucleotides)


def is_rna(seq: str) -> bool:
    """
    Checks whether sequence is RNA.
    :param seq: str, nucleotide sequence
    :return: bool
    """
    unique_chars = set(seq)
    nucleotides = set('AUGCaugc')
    return unique_chars.issubset(nucleotides)


def transcribe(*seqs: str) -> list:
    """
    Transcribes DNA into RNA. Changes T/t to U/u nucleotides.
    :param seqs: str, DNA sequences
    :return: str if one sequence, else returns list of sequences
    Raises ValueError('Not a DNA') if sequence is not a DNA
    """
    transcribed_seqs = []
    for seq in seqs:
        if not is_dna(seq):
            raise ValueError('Not a DNA')
        transcribed_seqs.append(seq.replace('T', 'U').replace('t', 'u'))
    return transcribed_seqs if len(transcribed_seqs) > 1 else transcribed_seqs[0]


def reverse_transcribe(*seqs: str) -> list:
    """
    Transcribes RNA into DNA. Changes U/u to T/t nucleotides.
    :param seqs: str, RNA sequences
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA') if sequence is not a RNA
    """
    reverse_transcribed_seqs = []
    for seq in seqs:
        if not is_rna(seq):
            raise ValueError('Not a RNA')
        reverse_transcribed_seqs.append(seq.replace('U', 'T').replace('u', 't'))
    return reverse_transcribed_seqs if len(reverse_transcribed_seqs) > 1 else reverse_transcribed_seqs[0]


def reverse(*seqs):
    """
    Return reverse sequence to a given DNA or RNA sequence.
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    reverse_seqs = []
    for seq in seqs:
        if not is_dna(seq) and not is_rna(seq):
            raise ValueError('Not a DNA or RNA')
        reverse_seqs.append(''.join(list(seq)[::-1]))
    return reverse_seqs if len(reverse_seqs) > 1 else reverse_seqs[0]


def complement(*seqs: str) -> list:
    """
    Returns complement to the input DNA or RNA sequence
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    complement_seqs = []
    for seq in seqs:
        if is_dna(seq):
            complement_dna = ''
            for nucleotide in seq:
                complement_dna += DNA_NUCLEOTIDES[nucleotide]
            complement_seqs.append(complement_dna)
        elif is_rna(seq):
            complement_rna = ''
            for nucleotide in seq:
                complement_rna += RNA_NUCLEOTIDES[nucleotide]
            complement_seqs.append(complement_rna)
        else:
            raise ValueError('Not a DNA or RNA')
    return complement_seqs if len(complement_seqs) > 1 else complement_seqs[0]


def reverse_complement(*seqs: str) -> list:
    """
    Returns reverse complement to the input DNA or RNA sequence
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    reverse_complement_seqs = []
    for seq in seqs:
        if is_dna(seq) or is_rna(seq):
            reverse_complement_seqs.append(''.join(complement(reverse(seq))))
        else:
            raise ValueError('Not a DNA or RNA')
    return reverse_complement_seqs if len(reverse_complement_seqs) > 1 else reverse_complement_seqs[0]


def compute_melting_temperature(*seqs: str) -> float:
    """
    Computes melting temperature of the input DNA sequence
    :param seqs: str, DNA sequence
    :return: float, melting temperature, rounded to 1 decimal place
    Raises ValueError('Not a DNA') if sequence is not a DNA.
    Raises ValueError('Wrong length') if sequence's length is less than 6 and more than 50 nucleotides.
    """
    melting_temperatures = []
    for seq in seqs:
        if not is_dna(seq):
            raise ValueError('Not a DNA')
        if 6 <= len(seq) < 14:
            melting_temperatures.append(round(float((seq.upper().count('A') + seq.upper().count('T')) * 2
                                                    + (seq.upper().count('G') + seq.upper().count('C')) * 3),
                                              ndigits=1))
        elif 14 <= len(seq) <= 50:
            melting_temperatures.append(round(64.9 + 41 * float((seq.upper().count('G') +
                                                                 seq.upper().count('C') - 16.4)) / len(seq),
                                              ndigits=1))
        else:
            raise ValueError('Wrong length')
    return melting_temperatures if len(melting_temperatures) > 1 else melting_temperatures[0]


def run_dna_rna_tools(*args: str) -> list:
    """
    Runs dna_rna_tools, available functions: transcribe, reverse_transcribe, reverse, reverse_complement,
    compute_melting_temperature.
    :param args: str, DNA or RNA sequences, last argument stand for desired function
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueErrors depending on the function applied.
    """
    list_of_functions = {'reverse': reverse, 'complement': complement, 'transcribe': transcribe,
                         'reverse_transcribe': reverse_transcribe, 'reverse_complement': reverse_complement,
                         'compute_melting_temperature': compute_melting_temperature}
    *seqs, function = args
    return list_of_functions[function](*seqs)


def compute_gc_content(seq: str) -> float:
    """
    Computes GC-content of the input sequence.
    :param seq: str, nucleotide sequence
    :return: float, result of computation
    """
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def compute_nucleotide_quality(seq: str) -> float:
    """
    Computes average nucleotide phred33 quality.
    :param seq: str, quality sequence
    :return: int, computed average quality
    """
    quality = 0
    for nucleotide in list(seq):
        quality += ord(nucleotide) - 33
    return quality / len(seq)


def filter_gc_content(seqs: dict, gc_bounds=None) -> dict:
    """
    Filters fastq dictionary by GC-content.
    :param seqs: dict, fastq dictionary.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input dictionary.
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    if gc_bounds is None:
        return seqs
    seqs_filtered: dict = {}
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    for name, seq in seqs.items():
        if gc_bounds[0] <= compute_gc_content(seq[0]) <= gc_bounds[1]:
            seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def filter_length(seqs: dict, length_bounds=(0, 2 ** 32)) -> dict:
    """
    Filters fastq dictionary by length.
    :param seqs: dict, fastq dictionary.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    seqs_filtered: dict = {}
    for name, seq in seqs.items():
        if isinstance(length_bounds, tuple):
            if length_bounds[0] <= len(seq[0]) <= length_bounds[1]:
                seqs_filtered[name] = seq
        else:
            if len(seq[0]) <= length_bounds:
                seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def filter_quality(seqs: dict, quality_threshold=0) -> dict:
    """
    Filters fastq dictionary by quality.
    :param seqs: dict, fastq dictionary.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    seqs_filtered: dict = {}
    for name, seq in seqs.items():
        if compute_nucleotide_quality(seq[1]) >= quality_threshold:
            seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def read_fastq_file(input_path: str) -> dict:
    """
    Reads a Fastq file and converts it into a dictionary sequentially.
    :param input_path: str, path to the input Fastq file.
    :return: dict, Fastq data in dictionary format.
    """
    fastq_data = {}

    with open(input_path, 'r') as file:
        header = False
        sequence = False
        comment = False
        quality = False

        for line in file:
            line = line.strip()
            if not header:
                header = line.strip('@')
            elif not sequence:
                sequence = line
            elif not comment:
                comment = line
            elif not quality:
                quality = line
                fastq_data[header] = [sequence, quality]
                header = False
                sequence = False
                comment = False
                quality = False

    return fastq_data


def write_fastq(fastq_data: dict, output_filename: str):
    """
    Writes fastq dictionary into the fastq file
    :param fastq_data: dict, fastq data
    :param output_filename: str, name of the output fastq file.
    """
    with open(output_filename, 'w') as fastq_file:
        for key, value in fastq_data.items():
            name = key
            sequence = value[0]
            quality = value[1]
            fastq_file.write(f'@{name}\n')
            fastq_file.write(f'{sequence}\n')
            fastq_file.write('+\n')
            fastq_file.write(f'{quality}\n')


def run_fastq_tools(input_path: str, output_filename=None,
                    gc_bounds=None, length_bounds=(0, 2 ** 32), quality_threshold=0):
    """
    Filters fastq file by GC content, length, and quality.
    :param input_path: str, path to fastq_file.
    :param output_filename: str, name of output fastq file with filtered data.
    Optional, takes input file name if not mentioned otherwise.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input dictionary.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: fastq file in fastq_filtrator_results folder.
    Raises ValueError("Too strict conditions") if the return dictionary in any of the functions is empty.
    """
    input_path = os.path.abspath(input_path)

    if not os.path.exists('fastq_filtrator_results'):
        os.mkdir('fastq_filtrator_results')

    if output_filename is None:
        output_filename = os.path.basename(input_path)

    output_path = os.path.join('fastq_filtrator_results', output_filename + '.fastq')

    fastq_dict = read_fastq_file(input_path=input_path)
    filtered_fastq = filter_quality(filter_gc_content(filter_length(fastq_dict, length_bounds=length_bounds),
                                                      gc_bounds=gc_bounds), quality_threshold=quality_threshold)

    write_fastq(filtered_fastq, output_path)


PROTEIN_ALPHABET = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                    "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}


RNA_ALPHABET = {"A", "U", "G", "C"}


AMINO_ACID_MASSES = {"A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "Q": 128.05858, "E": 129.04259, "G": 57.02146, "H": 137.05891,
    "I": 113.08406, "L": 113.08406, "K": 128.09496, "M": 131.04049, "F": 147.06841,
    "P": 97.05276, "S": 87.03203, "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841}


HYDROPHOBIC_AMINOACIDS = {"A", "V", "L", "I", "P", "F", "W", "M"}


DNA_CODONS = {"A": ["GCT", "GCC", "GCA", "GCG"], "C": ["TGT", "TGC"], "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"], "F": ["TTT", "TTC"], "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"], "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "M": ["ATG"], "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "Q": ["CAA", "CAG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "T": ["ACT", "ACC", "ACA", "ACG"],
    "V": ["GTT", "GTC", "GTA", "GTG"], "W": ["TGG"], "Y": ["TAT", "TAC"], "*": ["UAA", "UAG", "UGA"]}

RNA_CODONS = {"F": ["UUC", "UUU"], "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "I": ["AUU", "AUC", "AUA"], "M": ["AUG"], "V": ["GUU", "GUC", "GUA", "GUG"],
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"], "P": ["CCU", "CCC", "CCA", "CCG"],
    "T": ["ACU", "ACC", "ACA", "ACG"], "A": ["GCU", "GCC", "GCA", "GCG"], "Y": ["UAC", "UAU"],
    "*": ["UAA", "UAG", "UGA"], "H": ["CAU", "CAC"], "Q": ["CAA", "CAG"], "N": ["AAU", "AAC"],
    "K": ["AAA", "AAG"], "D": ["GAU", "GAC"], "E": ["GAA", "GAG"], "C": ["UGU", "UGC"],
    "W": ["UGG"], "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"], "G": ["GGU", "GGC", "GGA", "GGG"]}


def is_protein(seq: str) -> bool:
    """
    Checks if input sequence consists of aminoacids, returns boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= PROTEIN_ALPHABET


def is_rna(seq: str) -> bool:
    """
    Checks if input sequence is RNA, returns boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= RNA_ALPHABET


def compute_molecular_weight(protein: str) -> float:
    """
    Compute molecular weight (g/mol) of protein sequence.

    Argument:
    - proteins (str): protein sequence.

    Return:
    - float, computed molecular weight (float rounded to 3 decimal places).
    """
    molecular_weight = 0
    for amino_acid in protein.upper():
        molecular_weight += AMINO_ACID_MASSES[amino_acid]
    return round(molecular_weight, 3)


def compute_length(protein: str) -> int:
    """
    Compute the length of the input protein sequence.

     Argument:
    - proteins (str): protein sequence.

    Return:
    - int, computed length
    """
    return len(protein)


def protein_to_dna(protein: str) -> str:
    """
    Returns possible variants of DNAs for a given protein sequence.

    Argument:
    - protein (str): protein sequence.

    Return:
    - string, variants of nucleic acids.
    If several codons correspond to a given amino acid they are displayed with a '/'.

    Does not distinguish between lowercase and uppercase letters.

    Examples:

    -'MACDRS' -> 'ATG GCT/GCC/GCA/GCG TGT/TGC GAT/GAC CGT/CGC/CGA/CGG/AGA/AGG TCT/TCC/TCA/TCG/AGT/AGC'
    -'MaCdrS' -> 'ATG GCT/GCC/GCA/GCG TGT/TGC GAT/GAC CGT/CGC/CGA/CGG/AGA/AGG TCT/TCC/TCA/TCG/AGT/AGC'
    """
    nucleic_acid_seq = []
    for aa in protein.upper():
        codons = "/".join(DNA_CODONS[aa])
        nucleic_acid_seq.append(codons)
    return " ".join(nucleic_acid_seq)


def get_frequencies_of_amino_acids(protein: str) -> dict:
    """
    Calculates the number of each aminoacid in a given protein sequence.

    Argument:
    - protein (str): protein sequence.

    Return:
    - dictionary, where a key is the aminoacid letter and value is number of this aminoacid.

    Does not distinguish between lowercase and uppercase letters.

    Examples:

    -'MACDRS' -> {'M': 1, 'A': 1, 'C': 1, 'D': 1, 'R': 1, 'S': 1}
    -'MaCdrS' -> {'M': 1, 'A': 1, 'C': 1, 'D': 1, 'R': 1, 'S': 1}
    """
    amino_acids_dict = {}
    for aa in protein.upper():
        if aa in amino_acids_dict:
            amino_acids_dict[aa] += 1
        else:
            amino_acids_dict[aa] = 1
    return amino_acids_dict


def compute_hydrophobicity(protein: str) -> float:
    """
    Compute the percentage of hydrophobic aminoacids in protein sequence.

    Argument:
    - protein (str): protein sequence. Includes hydrophobic
    and hydrophilic aminoacids.

    Return:
    - tuple with protein sequence and computed percentage
    of hydrophobic aminoacids.
    """
    count_of_hydrophobic = 0
    for aa in protein:
        if aa in HYDROPHOBIC_AMINOACIDS:
            count_of_hydrophobic += 1
    return round(count_of_hydrophobic / len(protein) * 100, 3)


def translate_rna(rna: str) -> str:
    """
    Perform the translation of mRNA seguence into protein sequence.

    Argument:
    - rna (str): mRNA sequence. Must contain start-codon and one of
    the stop-codons.

    Return:
    - str, protein sequence after translation.
    Always starts with "M" and ends with "*".
    """
    triplets = [rna[i:i + 3].upper() for i in range(0, len(rna), 3)]
    protein = []
    for triplet in triplets:
        for aminoacid in RNA_CODONS.keys():
            if triplet in RNA_CODONS[aminoacid]:
                protein.append(aminoacid)

    if protein[-1] != "*":
        raise ValueError("Stop-codon (*) is absent in mRNA")
    if protein[0] != "M":
        raise ValueError("Start-codon (M) is absent in mRNA")

    start = protein.index("M")
    stop = protein.index("*")
    return "".join(protein[start:stop + 1])


def check_mutations(rna: str, protein: str) -> str:
    """
    Check missense mutations in the protein sequence after translation.

    Uses additional function "translate_rna(seq)".

    Arguments:
    - rna (str): sequence of mRNA with/without mutations.
    Must contain start-codon and one of the stop-codons.
    - protein (str): protein sequence translated from mRNA.
    Must start with "M" and ends with "*" (stop-codon).

    Note: is_protein(seq) doesn't see "*", but it's used in the other part of function.

    Return:
    - str, if mRNA without mutations return "Protein without mutations."
    If there are mutations in protein, returns aminoacid(s) and their position(s)

    Examples:
    - "AUGGUAGGGAAAUUUUGA", "MVGKF*" ->  "Protein without mutations."
    - "AUGGUAGGGAAAUUUUGA", "MGGVF*" -> "Mutations:G2, V4."
    - "AUGGUAGGGAAAUUUUGA", "MGGKF" –> "ValueError: Stop (*) is absent"
    - "AUGGUAGGGAAAUUUUGA", "GGKF*" –> "ValueError: Start (M) is absent"
    - "AUGAAAAAAUGA", "MK*" -> "ValueError: Different length of translated protein and protein. May be splicing."
    """
    correct_protein = translate_rna(rna)
    mutations = []

    if is_protein(protein[:-1]) is not True:
        raise ValueError("Invalid protein sequence")
    if is_rna(rna) is not True:
        raise ValueError("Invalid RNA sequence")
    if protein[-1] != "*":
        raise ValueError("Stop (*) is absent")
    if protein[0] != "M":
        raise ValueError("Start (M) is absent")
    if len(protein) != len(rna) / 3:
        raise ValueError("Different length of translated protein and protein. May be splicing.")

    for i in range(len(correct_protein)):
        if correct_protein[i] != protein[i]:
            mutations.append(f"{protein[i]}{i + 1}")

    if len(mutations) == 0:
        return "Protein without mutations."
    else:
        return "Mutations: " + ", ".join(mutations) + "."


def run_protein_tools(*args: str):
    """
    Function containing methods for protein analysis.

    Takes arbitrary number of arguments with protein sequences
    and the name of the procedure to be performed (always the last
    argument). Returns the result of the procedure as string, tuple
    or dictionary if one sequence is submitted or list if several.

    Note: if procedure 'check_mutations' is used then input must
    contain only three arguments: RNA sequence, protein sequence
    and the name of procedure itself.
    """
    *seqs, procedure = args
    results = []
    functions = {"compute_molecular_weight": compute_molecular_weight, "compute_length": compute_length,
        "compute_hydrophobicity": compute_hydrophobicity, "get_frequencies_of_amino_acids": get_frequencies_of_amino_acids,
        "protein_to_dna": protein_to_dna}
    if procedure == "check_mutations":
        results.append(check_mutations(seqs[0], seqs[1]))
    else:
        for seq in seqs:
            if is_protein(seq) is not True:
                raise ValueError("Invalid protein sequence")
            if procedure not in functions:
                raise ValueError("Wrong procedure name")
            else:
                results.append(functions[procedure](seq))
    if len(results) == 1:
        return results[0]
    else:
        return results

