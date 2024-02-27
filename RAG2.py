from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq

class BiologicalSequence(ABC):

    def __init__(self, sequence: str):
        self.sequence = sequence
 
    def __len__(self) -> int:
        return len(self._sequence)

    @abstractmethod
    def is_valid(self) -> bool:
        """
        Checks if the sequence matches the specified alphabet.
        Returns:
            bool: True if the sequence is correct, otherwise False.
        """
        pass

    def __getitem__(self, key):
        return self.sequence[key]
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.sequence}')"



class NucleicAcidSequence(BiologicalSequence):
    complement_map = {}

    def is_valid(self) -> bool:
        return all(nucleotide in self.complement_map for nucleotide in self.sequence)

    def complement(self):
        return ''.join(self.complement_map[nucleotide] for nucleotide in self.sequence)

    def gc_content(self):
        gc_content = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) if self.sequence else 0
        return gc_content

class DNASequence(NucleicAcidSequence):
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))

class RNASequence(NucleicAcidSequence):
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):

    def is_valid(self) -> bool:
        amino_acids = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy"
        return all(aa in amino_acids for aa in self.sequence)
    

    def one_to_three_letter_code(self) -> str:
        """
        This function converts a protein sequence from one-letter amino acid code to three-letter code.
    
        Args:
            sequence (str): The input protein sequence in one-letter code.
        
        Returns:
            str: The converted protein sequence in three-letter code.
        """
        AMINO_ACIDS = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
               'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
               'T': 'Thr', 'V': 'Val',
               'W': 'Trp', 'Y': 'Tyr'}
        three_letter_code = [AMINO_ACIDS.get(aa.upper()) for aa in self.sequence]
        return '-'.join(three_letter_code)


def filter_fastq(input_path: str, output_filename: str = None, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> dict:
    """
        Filters a dictionary of FASTQ sequences based on specified criteria. Saves the output FASTAQ file.

        Args:
            input_path (str): Path to the input FASTQ file.
            output_filename (str, optional): Name of the output FASTQ file. If not provided,
                the filtered sequences will not be saved to a file (default is None).
            gc_bounds (tuple or float, optional): Tuple with lower and upper bounds or a single float
                representing the upper bound for GC content (default is (0, 100)).
            length_bounds (tuple or int, optional): Tuple with lower and upper bounds or a single integer
                representing the upper bound for sequence length (default is (0, 2**32)).
            quality_threshold (int, optional): The threshold for average quality (default is 0).

        Returns:
            dict: Filtered dictionary of sequences. Keys are sequence names, values are tuples of
            (0) sequence and (1) quality scores.
        """
    seqs = read_fastaq(input_path)
    filtered_seqs = {}

    for seq_name, (sequence, quality) in seqs.items():
        if not check_gc_content(sequence, gc_bounds):
            continue

        if not check_length(sequence, length_bounds):
            continue

        if not check_quality(quality, quality_threshold):
            continue

        filtered_seqs[seq_name] = (sequence, quality)

    save_fastq_from_dict(filtered_seqs, output_filename)

    return filtered_seqs


def filter_fastq_with_Bio(input_path: str, output_filename: str = None, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32), quality_threshold: int = 0) -> dict:
    """
    Filters FASTQ sequences from fastq format file based on specified criteria. Saves the output FASTAQ file. Uses Biopython libraries
    """
    filtered_seqs = {}
    
    for record in SeqIO.parse(input_path, "fastq"):
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]
        gc_content =  gc_fraction(record.seq)*100

        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        
        if not (length_bounds[0] <= len(sequence) <= length_bounds[1]):
            continue
        
        if not check_quality(quality_scores, quality_threshold):
            continue
        
        filtered_seqs[record.id] = (sequence, quality_scores)
    
    if output_filename:
        with open(output_filename, "w") as output_handle:
            SeqIO.write((SeqIO.SeqRecord(Seq(seq), id=seq_id, description="", letter_annotations={"phred_quality": quality}) for seq_id, (seq, quality) in filtered_seqs.items()), output_handle, "fastq")
    
    return filtered_seqs


def run_dna_rna_tools(*arguments):
    """
    Executes DNA/RNA sequence manipulation procedures.

    Args:
        *arguments (tuple): Variable-length argument list containing sequences and procedure.

    Returns:
        str or list of str: Result of the selected procedure.
    """
    procedure = arguments[-1]
    sequences = arguments[:-1]
    if not check_valid_sequence(sequences):
        raise ValueError("At least one of your sequences does not correspond to either DNA or RNA")
    if contains_T_and_U_at_the_same_time(sequences):
        raise ValueError(
            "One of your sequences contains both thymine and uracil at the same time, which is not possible((((")
    if procedure == "transcribe":
        return transcribe(sequences)
    elif procedure == "reverse":
        return reverse(sequences)
    elif procedure == "complement":
        return complement(sequences)
    elif procedure == "reverse_complement":
        return reverse_complement(sequences)
    else:
        return "Something went wrong, please, verify the chosen procedure is written correctly"


def run_amino_analyzer(sequence: str, procedure: str, *, weight_type: str = 'average', enzyme: str = 'trypsin'):
    """
    This is the main function to run the amino-analyzer.py tool.
    
    Args:
        sequence (str): The input protein sequence in one-letter code.
        procedure (str): amino-analyzer.py tool has 5 functions at all:
            1. aa_weight - Calculate the amino acids weight in a protein sequence. Return float weight
            weight_type = 'average': default argument for 'aa_weight' function. weight_type = 'monoisotopic' can be
            used as a second option.
            2. count_hydroaffinity - Count the quantity of hydrophobic and hydrophilic amino acids in a protein
            sequence. Return list in order: hydrophobic, hydrophilic
            3. peptide_cutter - This function identifies cleavage sites in a given peptide sequence using a specified
            enzyme. Return list of cleavage sites enzyme = 'trypsin': default argument for 'peptide_cutter' function.
            enzyme = 'chymotrypsin' can be used as a second option.
            4. one_to_three_letter_code - This function converts a protein sequence from one-letter amino acid code
            to three-letter code. Return string of amino acids in three-letter code
            5. sulphur_containing_aa_counter - This function counts sulphur-containing amino acids in a protein
            sequence. Return quantity of sulphur-containing amino acids.

    Returns:
        The result of the specified procedure.

    Raises:
        ValueError: If the procedure is not recognized or if the input sequence contains non-amino acid characters.

    Note: - Supported amino acid characters: V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e,
    q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c. - Make sure to provide a valid procedure name and sequence for
    analysis. :param enzyme: :param sequence: :param procedure: :param weight_type:
    """
    procedures = ['aa_weight', 'count_hydroaffinity', 'peptide_cutter', 'one_to_three_letter_code',
                  'sulphur_containing_aa_counter']
    if procedure not in procedures:
        raise ValueError(f"Incorrect procedure. Acceptable procedures: {', '.join(procedures)}")

    if not is_aa(sequence):
        raise ValueError("Incorrect sequence. Only amino acids are allowed (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, "
                         "T, M, A, G, P, C, v, i, l, e, q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c).")
    result = ''
    if procedure == 'aa_weight':
        result = aa_weight(sequence, weight_type)
    elif procedure == 'count_hydroaffinity':
        result = count_hydroaffinity(sequence)
    elif procedure == 'peptide_cutter':
        result = peptide_cutter(sequence, enzyme)
    elif procedure == 'one_to_three_letter_code':
        result = one_to_three_letter_code(sequence)
    elif procedure == 'sulphur_containing_aa_counter':
        result = sulphur_containing_aa_counter(sequence)
    return result



import os


def check_gc_content(sequence: str, gc_bounds: tuple or float) -> bool: # type: ignore
    """
    Checks if the GC content of a sequence is within the specified bounds.
    Args:
        sequence (str): The input DNA sequence.
        gc_bounds (tuple or float): Tuple with lower and upper bounds or a single float representing the upper bound.
    Returns:
        bool: True if GC content is within bounds, False otherwise.
    """
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    if isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= gc_content <= gc_bounds[1]
    else:
        return gc_content <= gc_bounds


def check_length(sequence: str, length_bounds: tuple or int) -> bool: # type: ignore
    """
    Checks if the length of a sequence is within the specified bounds.

    Args:
        sequence (str): The input DNA sequence.
        length_bounds (tuple or int): Tuple with lower and upper bounds or a single integer representing the upper bound

    Returns:
        bool: True if length is within bounds, False otherwise.
    """
    seq_length = len(sequence)
    if isinstance(length_bounds, tuple):
        return length_bounds[0] <= seq_length <= length_bounds[1]
    else:
        return seq_length <= length_bounds


def check_quality(quality_scores, quality_threshold: int) -> bool:
    """
    Checks the average quality of a sequence, accepting both preprocessed numerical quality scores
    and raw ASCII character quality scores.

    This function allows for flexible handling of quality scores, whether they come directly from FASTQ files
    as ASCII characters or have been preprocessed into numerical scores. It calculates the average quality
    and compares it to a specified threshold to determine if the sequence meets the quality criteria.

    Args:
        quality_scores: Numerical list of quality scores or a string of ASCII quality characters.
        quality_threshold (int): The threshold for average quality.

    Returns:
        bool: True if the average quality is above the threshold, False otherwise.

    Raises:
        ValueError: If `quality_scores` is neither a string nor a list/tuple.
    """
    # If quality_scores is a string, assume these are ASCII characters, raw data from a FASTQ file
    if isinstance(quality_scores, str):
        avg_quality = sum(ord(score) - 33 for score in quality_scores) / len(quality_scores)
    # If quality_scores is a list or tuple, assume these are numerical quality scores
    elif isinstance(quality_scores, (list, tuple)):
        avg_quality = sum(quality_scores) / len(quality_scores)
    else:
        raise ValueError("quality_scores must be either a string or a list/tuple")

    return avg_quality >= quality_threshold




def read_fastaq(input_path: str) -> dict:
    """
    Reads a FASTQ file and returns a dictionary.

    Args:
        input_path (str): The path to the FASTQ file.

    Returns:
        dict: A dictionary where keys are sequence names (starting with "@")
              and values are tuples of (sequence (line 2), quality (line 4)).

    Example:
        Given a FASTQ file like this:

        @Sequence1
        AGCTAGCTAGCTAGCT
        +
        !@#$!@#$!@#$!@#$
        @Sequence2
        CGATCGATCGATCGAT
        +
        !@#$!@#$!@#$!@#$

        The function will return:
        {'@Sequence1': ('AGCTAGCTAGCTAGCT', '!@#$!@#$!@#$!@#$'),
         '@Sequence2': ('CGATCGATCGATCGAT', '!@#$!@#$!@#$!@#$')}

    """
    seqs = {}
    with open(input_path) as fastaq:
        lines = fastaq.readlines()
        number_of_line = 0
        while number_of_line < len(lines):
            if lines[number_of_line].startswith("@"):
                name = lines[number_of_line].strip()
                sequence = lines[number_of_line + 1].strip()
                quality = lines[number_of_line + 3].strip()
                seqs[name] = (sequence, quality)
                number_of_line += 4
            else:
                number_of_line += 1
    return seqs


def save_fastq_from_dict(filtered_seqs: dict, output_filename=None) -> None:
    """
    Save sequences from a dictionary to a FASTQ file.

    Args:
        filtered_seqs (dict): A dictionary where keys are sequence names and values
                              are tuples of (sequence, quality).
        output_filename (str, optional): The output filename. If not provided,
                                         the default is 'fastq_filtrator_results/filtered_data.fastq'.

    Returns:
        None

    """
    if not output_filename:
        output_filename = 'fastq_filtrator_results/filtered_data.fastq'
    else:
        output_filename = f'fastq_filtrator_results/{output_filename}.fastq'

    output_folder = os.path.dirname(output_filename)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(output_filename, 'w') as output_file:
        for name, (sequence, quality) in filtered_seqs.items():
            output_file.write(f'{name}\n{sequence}\n+\n{quality}\n')

def check_valid_sequence(sequences: tuple) -> bool:
    """
    Checks if the input sequences consist only of valid DNA or RNA characters.

    Args:
        sequences (iterable of str): List of sequences to be validated.
    Returns:
        bool: True if all sequences are valid, False otherwise.
    """
    allowed_characters = set('ATGCUatgcu')
    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide not in allowed_characters:
                return False
    return True


def contains_T_and_U_at_the_same_time(sequences: tuple) -> bool:
    """
    Checks if any sequence in the input list contains both 'T' and 'U' nucleotides simultaneously.
    Args:
        sequences (iterable of str): List of sequences to be checked.
    Returns:
        bool: True if any sequence contains both 'T' and 'U', False otherwise.
    """
    for sequence in sequences:
        sequence = sequence.upper()
        if sequence.count("T") and sequences.count("U"):
            return True
    return False


def get_first_sequence(my_tuple: tuple or list[str]) -> str: # type: ignore
    """
    Extracts the first sequence from the input tuple, if applicable.

    Args:
        my_tuple (str or tuple): Input tuple of sequences.

    Returns:
        str or None: The first sequence if available, otherwise None.

    """
    if isinstance(my_tuple, str):
        return my_tuple
    elif len(my_tuple) == 1:
        return str(my_tuple[0])


def is_dna(sequences: tuple) -> bool:
    """
    Checks if all sequences in the input list consist only of valid DNA characters.

    Args:
        sequences (iterable of str): List of sequences to be validated.

    Returns:
        bool: True if all sequences are valid DNA, False otherwise.
    """
    allowed_characters = set('ATGCatgc')
    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide not in allowed_characters:
                return False
    return True


def transcribe(sequences: tuple) -> str or list[str]: # type: ignore
    """
    Transcribes DNA sequences to RNA sequences.

    Args:
        sequences (iterable of str): List of DNA sequences to be transcribed.

    Returns:
        str or list of str: Transcribed RNA sequence(s).
    """
    for sequence in sequences:
        if not is_dna(sequence):
            raise ValueError("At least one of your sequences is RNA instead of DNA, and RNA can not be transcribed")
    first_sequence = get_first_sequence(sequences)
    if first_sequence:
        return first_sequence.replace("T", "U").replace('t', 'u')
    else:
        return [sequence.replace("T", "U").replace('t', 'u') for sequence in sequences]


def reverse(sequences: tuple) -> str or list[str]: # type: ignore
    """
    Reverses the input sequences.

    Args:
        sequences (iterable of str): List of sequences to be reversed.

    Returns:
        str or list of str: Reversed sequence(s).
    """
    first_sequence = get_first_sequence(sequences)
    if first_sequence:
        return first_sequence[::-1]
    else:
        return [sequence[::-1] for sequence in sequences]


def complement(sequences: tuple or list[str]) -> str or list[str]: # type: ignore
    """
    Finds the complement of DNA or RNA sequences.

    Args:
        sequences (str or iterable of str): Input sequence(s).

    Returns:
        str or list of str: Complemented sequence(s).
    """
    if type(sequences) == str:
        sequences = () + (sequences,)
    complement_sequences = []
    for sequence in sequences:
        complement_seq = ''
        if sequence.count("U"):
            nucl_complement_map = {"A": "U", "C": "G", "U": "A", "G": "C", 'a': 'u', 'c': 'g', 'u': 'a', 'g': 'c'}
        else:
            nucl_complement_map = {"A": "T", "C": "G", "T": "A", "G": "C", 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
        for nucleotide in sequence:
            complement_seq += nucl_complement_map[nucleotide]
        complement_sequences.append(complement_seq)
    first_sequence = get_first_sequence(sequences)
    if first_sequence:
        return get_first_sequence(complement_sequences)
    else:
        return complement_sequences


def reverse_complement(sequences: tuple) -> str or list[str]: # type: ignore
    """
    Finds the reverse complement of DNA or RNA sequences.

    Args:
        sequences (str or iterable of str): Input sequence(s).

    Returns:
        str or list of str: Reverse complemented sequence(s).
    """
    return complement(reverse(sequences))

AA_SET = {'V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C',
          'v', 'i', 'l', 'e', 'q', 'd', 'n', 'h', 'w', 'f', 'y', 'r', 'k', 's', 't', 'm', 'a', 'g', 'p', 'c'}
HYDROPHOBIC_AA = ['A', 'V', 'L', 'I', 'P', 'F', 'W', 'M']
HYDROPHILIC_AA = ['R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'K', 'S', 'T', 'Y']
AMINO_ACIDS = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
               'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
               'T': 'Thr', 'V': 'Val',
               'W': 'Trp', 'Y': 'Tyr'}


def is_aa(seq: str) -> bool:
    """
    Check if a sequence contains only amino acids.

        Args:
        seq (str): The input sequence to be checked.

    Returns:
        bool: True if the sequence contains only amino acids, False otherwise.
    """
    unique_chars = set(seq)
    return unique_chars <= AA_SET


def choose_weight(weight: str) -> dict:
    """
    Choose the weight type of amino acids - average or monoisotopic.

    Args:
        weight (str): The type of weight to choose, either 'average' or 'monoisotopic'.

    Returns:
        dict: A dictionary mapping amino acids to their weights based on the chosen type.
    """
    if weight == 'average':
        weights = {
            'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
            'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
            'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
            'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
            }
    elif weight == 'monoisotopic':
        weights = {
            'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
            'E': 129.04259, 'Q': 128.05858, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
            'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
            'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
            }
    else:
        raise ValueError(f"I do not know what '{weight}' is :( \n Read help or just do not write anything except your "
                         f"sequence")

    return weights


def aa_weight(seq: str, weight: str = 'average') -> float:
    """
    Calculate the amino acids weight in a protein sequence.

    Args:
        seq (str): The amino acid sequence to calculate the weight for.
        weight (str, optional): The type of weight to use, either 'average' or 'monoisotopic'. Default is 'average'.

    Returns:
        float: The calculated weight of the amino acid sequence.
    """
    weights_aa = choose_weight(weight)
    final_weight = 0
    for aa in seq.upper():
        final_weight += weights_aa[aa]
    return round(final_weight, 3)


def count_hydroaffinity(seq: str) -> list:
    """
    Count the quantity of hydrophobic and hydrophilic amino acids in a protein sequence.

    Args:
        seq (str): The protein sequence for which to count hydrophobic and hydrophilic amino acids.

    Returns:
        tuple: A tuple containing the count of hydrophobic and hydrophilic amino acids, respectively.
    """
    hydrophobic_count = 0
    hydrophilic_count = 0
    seq = seq.upper()

    for aa in seq:
        if aa in HYDROPHOBIC_AA:
            hydrophobic_count += 1
        elif aa in HYDROPHILIC_AA:
            hydrophilic_count += 1

    return [hydrophobic_count, hydrophilic_count]


def peptide_cutter(sequence: str, enzyme: str = "trypsin") -> str:
    """
    This function identifies cleavage sites in a given peptide sequence using a specified enzyme.
    
    Args: sequence (str): The input peptide sequence. enzyme (str): The enzyme to be used for cleavage. Choose
    between "trypsin" and "chymotrypsin". Default is "trypsin".
        
    Returns: str: A message indicating the number and positions of cleavage sites, or an error message if an invalid
    enzyme is provided.
    """
    cleavage_sites = []
    if enzyme not in ("trypsin", "chymotrypsin"):
        return "You have chosen an enzyme that is not provided. Please choose between trypsin and chymotrypsin."

    if enzyme == "trypsin":  # Trypsin cuts peptide chains mainly at the carboxyl side of the amino acids lysine or
        # arginine.
        for aa in range(len(sequence) - 1):
            if sequence[aa] in ['K', 'R', 'k', 'r'] and sequence[aa + 1] not in ['P', 'p']:
                cleavage_sites.append(aa + 1)

    if enzyme == "chymotrypsin":  # Chymotrypsin preferentially cleaves at Trp, Tyr and Phe in position P1(high
        # specificity)
        for aa in range(len(sequence) - 1):
            if sequence[aa] in ['W', 'Y', 'F', 'w', 'y', 'f'] and sequence[aa + 1] not in ['P', 'p']:
                cleavage_sites.append(aa + 1)

    if cleavage_sites:
        return f"Found {len(cleavage_sites)} {enzyme} cleavage sites at positions {', '.join(map(str, cleavage_sites))}"
    else:
        return f"No {enzyme} cleavage sites were found."


def one_to_three_letter_code(sequence: str) -> str:
    """
    This function converts a protein sequence from one-letter amino acid code to three-letter code.
    
    Args:
        sequence (str): The input protein sequence in one-letter code.
        
    Returns:
        str: The converted protein sequence in three-letter code.
    """
    three_letter_code = [AMINO_ACIDS.get(aa.upper()) for aa in sequence]
    return '-'.join(three_letter_code)


def sulphur_containing_aa_counter(sequence: str) -> str:
    """
    This function counts sulphur-containing amino acids (Cysteine and Methionine) in a protein sequence.
    
    Args:
        sequence (str): The input protein sequence in one-letter code.
        
    Returns:
        str: The number of sulphur-containing amino acids in a protein sequence.
    """
    counter = 0
    for aa in sequence:
        if aa == 'C' or aa == 'M':
            counter += 1
    answer = str(counter)
    return 'The number of sulphur-containing amino acids in the sequence is equal to ' + answer
