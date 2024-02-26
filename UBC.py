import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

SHORT_CODE = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O',
                'a', 'r', 'n', 'd', 'c', 'e', 'q', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'u', 'o']
LONG_CODE = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu',
            'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', 'U': 'Sec', 'O': 'Pyl',
            'a': 'Ala', 'r': 'Arg', 'n': 'Asn', 'd': 'Asp', 'c': 'Cys', 'e': 'Glu', 'q': 'Gln', 'g': 'Gly', 'h': 'His', 'i': 'Ile', 'l': 'Leu',
            'k': 'Lys', 'm': 'Met', 'f': 'Phe', 'p': 'Pro', 's': 'Ser', 't': 'Thr', 'w': 'Trp', 'y': 'Tyr', 'v': 'Val', 'u': 'Sec', 'o': 'Pyl'}
MASS = {'A': 71.08, 'R': 156.2, 'N': 114.1, 'D': 115.1, 'C': 103.1, 'E': 129.1, 'Q': 128.1, 'G': 57.05, 'H': 137.1, 'I': 113.2, 'L': 113.2,
        'K': 128.2, 'M': 131.2, 'F': 147.2, 'P': 97.12, 'S': 87.08, 'T': 101.1, 'W': 186.2, 'Y': 163.2, 'V': 99.13, 'U': 168.05, 'O': 255.3,
        'a': 71.08, 'r': 156.2, 'n': 114.1, 'd': 115.1, 'c': 103.1, 'e': 129.1, 'q': 128.1, 'g': 57.05, 'h': 137.1, 'i': 113.2, 'l': 113.2,
        'k': 128.2, 'm': 131.2, 'f': 147.2, 'p': 97.12, 's': 87.08, 't': 101.1, 'w': 186.2, 'y': 163.2, 'v': 99.13, 'u': 168.05, 'o': 255.3}
DNA_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
RNA_COMPLEMENT = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}


def molecular_weight(seq: str) -> float:
    """
    Function calculates molecular weight of the amino acid chain
        Parameters:
            seq (str): each letter refers to one-letter coded proteinogenic amino acids
    Returns:
        (float) Molecular weight of tge given amino acid chain in Da
    """
    m = 0
    for acid in seq:
        m += MASS[acid]
    return m


def three_letter_code(seq: str) -> str:
    """
    Function converts single letter translations to three letter translations
        Parameters:
            seq (str): each letter refers to one-letter coded proteinogenic amino acids
        Returns:
            (str) translated in three-letter code
    """
    recording = seq.maketrans(LONG_CODE)
    return seq.translate(recording)


def show_length(seq: str) -> int:
    """
    Function counts the number of amino acids in the given sequence
        Parameters:
            seq (str): amino acid sequence
        Returns:
            (int): integer number of amino acid residues
    """
    return len(seq)


def folding(seq: str) -> str:
    """
    Counts the number of amino acids characteristic separately for alpha helixes and beta sheets,
    and gives out what will be the structure of the protein more.
    This function has been tested on proteins such as 2M3X, 6DT4 (PDB ID) and MHC, CRP.
    The obtained results corresponded to reality.
        Parameters:
            seq (str): amino acid sequence
        Returns:
            (str): overcoming structure ('alfa_helix', 'beta_sheet', 'equally')
    """
    alfa_helix = ['A', 'E', 'L', 'M', 'G', 'Y', 'S', 'a', 'e', 'l', 'm', 'g', 'y', 's']
    beta_sheet = ['Y', 'F', 'W', 'T', 'V', 'I', 'y', 'f', 'w', 't', 'v', 'i']
    alfa_helix_counts = 0
    beta_sheet_counts = 0
    for amino_acid in seq:
        if amino_acid in alfa_helix:
            alfa_helix_counts += 1
        elif amino_acid in beta_sheet:
            beta_sheet_counts += 1
    if alfa_helix_counts > beta_sheet_counts:
        return 'alfa_helix'
    elif alfa_helix_counts < beta_sheet_counts:
        return 'beta_sheet'
    elif alfa_helix_counts == beta_sheet_counts:
        return 'equally'


def aminoacid_seqs_only(seqs: list) -> list:
    """
    Leaves only the amino acid sequences from the fed into the function.
        Parameters:
            seqs (list): amino acid sequence list
        Returns:
            aminoacid_seqs (list): amino acid sequence list without non amino acid sequence
    """
    aminoacid_seqs = []
    for seq in seqs:
        unique_chars = set(seq)
        amino_acid = set(SHORT_CODE)
        if unique_chars <= amino_acid:
            aminoacid_seqs.append(seq)
    return aminoacid_seqs


def dna_or_rna_only(seqs):
    """
        The function leaves only DNA and RNA sequences
            Parameters:
                seqs - sequences fed to the input of the main function
            Return:
                list of DNA and RNA sequences only
    """
    dna_seqs = []
    for seq in seqs:
        unique_chars = set(seq)
        nuc = set('ATGCatgcUu')
        if unique_chars <= nuc:
            dna_seqs.append(seq)
    return dna_seqs


def transcribe(seq):
    """
        The function returns the transcribed sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) transcript sequence
    """
    transcript = seq.replace('T', 'U').replace('t', 'u')
    return transcript


def reverse(seq):
    """
        The function produces a reverse sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) reverse sequence
    """
    reverse_seq = seq[::-1]
    return reverse_seq


def complement(seq):
    """
        The function produces a complementary sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) complementary sequence
    """
    complement_seq = []
    if 'U' in set(seq):
        using_dict = RNA_COMPLEMENT
    else:
        using_dict = DNA_COMPLEMENT
    for nucl in seq:
        complement_seq.append(using_dict[nucl])
    return ''.join(complement_seq)


def reverse_complement(seq):
    """
        The function produces a reverse and complementary sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) reverse and complementary sequence
    """
    reverse_compl_seq = complement(reverse(seq))
    return reverse_compl_seq


def filter_fastq(input_path, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0, output_filename=None):
    """
    Performs functions for working with fastq.
        Parameters:
            - input_path - takes as input the path to the FASTQ-file.
                        A DNA or RNA sequences can consist of either uppercase or lowercase letters.
            - gc_bounds - GC composition interval (in percent) for filtering, by default it is (0, 100).
                        If you pass one number as an argument, it is considered to be the upper limit.
                        Examples: gc_bounds = (20, 80) - save only reads with GC content from 20 to 80%,
                        gc_bounds = 44.4 - save reads with GC content less than 44.4%.
            - length_bounds - length interval for filtering, everything is similar to gc_bounds,
                        but by default it is equal to (0, 2**32).
            - quality_threshold - threshold value of average read quality for filtering, default is 0 (phred33 scale).
                        Reads with average quality across all nucleotides below the threshold are discarded.
            - output_filename - name of the file with the filtering result, if not specified, the name of the input file is assigned by default.
                        Important: specify the file extension.
        Example input:
            filter_seqs(seqs = 'example_fastq.fastq', gc_bounds = (20, 80), length_bounds = (0, 89), quality_threshold = 34), output_filename = 'my_out_file'
        Return:
            The function returns a file consisting of only those sequences that satisfy all conditions.
            It puts this file into the "fastq_filtrator_resuls" folder and creates it, if it does not exist.
            All described intervals include both upper and lower boundaries.
            Depending on the function being performed, the following returns will occur:
                If the sequences in the file are RNA, then there will be no filtering by gc composition.
                If you supply a tuple of more than 2 values for the gc_bounds and length_bounds arguments,
                you will receive the errors "Incorrect gc_bounds input" and "Incorrect length_bounds input" respectively.
    """
    if type(gc_bounds) is tuple and len(gc_bounds) == 2:
        min_gc_bound = gc_bounds[0]
        max_gc_bound = gc_bounds[1]
    elif isinstance(gc_bounds, (int, float)):
        min_gc_bound = 0
        max_gc_bound = gc_bounds
    else:
        return print("Incorrect gc_bounds input")

    if type(length_bounds) is tuple and len(length_bounds) == 2:
        min_len_bound = length_bounds[0]
        max_len_bound = length_bounds[1]
    elif type(length_bounds) is tuple and len(length_bounds) == 0:
        min_len_bound = 0
        max_len_bound = 2**32
    elif isinstance(length_bounds, (int, float)):
        min_len_bound = 0
        max_len_bound = length_bounds
    else:
        return print("Incorrect length_bounds input")

    filtered_dict_key = []
    with open(input_path, 'r') as input_file:
        for record in SeqIO.parse(input_file, 'fastq'):
            seq_length = len(record.seq)
            seq_gc_content = gc_fraction(record.seq)
            seq_quality = sum(record.letter_annotations["phred_quality"]) / len(record)

            if (min_len_bound <= seq_length <= max_len_bound) and (min_gc_bound <= seq_gc_content <= max_gc_bound) and (seq_quality >= quality_threshold):
                filtered_dict_key.append(record)

    if output_filename is None:
        output_filename = input_path.split('.')[0] + '_filtered.fastq'
    with open(output_filename, 'w') as output_handle:
        SeqIO.write(filtered_dict_key, output_handle, 'fastq')

    print(f"Filtered FastQ sequences saved to '{output_filename}'")


def amino_acid_tools(*args: str):
    """
    Performs functions for working with protein sequences.
        Parameters:
            The function must accept an unlimited number of protein sequences (str) as input,
            the last  variable must be the function (str) you want to execute.
            The amino acid sequence can consist of both uppercase and lowercase letters.
        Input example:
            amino_acid_tools('PLPKVEL','VDviRIkLQ','PPDFGKT','folding')
        Functions:
            molecular_weight: calculates molecular weight of the amino acid chain
            three_letter_code: converts single letter translations to three letter translations
            show_length: counts the number of amino acids in the given sequence
            folding: counts the number of amino acids characteristic separately for alpha helixes and beta sheets,
                    and gives out what will be the structure of the protein more
            seq_charge: evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7)
        Returns:
            If one sequence is supplied, a string with the result is returned.
            If several are submitted, a list of strings is returned.
            Depending on the function performed, the following returns will occur:
                molecular_weight (int) or (list): amino acid sequence molecular weight number or list of numbers
                three_letter_code (str) or (list): translated sequence from one-letter in three-letter code
                show_length (int) or (list): integer number of amino acid residues
                folding (str) or (list): 'alpha_helix', if there are more alpha helices
                                        'beta_sheet', if there are more beta sheets
                                        'equally', if the probability of alpha spirals and beta sheets are the same
                seq_charge(str) or (list): "positive", "negative" or "neutral"
    """
    *seqs, function = args
    d_of_functions = {'molecular_weight': molecular_weight,
                      'three_letter_code': three_letter_code,
                      'show_length': show_length,
                      'folding': folding}
    answer = []
    aminoacid_seqs = aminoacid_seqs_only(seqs)
    for sequence in aminoacid_seqs:
        if function in d_of_functions.keys():
            answer.append(d_of_functions[function](sequence))
        elif function == 'seq_charge':
            ac_seq = AminoAcidSequence(sequence)
            ac_seq_charge = ac_seq.seq_charge()
            answer.append(ac_seq_charge)
    if len(answer) == 1:
        return answer[0]
    return answer


def run_dna_rna_tools(*args: str):
    """
        Performs functions for processing one or several DNA and/or RNA sequences.
           Parameters:
               The function must accept an unlimited number of sequences (str) as input.
               the last variable should be the function (str) you want to execute.
               The sequence can consist of both uppercase and lowercase letters.
           Example input:
               def run_dna_rna_tools("ATgAAaC", "cUgAuaC", "reverse")
           Functions:
               transcribe — print the transcribed sequence
               reverse — print the reversed sequence
               complement — print the complementary sequence
               reverse_complement — print the reverse complementary sequence
           Return:
               If a single sequence is specified, a string containing the result is returned.
               If multiple strings are sent, a list of strings is returned.
               Depending on the function being performed, the following returns will occur:
                   - if the sequences you pass are not DNA or RNA, then the result of the function will be a list without them
                   - if the sequence is RNA, then the 'transcribe' procedure will produce the unchanged sequence
    """
    *seqs, function = args
    d_of_functions = {'transcribe': transcribe, 
                      'reverse': reverse,
                      'reverse_complement': reverse_complement,
                      'complement': complement,
                     }
    answer = []
    dna_rna_seqs = dna_or_rna_only(seqs)
    for sequence in dna_rna_seqs:
        answer.append(d_of_functions[function](sequence))
    if len(answer) == 1:
        return answer[0]
    return answer


class BiologicalSequence():
    def __init__(self, seq):
        self.seq = seq
    
    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]
        elif isinstance(index, (list, tuple)) and len(index) == 2:
            return self.seq[index[0]:index[1]]

    def __str__(self):
        return str(self.seq)
    
    def check_alphabet(self):
        return set(str(self.seq)) <= set("ACGTUacgtu")


class NucleicAcidSequence(BiologicalSequence):
    def complement(self):
        complement_seq = []
        DNA_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        RNA_COMPLEMENT = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}
        if 'U' in set(self.seq):
            using_dict = RNA_COMPLEMENT
        else:
            using_dict = DNA_COMPLEMENT
        for nucl in self.seq:
            complement_seq.append(using_dict[nucl])
        return ''.join(complement_seq)
    
    def gc_content(self):
        return gc_fraction(self.seq)


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        return RNASequence(self.complement().replace('T', 'U').replace('t', 'u'))


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq):
        super().__init__(seq)
        self.complement_seq = self.complement()

class AminoAcidSequence(BiologicalSequence):
    def seq_charge(self):
        aminoacid_charge = {'R': 1, 'D': -1, 'E': -1, 'K': 1, 'O': 1, 'r': 1, 'd': -1, 'e': -1, 'k': 1, 'o': 1}
        charge = 0
        for aminoacid in self.seq:
            if aminoacid in aminoacid_charge:
                charge += aminoacid_charge[aminoacid]
        if charge > 0:
            return 'positive'
        elif charge < 0:
            return 'negative'
        else:
            return 'neutral'