ALPHABET_FOR_DNA = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
ALPHABET_FOR_RNA = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
ALPHABET_FOR_PROTEIN = set('FLIMVSPTAYHQNKDECWRG')
COMPLEMENT_ALPHABET_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                           'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
COMPLEMENT_ALPHABET_RNA = {'U': 'A', 'A': 'U', 'G': 'C', 'C': 'G',
                           'u': 'a', 'a': 'u', 'g': 'c', 'c': 'g'}
DICT_RNA_TO_PROTEIN = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L',
    'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I',
    'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S',
    'AGC': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y',
    'UAA': 'stop', 'UAG': 'stop', 'UGA': 'stop', 'CAU': 'H',
    'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'UGG': 'W',
    'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C', 'CGU': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'AUG': 'M'
}

DICT_MOLECULAR_MASS = {
    'G': 75, 'A': 89, 'V': 117, 'L': 131, 'I': 131, 'P': 115,
    'F': 165, 'Y': 181, 'W': 204, 'S': 105, 'T': 119, 'C': 121,
    'M': 149, 'N': 132, 'Q': 146, 'D': 133, 'E': 147, 'K': 146,
    'R': 174, 'H': 155
}


class BiologicalSequence(str):
    """
    Abstract class for working with biological sequences
    Attributes of BiologicalSequence class:
    -   seq (str): biological sequence
    -   seq_type (str): which sequence DNA, RNA or protein the sequence is
    """

    def __init__(self, seq):
        self.seq = seq
        self.seq_type = None

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, idx):
        return self.seq[idx]

    def __getslice__(self, start, end):
        return self.seq[start: end]

    def __repr__(self):
        return f'The sequence is: {self.seq}, type is {self.seq_type}'

    def __str__(self):
        return self.seq

    def check_seq_type(self, check_type):
        """
        Check if sequence is check_type of molecule DNA, RNA or protein
        Input:
        -   check_type (str): type of molecule DNA, RNA, protein
        Return:
        -   bool: True is type match with check_type, False if not
        """

        seq_types = {
            'DNA': ALPHABET_FOR_DNA,
            'RNA': ALPHABET_FOR_RNA,
            'Protein': ALPHABET_FOR_PROTEIN
        }

        if check_type not in seq_types:
            raise ValueError(f'There is not such sequence type as {check_type}!')
        return set(self.seq) <= seq_types[check_type]

    def type_definition(self):
        """
        Set sequence type (DNA, RNA, protein).
        Input:
        -   seq (str): sequence that we need to get type:
        Return:
        -   type (str): type that sequence are
        """
        if self.check_seq_type('DNA'):
            self.seq_type = 'DNA'
        elif self.check_seq_type('RNA'):
            self.seq_type = 'RNA'
        elif self.check_seq_type('Protein'):
            self.seq_type = 'Protein'
        else:
            raise ValueError(f'Sequence {self.seq} can not be analysed!')


class NucleicAcidSequence(BiologicalSequence):
    """
    Methods for nucleic acid sequences: DNA & RNA
    Attributes:
    -   seq (str): nucleic acid sequence
    -   seq_type (str): type of nucleic acid sequence (DNA or RNA)
    -   complement_alphabet (dict): dictionary with complement pairs for DNA or RNA
    """

    def __init__(self, seq):
        super().__init__(seq)
        super().type_definition()
        if super().check_seq_type('DNA'):
            self.complement_alphabet = COMPLEMENT_ALPHABET_DNA
        elif super().check_seq_type('RNA'):

            self.complement_alphabet = COMPLEMENT_ALPHABET_RNA
        else:
            raise ValueError(f'Sequence {self.seq} is not nucleic acid!')

    def complement(self):
        """
        Complement DNA or RNA sequences to RNA or DNA sequence
        Return:
        -   complement sequence
        """
        output_seq = []
        for nucleotide in self.seq:
            output_seq.append(self.complement_alphabet[nucleotide])
        return NucleicAcidSequence(''.join(output_seq))

    def gc_content(self):
        """
        Calculates GC composition for DNA or RNA
        Return:
        -   sequence gc content, %
        """
        gc_result = (self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100
        return round(gc_result, 3)


class DNASequence(NucleicAcidSequence):
    """
    DNA sequence.
    Attributes:
    -   seq (str): nucleic acid sequence
    -   seq_type (str): type of nucleic acid sequence (DNA or RNA)
    -   complement_alphabet (dict): dictionary with complement pairs for DNA or RNA
    """

    def __init__(self, seq):
        super().__init__(seq)
        super().type_definition()
        if not super().check_seq_type('DNA'):
            raise ValueError(f'Sequence {self.seq} is not DNA!')
        else:
            super().type_definition()

    def transcribe(self):
        """
        Transcribe DNA sequence to RNA
        return:
        -   RNA sequence
        """
        output_seq = self.seq.replace('T', 'U').replace('t', 'u')
        return RNASequence(output_seq)


class RNASequence(NucleicAcidSequence):
    """
    RNA sequence.
    Attributes:
    -   seq (str): nucleic acid sequence
    -   seq_type (str): type of nucleic acid sequence (DNA or RNA)
    -   complement_alphabet (dict): dictionary with complement pairs for DNA or RNA
    """

    def __init__(self, seq):
        super().__init__(seq)
        super().type_definition()
        if not super().check_seq_type('RNA'):
            raise ValueError(f'Sequence {self.seq} is not RNA!')

    def transcribe_reverse(self):
        """
        Transcribe and then reverse RNA sequence to DNA
        return:
        -   DNA sequence
        """
        output_seq = self.seq.replace('U', 'T').replace('u', 't')
        return DNASequence(output_seq[::-1])


class AminoAcidSequence(BiologicalSequence):
    """
    Protein (amino acid) sequence.
    Attributes:
    -   seq (str): protein sequence
    -   seq_type (str): The type of molecule in the sequence
    """

    def __init__(self, seq):
        super().__init__(seq)
        super().type_definition()
        if not super().check_seq_type('Protein'):
            raise ValueError(f'Sequence {self.seq} is not protein')

    def counting_molecular_weight(self):
        """
        Counts the molecular mass of a protein sequence seq
        Arguments:
        - seq (str): sequence to count the molecular weight
        Return:
        - output (int): molecular weight value
        """
        output = 0
        for amino_acid in self.seq:
            output += DICT_MOLECULAR_MASS[amino_acid]
        return output - 18 * (len(self.seq) - 1)
