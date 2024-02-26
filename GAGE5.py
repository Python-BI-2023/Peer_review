from Bio import SeqIO
from Bio.SeqUtils import GC
from abc import ABC, abstractmethod

def filter_fastq(input_path, gc_bounds=(0, 100), length_bounds=(0, float('inf')), quality_threshold=0, output_filename=None):
    """
    Filters a FASTQ file based on GC content, sequence length, and quality threshold using Biopython.

    Args:
    - input_path (str): Path to the input FASTQ file.
    - gc_bounds (tuple): Tuple specifying the minimum and maximum GC content for filtering. Default is (0, 100).
    - length_bounds (tuple): Tuple specifying the minimum and maximum sequence length for filtering. Default is (0, infinity).
    - quality_threshold (float): Minimum quality score for filtering. Default is 0.
    - output_filename (str): Name of the output file. If None, the default filename will be used.

    Returns:
    - str: Message indicating the success of the filtering process.
    """
    filtered_seqs = []

    with open(input_path, 'r') as fastq_file:
        for record in SeqIO.parse(fastq_file, 'fastq'):
            gc_content = GC(record.seq)
            seq_length = len(record.seq)
            quality_score = sum(record.letter_annotations["phred_quality"]) / seq_length

            if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= seq_length <= length_bounds[1] and
                quality_score >= quality_threshold):
                filtered_seqs.append(record)

    if output_filename is None:
        output_filename = f"filtered_{input_path}"
    elif not output_filename.endswith('.fastq'):
        output_filename += '.fastq'

    with open(output_filename, 'w') as output_file:
        SeqIO.write(filtered_seqs, output_file, 'fastq')

    return "Filtered data was saved into output file"


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences."""

    def __init__(self, sequence):
        """Initialize a BiologicalSequence object with a given sequence."""
        self.sequence = sequence

    def __len__(self):
        """Return the length of the sequence."""
        return len(self.sequence)

    def __getitem__(self, index):
        """Return the item at the specified index."""
        return self.sequence[index]

    def __str__(self):
        """Return the string representation of the sequence."""
        return self.sequence

    @abstractmethod
    def check_alphabet(self):
        """Check if the sequence contains valid alphabet characters."""
        pass

class NucleicAcidSequence(BiologicalSequence):
    """Abstract base class"""

    DNA_LETTERS = set("ATGCatgc")
    RNA_LETTERS = set("AUGCaugc")
    AMINO_ACID_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self):
        """Check if the sequence contains valid nucleic- or aminoacid alphabet characters."""
        return set(self.sequence).issubset(self.DNA_LETTERS | self.RNA_LETTERS | self.AMINO_ACID_LETTERS)

    def complement(self):
        """Return the complement sequence."""
        comp_map_dna = {"A": "T", "G": "C", "T": "A", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
        return ''.join(comp_map_dna.get(base, base) for base in self.sequence)

    def gc_content(self):
        """Return the GC content of the sequence."""
        gc_count = (self.sequence.upper().count('G') + self.sequence.upper().count('C')) / len(self.sequence) * 100
        return gc_count

class DNASequence(NucleicAcidSequence):
    """Class representing a DNA sequence."""

    TRANSCRIBE_DICT = {
        'T': 'U',
        't': 'u'
    }

    def __init__(self, sequence):
        """Initialize a DNASequence object with a given sequence."""
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Invalid DNA sequence")

    def transcribe(self):
        """Transcribe the DNA sequence into an RNA sequence."""
        transcribed_seq = ''.join(self.TRANSCRIBE_DICT.get(base, base) for base in self.sequence)
        return transcribed_seq

class RNASequence(NucleicAcidSequence):
    """Class representing an RNA sequence."""

    def __init__(self, sequence):
        """Initialize an RNASequence object with a given sequence."""
        super().__init__(sequence)
        if not self.check_alphabet():
            raise ValueError("Invalid RNA sequence")

    def reverse(self):
        """Return the reverse of the RNA sequence."""
        return RNASequence(self.sequence[::-1])

class AminoAcidSequence(BiologicalSequence):
    """Class representing an amino acid sequence."""

    def check_alphabet(self):
        """Check if the sequence contains valid amino acid alphabet characters."""
        if not set(self.sequence).issubset(NucleicAcidSequence.AMINO_ACID_LETTERS):
            raise ValueError("Invalid amino acid sequence")

    def amino_acid_profile(self):
        """Return the profile of the amino acid sequence."""
        self.check_alphabet()

        aa_biochemistry = {
            'hydrophobic': ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'],
            'polar': ['S', 'T', 'C', 'N', 'Q', 'Y'],
            '- charged': ['E', 'D'],
            '+ charged': ['K', 'H', 'R']
        }
        profile = {}

        for group in aa_biochemistry:
            profile[group] = 0.0

        for amino_acid in self.sequence:
            for group_name, group_list in aa_biochemistry.items():
                if amino_acid.upper() in group_list:
                    profile[group_name] += 1

        total_length = len(self.sequence)
        for group, count in profile.items():
            profile[group] = round((count / total_length), 2)
        return profile
