from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class NuclAcidnucleotideError(ValueError):
    """Custom error for NucleicAcidSequence classes"""
    pass


class BiologicalSequence(str):
    def __init__(self, sequence):
        self.sequence = sequence

    def check_seq(self):
        valid_nucleotide_symbols = {'A', 'C', 'G', 'T', 'U'}
        valid_prot_symbols = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}
        seq_symbols = set(self.sequence.upper())
        if seq_symbols.issubset(valid_nucleotide_symbols):
            print('NA sequence')
        elif seq_symbols.issubset(valid_prot_symbols):
            print('AA sequence')
        else:
            raise ValueError('Incorrect sequence input!')
        
    def __str__(self):
        return self.sequence
    
    def __repr__(self):
        return f'BiologicalSequence("{self.sequence}")'


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        self.complement_dict = None
        super().__init__(sequence)
    
    def complement(self):
        if self.complement_dict == None:
            raise NotImplementedError('It is a basic NA class. You should implement it for descendant class: DNASequence or RNASequence.')
        result = type(self)(''.join([self.complement_dict[nuc] for nuc in self.sequence]))
        return result
        
    def gc_content(self):
        gc_sum = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        return 100 * gc_sum / len(self.sequence)
        
    def __repr__(self):
        return f'NucleicAcidSequence("{self.sequence}")'


class DNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence)
        self.complement_dict = {
            'a': 't', 'A': 'T',
            't': 'a', 'T': 'A',
            'g': 'c', 'G': 'C',
            'c': 'g', 'C': 'G'
        }
        if 'U' in self.sequence.upper():
            raise NuclAcidnucleotideError('U-contain sequence is not proper DNA sequence')

    def transcribe(self):
        transcription_dict = {
        'a': 'a', 'A': 'A',
        't': 'u', 'T': 'U',
        'g': 'g', 'G': 'G',
        'c': 'c', 'C': 'C'
    }
        result = RNASequence(''.join([transcription_dict[nuc] for nuc in self.sequence]))
        return result

    def __repr__(self):
        return f'DNASequence("{self.sequence}")'


class RNASequence(NucleicAcidSequence):
    def __init__(self, sequence):
        super().__init__(sequence)
        self.complement_dict = {
            'a': 'u', 'A': 'U',
            'u': 'a', 'U': 'A',
            'g': 'c', 'G': 'C',
            'c': 'g', 'C': 'G'
        }
        if 'T' in self.sequence.upper():
            raise NuclAcidnucleotideError('T-contain sequence is not proper RNA sequence')
    
    def __repr__(self):
        return f'RNASequence("{self.sequence}")'


class AminoAcidSequence(BiologicalSequence):
    def gravy(self):
        """Calculate GRAVY (grand average of hydropathy) value"""
        gravy_aa_values = {'L': 3.8,
                           'K': -3.9,
                           'M': 1.9,
                           'F': 2.8,
                           'P': -1.6,
                           'S': -0.8,
                           'T': -0.7,
                           'W': -0.9,
                           'Y': -1.3,
                           'V': 4.2,
                           'A': 1.8,
                           'R': -4.5,
                           'N': -3.5,
                           'D': -3.5,
                           'C': 2.5,
                           'Q': -3.5,
                           'E': -3.5,
                           'G': -0.4,
                           'H': -3.2,
                           'I': 4.5}
        gravy_aa_sum = 0
        for amino_ac in self.sequence.upper():
            gravy_aa_sum += gravy_aa_values[amino_ac]
        return round(gravy_aa_sum / len(self.sequence), 3)
    
    def __repr__(self):
        return f'AminoAcidSequence("{self.sequence}")'



def make_thresholds(threshold: int | float | tuple) -> tuple:
    """Check thresholds input and convert single value to tuple"""
    if isinstance(threshold, int) or isinstance(threshold, float):
        lower = 0
        upper = threshold
    else:
        lower = threshold[0]
        upper = threshold[1]
    return lower, upper


def filter_fastq(input_path: str,
                 gc_thresholds: int | float | tuple = (20, 80),
                 len_thresholds: int | float | tuple = (0, 2 ** 32),
                 quality_threshold: int | float = 0,
                 output_path: str = 'filtered.fastq'):
    """Filters out sequences from fastq file by the specified conditions:
        - GC-content, inside interval include borders, or, if single value, not bigger than specified
        - length, inside interval include borders, or, if single value, not bigger than specified
        - average phred scores, not less than specified
        Default output file name 'filtered.fastq'
    """
    records = list(SeqIO.parse(input_path, 'fastq'))

    filtered_1_gc_idxs = []
    filtered_2_len_idxs = []
    filtered_3_phred_idxs = []
    filtered_results = []

    min_gc, max_gc = make_thresholds(gc_thresholds)
    min_len, max_len = make_thresholds(len_thresholds)

    for count, record in enumerate(records):
        gc_percent = gc_fraction(record.seq) * 100
        if min_gc <= gc_percent <= max_gc:
            filtered_1_gc_idxs.append(count)

    for idx in filtered_1_gc_idxs:
        if min_len <= len(records[idx].seq) <= max_len:
            filtered_2_len_idxs.append(idx)

    for idx in filtered_2_len_idxs:
        phred_values = records[idx].letter_annotations['phred_quality']
        if sum(phred_values) / len(phred_values) >= quality_threshold:
            filtered_3_phred_idxs.append(idx)

    for idx in filtered_3_phred_idxs:
        filtered_results.append(records[idx])

    with open(output_path, 'w') as file:
        SeqIO.write(filtered_results, file, 'fastq')
