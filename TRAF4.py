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

