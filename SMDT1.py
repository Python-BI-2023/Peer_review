import string
COMP_BASES_DNA: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't',
                        'g': 'c', 'c': 'g', 't': 'a'}
COMP_BASES_RNA: dict = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'a': 'u',
                        'g': 'c', 'c': 'g', 'u': 'a'}
AA_MASSES: dict = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147, 'Z': 147,
                   'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105,
                   'T': 119, 'W': 204, 'Y': 181, 'V': 117}

class BiologicalSequence:
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return str(self.seq)

    def __repr__(self):
        return f"{self.seq}"

    def __getitem__(self, key):
        if isinstance(key, int):
            item = str(self).__getitem__(key)
            return item
        elif isinstance(key, slice):
            # Get the start, stop, and step from the slice
            piece = [self[i] for i in range(*key.indices(len(self)))]
            return ''.join(piece)
    def check_alphabet(self):
        return set(self.seq) <= set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')

class NucleicAcidSequence(BiologicalSequence):

    def complement(self):
        if 'T' in set(self.seq.upper()):
            comp_seq = [COMP_BASES_DNA[nuc] for nuc in list(self.seq)]
        else:
            comp_seq = [COMP_BASES_RNA[nuc] for nuc in list(self.seq)]
        return ''.join(comp_seq)

    def gc_content(self):
        gc_score = (100 * (self.seq.count('G') + self.seq.count('C')) / len(self.seq))
        return gc_score

class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        rna_seq = self.seq.replace('t', 'u').replace('T', 'U')
        return rna_seq

class RNASequence(NucleicAcidSequence):
    pass

class AminoAcidSequence(BiologicalSequence):
    """

    Calculates the mass (Da) of a protein based on its amino acids sequence.
    Takes a string of amino acids, returns the molecular weight in Da.
    Amino acids in the string should be indicated as one-letter symbols.

    """
    def calculate_protein_mass(self) -> float:
        aa_seq = self.seq.upper()
        mass = 0
        for amino_acid in aa_seq:
            mass += AA_MASSES[amino_acid]
        return mass

