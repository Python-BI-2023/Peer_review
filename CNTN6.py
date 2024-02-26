class BiologicalSequence:
    def __init__(self, sequence):
        self.sequence = sequence
        self.ALPHABET = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, idx):
        return self.sequence[idx]

    def __get_slice__(self, start, end):
        if start < 0 or start > len(self.sequence):
            raise IndexError("start index out of range")

        if end < 0 or end > len(self.sequence):
            raise IndexError("end index out of range")

        if start > end:
            raise ValueError("start index must be less than or equal to end index")

        return self.sequence[start:end]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f'{self.__class__.__name__}({self.sequence!r})'

    def is_valid_bioseq(self):
        for char in self.sequence:
            if char not in self.ALPHABET:
                return False
        return True


class NucleicAcidSequence(BiologicalSequence):
    def complement(self, COMPLEMENT_NUCLEOTIDES=None):
        if COMPLEMENT_NUCLEOTIDES is None:
            COMPLEMENT_NUCLEOTIDES = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        complement_seq = ''
        for nucleotide in self.sequence:
            complement_seq += COMPLEMENT_NUCLEOTIDES.get(nucleotide)
        return complement_seq

    def gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C'))/len(self.sequence)


class RNASequence(NucleicAcidSequence):
    REV_TRANSCRTPTION_COMPLEMENT = {'A': 'T', 'U': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, sequence):
        super().__init__(sequence)
        self.ALPHABET = {'A','U','G','C'}

    def reverse_transribe(self):
        return DnaSequence(self.complement(self.REV_TRANSCRTPTION_COMPLEMENT))


class DnaSequence(NucleicAcidSequence):
    TRANSCRTPTION_COMPLEMENT = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    def __init__(self, sequence):
        super().__init__(sequence)
        self.ALPHABET = {'A', 'T', 'G', 'C'}

    def transcribe(self):
        return RNASequence(self.complement(self.TRANSCRTPTION_COMPLEMENT))


class AminoAcidSequence(BiologicalSequence):
    AMINOACID_ALPHABET_1TO3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                               'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                               'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                               'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}

    def __init__(self, sequence):
        super().__init__(sequence)

    def convert_1_to_str3(self) -> str:
        seq3 = ''
        if len(self.sequence) > 0:
            for aminoacid in self.sequence:
                if aminoacid in self.AMINOACID_ALPHABET_1TO3:
                    seq3 += self.AMINOACID_ALPHABET_1TO3[aminoacid]

        return seq3