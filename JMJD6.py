from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    @abstractmethod
    def seq_len():
        pass

    def seq_slices():
        pass

    def seq_indexes():
        pass

    def seq_to_string():
        pass

    def check_alphabet():
        pass


class NucleicAcidSequence(BiologicalSequence):
    COMPLEMENTARITY_DICT = {"A": "T", "T": "A", "C": "G", "G": "C"}
    alphabet = ["A", "T", "U", "G", "C"]

    def __init__(self, seq):
        self.seq = seq

    def seq_len(self):
        return len(self.seq)

    def seq_slices(self, start, stop):
        slice = self.seq[start:stop]
        return slice

    def seq_indexes(self, index):
        return self.seq[index]

    def seq_to_string(self):
        return str(self.seq)

    def check_alphabet(self):
        if set(self.seq) <= set(self.alphabet):
            return 'That`s a nucleic acid sequence'
        else:
            return 'That`s not a nucleic acid sequence'

    def complement(self):
        complement = ""
        for nucl in self.seq:
            complement += self.COMPLEMENTARITY_DICT.get(nucl)
        return complement

    def gc_content(self):
        return (self.seq.count('C') + self.seq.count('G'))/len(self.seq)


class DNASequence(NucleicAcidSequence):
    COMPLEMENTARITY_DICT = {"A": "T", "T": "A", "C": "G", "G": "C"}
    TRANSCRIBE_DICT = {"A": "U", "T": "A", "C": "G", "G": "C"}
    alphabet = ["A", "T", "G", "C"]

    def __init__(self, seq):
        self.seq = seq

    def transcribe(self):
        transcribed = ""
        for nucl in self.seq:
            transcribed += self.TRANSCRIBE_DICT.get(nucl)
        return transcribed


class RNASequence(NucleicAcidSequence):
    COMPLEMENTARITY_DICT = {"A": "U", "U": "A", "C": "G", "G": "C"}
    alphabet = ["A", "U", "G", "C"]

    def __init__(self, seq):
        self.seq = seq


class AminoAcidSequence(BiologicalSequence):
    alphabet = frozenset('ARNDCEQGHILKMFPSTWYV')

    AA_CHARGES = {"A": 0, "R": 1, "N": 0, "D": -1, "C": 0,
                  "Q": 0, "E": -1, "G": 0, "H": 1, "I": 0,
                  "L": 0, "K": 1, "M": 0, "F": 0, "P": 0,
                  "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0}

    def __init__(self, seq):
        self.seq = seq

    def seq_len(self):
        return len(self.seq)

    def seq_slices(self, start, stop):
        slice = self.seq[start:stop]
        return slice

    def seq_indexes(self, index):
        return self.seq[index]

    def seq_to_string(self):
        return str(self.seq)

    def check_alphabet(self):
        if set(self.seq) <= self.alphabet:
            return 'That`s a protein sequence'
        else:
            return 'That`s not a protein sequence'

    def aa_chain_charge(self):
        aa_charge = 0
        for amino in self.seq:
            aa_charge += self.AA_CHARGES.get(amino)
        return aa_charge
