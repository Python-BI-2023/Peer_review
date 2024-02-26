from abc import ABC, abstractmethod


class BiologicalSequence(ABC):
    @abstractmethod
    def __init__(self, sequence):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, slc):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def is_valid_alphabet(self):
        pass


class NucleicAcid(BiologicalSequence):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return self.sequence

    def is_valid_alphabet(self):
        alphabet = type(self).ALPHABET
        if set(self.sequence).issubset(alphabet):
            return True
        else:
            return False

    def complement(self):
        if type(self) == NucleicAcid:
            raise NotImplementedError("Cannot complement NucleicAcid instance")

        map_dict = type(self).MAP
        comp_seq = "".join([map_dict[base] for base in self.sequence])

        return type(self)(comp_seq)

    def gc_content(self):
        if type(self) == NucleicAcid:
            raise NotImplementedError("Cannot gc_content NucleicAcid instance")
        gc_count = sum([1 for base in self.sequence if base in ["C", "G"]])
        gc_content = (gc_count / len(self)) * 100

        return gc_content


class DNASequence(NucleicAcid):
    ALPHABET = set("ATGC")
    MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def transcribe(self):
        transcribed = self.sequence.replace("T", "U")

        return RNASequence(transcribed)


class RNASequence(NucleicAcid):
    ALPHABET = set("AUGC")
    MAP = {"A": "U", "U": "A", "C": "G", "G": "C"}


class AminoAcidSequence(BiologicalSequence):
    ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return self.sequence

    def is_valid_alphabet(self):
        alphabet = type(self).ALPHABET
        if set(self.sequence).issubset(alphabet):
            return True
        else:
            return False

    amino_acid_frequency = {}

    def calculate_aa_freq(self):
        """
        Calculates the frequency of each amino acid in a protein sequence or sequences.

        :param sequences: protein sequence or sequences
        :type sequences: str or list of str
        :return: dictionary with the frequency of each amino acid
        :rtype: dict
        """

        # Creating a dictionary with aminoacid frequencies:
        amino_acid_frequency = {}

        for amino_acid in self.sequence:
            # If the aminoacid has been already in:
            if amino_acid in amino_acid_frequency:
                amino_acid_frequency[amino_acid] += 1
            # If the aminoacid hasn't been already in:
            else:
                amino_acid_frequency[amino_acid] = 1

        return amino_acid_frequency
