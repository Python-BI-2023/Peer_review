from abc import ABC, abstractmethod


class BiologicalSequence(ABC):

    @abstractmethod
    def __getitem__(self, index) -> str:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq: str):
        self.seq = seq
        self.length = len(self.seq)
        self.check_alphabet()

    def __getitem__(self, slc) -> str:
        if isinstance(slc, int):
            return self.seq.__getitem__(slc - 1)
        elif isinstance(slc, slice):
            new_slice = slice(slc.start - 1, slc.stop, slc.step)
            return self.seq.__getitem__(new_slice)

    def __str__(self) -> str:
        return f'{self.seq}'

    def __repr__(self) -> str:
        return f'{self.seq}'

    def check_alphabet(self) -> bool:
        unique_char_seq = set(self.seq)
        return all(nuc in self.complement_dict for nuc in unique_char_seq)

    def complement(self):
        seq_list = list(self.seq)
        comp_seq = [self.complement_dict[nuc] for nuc in seq_list]
        return self.__class__(''.join(comp_seq))

    def gc_content(self) -> float:
        gc_count = sum(1 for base in self.seq if base.upper() in ['G', 'C'])
        gc_content = (gc_count / len(self.seq)) * 100
        return gc_content


class DNASequence(NucleicAcidSequence):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    def transcribe(self):
        return RNASequence(self.seq.replace('T', 'U').replace('t', 'u'))


class RNASequence(NucleicAcidSequence):
    complement_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


class AminoAcidSequence(BiologicalSequence):
    amino_acids = {'D', 'E', 'R', 'K', 'H', 'N', 'Q', 'S', 'T', 'Y', 'C', 'A', 'G', 'V', 'L', 'I', 'P', 'F',
                   'M', 'W'}
    rna_dict = {'F': 'UUY', 'L': 'YUN', 'I': 'AUH', 'M': 'AUG',
                'V': 'GUN', 'S': 'WSN', 'P': 'CCN', 'T': 'ACN',
                'A': 'GCN', 'Y': 'UAY', 'H': 'CAY', 'Q': 'CAR',
                'N': 'AAY', 'K': 'AAR', 'D': 'GAY', 'E': 'GAR',
                'C': 'UGY', 'R': 'MGN', 'G': 'GGN', 'W': 'UGG'}

    def __init__(self, seq: str):
        self.seq = seq
        self.length = len(self.seq)
        self.check_alphabet()

    def __getitem__(self, slc) -> str:
        if isinstance(slc, int):
            return self.seq.__getitem__(slc - 1)
        elif isinstance(slc, slice):
            new_slice = slice(slc.start - 1, slc.stop, slc.step)
            return self.seq.__getitem__(new_slice)

    def __str__(self) -> str:
        return f'{self.seq}'

    def __repr__(self) -> str:
        return f'{self.seq}'

    def check_alphabet(self) -> bool:
        unique_amino_acids = set(self.seq)
        return unique_amino_acids <= self.amino_acids

    def to_rna(self):
        result = ''.join(self.rna_dict[base] for base in self.seq)
        return RNASequence(result)
