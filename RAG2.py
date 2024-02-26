from abc import ABC, abstractmethod

class BiologicalSequence(ABC):

    def __init__(self, sequence: str):
        self.sequence = sequence
 
    def __len__(self) -> int:
        return len(self._sequence)

    @abstractmethod
    def is_valid(self) -> bool:
        """
        Checks if the sequence matches the specified alphabet.
        Returns:
            bool: True if the sequence is correct, otherwise False.
        """
        pass

    def __getitem__(self, key):
        return self.sequence[key]
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.sequence}')"



class NucleicAcidSequence(BiologicalSequence):
    complement_map = {}

    def is_valid(self) -> bool:
        return all(nucleotide in self.complement_map for nucleotide in self.sequence)

    def complement(self):
        return ''.join(self.complement_map[nucleotide] for nucleotide in self.sequence)

    def gc_content(self):
        gc_content = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) if self.sequence else 0
        return gc_content

class DNASequence(NucleicAcidSequence):
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))

class RNASequence(NucleicAcidSequence):
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):

    def is_valid(self) -> bool:
        amino_acids = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy"
        return all(aa in amino_acids for aa in self.sequence)
    

    def one_to_three_letter_code(self) -> str:
        """
        This function converts a protein sequence from one-letter amino acid code to three-letter code.
    
        Args:
            sequence (str): The input protein sequence in one-letter code.
        
        Returns:
            str: The converted protein sequence in three-letter code.
        """
        AMINO_ACIDS = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
               'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
               'T': 'Thr', 'V': 'Val',
               'W': 'Trp', 'Y': 'Tyr'}
        three_letter_code = [AMINO_ACIDS.get(aa.upper()) for aa in self.sequence]
        return '-'.join(three_letter_code)