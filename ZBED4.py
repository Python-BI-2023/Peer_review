import Bio.SeqIO as SeqIO
import Bio.SeqUtils as SeqUt

class BiologicalSequence:
    """
    A class that implements work with biological sequences (DNA, RNA, protein).
    It has standard attributes: 

    Methods:
    + __len__(self) : count length of correct sequence
    + seq_to_print(self) : removes \t, \n and " " in sequence
    + get_subseq(self, start, stop) : make slice of sequence or return one letter
    + is_type(self) : checks the correctness of the correspondence between type and sequence

    NOTE: correct sequence maens that is_type(self) is True
    """
    def __init__(self, seq:str, type:str) -> None:
        """
        seq (str): biological sequence
        type (str): type of biological sequence (DNA, RNA, protein)
        DNA_ALPHABET (set): DNA nucleotides
        RNA_ALPHABET (set): RNA nucleotides
        PROTEIN_ALPHABET (set): proteinogenic aminoacids

        NOTE: DNA_ALPHABET, RNA_ALPHABET, PROTEIN_ALPHABET are PRIVATE variables
        """
        self.seq = seq
        self.type = type
        self.__DNA_ALPHABET = {"A", "T", "G", "C"}
        self.__RNA_ALPHABET = {"A", "U", "G", "C"}
        self.__PROTEIN_ALPHABET = {"A","C","D","E","F","G","H","I","K","L","M",
                    "N","P","Q","R","S","T","V","W","Y"}

    def __len__(self) -> int:
        """Count length of correct sequence."""
        if self.is_type():
            return len(self.seq)
    
    def seq_to_print(self) -> str:
        """Removes \t, \n and " " in sequence."""
        return self.seq.replace('\n', '').replace('\t', '').replace(" ", '').upper()

    def get_subseq(self, start:int, stop=None) -> list:
        """Make slice of sequence or return one letter if stop is None as default.
        """
        if self.is_type():
            if stop is None:
                return self.seq[start]
            return self.seq[start:stop:1]
    
    def is_type(self) -> bool:
        """Checks the correctness of the correspondence between type and sequence."""
        if self.type == 'DNA':
            return set(self.seq.upper()) <= self.__DNA_ALPHABET
        elif self.type == 'RNA':
            return set(self.seq.upper()) <= self.__RNA_ALPHABET
        elif self.type == 'protein':
            return set(self.seq.upper()) <= self.__PROTEIN_ALPHABET
        elif self.type not in ('DNA', 'RNA', 'protein'):
            raise ValueError(f'Unknown type ({self.type}) of biological sequence!')    



class NucleicAcidSequence(BiologicalSequence):
    """
    A child class of "BiologicalSequence" that implements work with nucleic acids.

    Methods:
    + gc_content(self) : count G and C in sequence and return result in %
    + complement(self) : find compliment sequence

    NOTE: correct sequence maens that is_type(self) is True
    """
    def __init__(self, seq:str, type:str) -> None:
        """
        RNA_COMPLEMENT_NUCLS (dict): compliment nucliotides of RNA in upper and lower cases
        DNA_COMPLEMENT_NUCLS (dict): compliment nucliotides of DNA in upper and lower cases
        """
        super().__init__(seq, type)
        self.__RNA_COMPLEMENT_NUCLS = { "A":"U", "U":"A", "a":"u", "u":"a", "G":"C", "C":"G", "c":"g", "g":"c"}
        self.__DNA_COMPLEMENT_NUCLS = {"A":"T", "T":"A", "a": "t", "t":"a", "G":"C", "C":"G", "g":"c", "c":"g"}

    def gc_content(self) -> float:
        """Count G and C in sequence and return result in %."""
        if self.is_type():
            gc_amount = 0
            for nucl in self.seq:
                if nucl == "G" or nucl == "C":
                    gc_amount += 1
            
            cur_gc_level = gc_amount / self.seq.__len__() * 100
            return cur_gc_level
        
        raise ValueError(f"Sequence doesn't match the current type ({self.type}).")
    
    def complement(self) -> str:
        """Find compliment sequence."""
        if self.is_type():
            compl_nucls_map = {}
            if self.type == 'DNA':
                compl_nucls_map = self.__DNA_COMPLEMENT_NUCLS
            elif self.type == 'RNA':
                compl_nucls_map = self.__RNA_COMPLEMENT_NUCLS
                
            trans_seq = ""
            for nucl in self.seq:
                trans_seq += compl_nucls_map[nucl]
            return trans_seq
            
        raise ValueError(f"Sequence doesn't match the current type ({self.type}).")
    

class DNASequence(NucleicAcidSequence):
    """
    Class has functionality for working with DNA and is also
    a child of the class NucleicAcidSequence.

    Methods:
    + transcribe(self): transcribes DNA into RNA
    """
    def __init__(self, seq: str, type: str) -> None:
        super().__init__(seq, type)

    def transcribe(self) -> str:
        """
        Transcribes DNA into RNA
        """
        if self.is_type() and self.type == 'DNA':
            return self.seq.replace("T", "U").replace("t", "u")
        elif self.is_type() and self.type == 'RNA':
            raise ValueError(f"Sequence is RNA.")
        raise ValueError(f"Sequence isn't RNA or DNA.")

class RNASequence(NucleicAcidSequence):
    """
    Class has functionality for working with RNA and is also
    a child of the class NucleicAcidSequence.
    """
    def __init__(self, seq: str, type: str) -> None:
        super().__init__(seq, type)
        self.__RNA_CODONS = {
            "F": ["UUC", "UUU"],
            "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
            "I": ["AUU", "AUC", "AUA"],
            "M": ["AUG"],
            "V": ["GUU", "GUC", "GUA", "GUG"],
            "S": ["UCU", "UCC", "UCA", "UCG"],
            "P": ["CCU", "CCC", "CCA", "CCG"],
            "T": ["ACU", "ACC", "ACA", "ACG"],
            "A": ["GCU", "GCC", "GCA", "GCG"],
            "Y": ["UAC", "UAU"],
            "*": ["UAA", "UAG", "UGA"],
            "H": ["CAU", "CAC"],
            "Q": ["CAA", "CAG"],
            "N": ["AAU", "AAC"],
            "K": ["AAA", "AAG"],
            "D": ["GAU", "GAC"],
            "E": ["GAA", "GAG"],
            "C": ["UGU", "UGC"],
            "W": ["UGG"],
            "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
            "S": ["AGU", "AGC"],
            "G": ["GGU", "GGC", "GGA", "GGG"],
        }

    def translation(self) -> str:
        """
        Perform the translation of mRNA seguence into protein sequence.

        Argument:
        self.seq: mRNA sequence. Must contain start-codon and one of
        the stop-codons.

        Return:
        - str, protein sequence after translation.
        Always starts with "M" and ends with "*".
        """
        if self.is_type() and self.type == 'RNA':
            triplets = [self.seq[i : i + 3].upper() for i in range(0, self.seq.__len__(), 3)]
            protein = []
            for triplet in triplets:
                for aminoacid in self.__RNA_CODONS.keys():
                    if triplet in self.__RNA_CODONS[aminoacid]:
                        protein.append(aminoacid)

            if self.seq.__len__() % 3 != 0:
                raise ValueError("The length of the sequence is a multiple of 3!")
            elif protein[-1] != "*":
                raise ValueError("Stop-codon (*) is absent in mRNA!")
            if protein[0] != "M":
                raise ValueError("Start-codon (M) is absent in mRNA!")

            start = protein.index("M")
            stop = protein.index("*")
            return "".join(protein[start : stop + 1])
        raise ValueError("Input sequence isn't RNA!")    

class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq: str, type: str) -> None:
        super().__init__(seq, type)
        self.__GYDROPHOBIC_AMINOACIDS = {"A", "V", "L", "I", "P", "F", "W", "M"}

    def compute_hydrophobicity(self) -> tuple:
        """
        Compute the percentage of gydrophobic aminoacids in protein sequence.

        Argument:
        - protein (str): protein sequence. Includes hydrophobic
        and hydrophilic aminoacids.

        Return:
        - percentage of gydrophobic aminoacids.
        """
        count_of_gydrophobic = 0
        for aa in self.seq:
            if aa in self.__GYDROPHOBIC_AMINOACIDS:
                count_of_gydrophobic += 1

        percentage = round(count_of_gydrophobic / self.seq.__len__() * 100, 3)

        return percentage


def filter_fastq(input_file, output_file, gc_bounds=(60, 100), length=(0, 2**32), quality_threshold=0):
    """
    Filters Fastq format files by three attributes:
    length, GC-composition, quality. Output filtered sequences to a new file.

    Arguments:
    - input_file (str): name of file with sequences
    - out_put (str): name of the file with filtered sequences
    - gc_bounds (tuple): permissible range of GC
    - length (tuple): permissible range of length
    - quality_threshold (int): minimum permissible value of quality
    """    
    filtered_records = []

    left_gc_bound = gc_bounds[0]
    right_gc_bound = gc_bounds[1]

    left_len_bound = length[0]
    right_len_bound = length[1]

    with open(output_file, "w") as output_name:

        for record in SeqIO.parse(input_file, "fastq"):
            gc_content = SeqUt.GC(record.seq)
            if len(record.seq) >= left_len_bound and len(record.seq) <= right_len_bound:
                if min(record.letter_annotations["phred_quality"]) >= quality_threshold:
                    if gc_content >= left_gc_bound and gc_content <= right_gc_bound:
                        filtered_records.append(record)
                        SeqIO.write(record, output_name, "fastq")

    print(f'Filtered sequences were recorded in {output_file}')    



