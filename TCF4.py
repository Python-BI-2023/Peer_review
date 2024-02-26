from Bio import SeqIO
from Bio.SeqUtils import GC

def filter_fastq(input_path: str, quality_threshold: int, output_filename="final_filtered.fastq",  gc_bounds=(40, 60), length_bounds=(50, 350)):
    filename = input_path
    records = SeqIO.parse(filename, "fastq")
    ###quality filter
    good_reads = (rec for rec in records if min(rec.letter_annotations["phred_quality"]) >= quality_threshold)
    result_quality = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
    result_quality_GC = SeqIO.parse("good_quality.fastq", "fastq")
    ###GC content filter
    min_gc_content = gc_bounds[0]
    max_gc_content = gc_bounds[1]
    GC_quality_filt = []
    
    for sequence in result_quality_GC:
        if min_gc_content <= GC(sequence.seq) <= max_gc_content:
            GC_quality_filt.append(sequence)
            
    result_quality = SeqIO.write(GC_quality_filt, "good_quality_GC.fastq", "fastq")
    result_quality_GC_length = SeqIO.parse("good_quality_GC.fastq", "fastq")
    
    ##length filter
    filtered_GC_quality_length = []
    
    for sequence in result_quality_GC_length:
        if len(sequence.seq) >= length_bounds[0] and len(sequence.seq) <= length_bounds[1]:
            filtered_GC_quality_length.append(sequence)
            
    result_quality = SeqIO.write(filtered_GC_quality_length, output_filename, "fastq")
    
    print(result_quality)

#filter_fastq("example_fastq.fastq", 15)


from abc import ABC, abstractmethod

class InvalidInputError(ValueError):
    pass

class BiologicalSequence(ABC, str):
    @abstractmethod
    def __init__(self, seq):
        self.seq = seq
        
    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, index):
        return self.seq[int(index)]
    
    def __repr__(self):
        return __str__(self.seq)
    
    def check_nucleic_acid(self):
        unique_chars = set(self.seq)
        nucleotides_dna = set('ATGCatgc')
        nucleotides_rna = set('AUGCaugc')
        if unique_chars <= nucleotides_dna:
            seq = 'dna'
        elif unique_chars <= nucleotides_rna:
            seq = 'rna'
        else:
            raise InvalidInputError()
            return seq_type
        
class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        super().__init__(seq)
        self.check_nucleic_acid()
        self.length = len(self.seq)
        
    def complement(self):
        list_input = list(self.seq)
        for i in range(len(self.seq)):
            if list_input[i] in self.complement_dict:
                list_input[i] = self.complement_dict[list_input[i]]
        return "".join(list_input)
        
class DNASequence(NucleicAcidSequence):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    def __init__(self, seq):
        super().__init__(seq)
        #self.complement()
    
    def transcribe(self):
        list_input = list(self.seq)
        for i in range(len(self.seq)):
            if (list_input[i] == 'T'):
                list_input[i] = 'U'
            elif (list_input[i] == 't'):
                list_input[i]='u'
        return "".join(list_input)

class RNASequence(NucleicAcidSequence):
    complement_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
    def __init__(self, seq):
        super().__init__(seq)
        #self.complement()
    
class AminoAcidSequence(BiologicalSequence):
    def __init__(self, seq):
        self.seq = seq
        
    def amino_acid_frequency(self):
        """Calculates molecular weight of a protein
    Arguments:
    - seq (str) 1-letter coded protein sequence
    Return:
    - int, molecular weight (g/mol) rounded to integer"""
        freq_dict = {}
        for letter in self.seq:
            if letter in freq_dict:
                freq_dict[letter] += 1
            else:
                freq_dict[letter] = 1
        for letter in freq_dict:
            freq_dict[letter] = round(freq_dict[letter] / len(self.seq) * 100, 2)
        return freq_dict

