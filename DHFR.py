from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np

class BiologicalSequence:
    def __init__(self, sequence):
        self.seq = Seq(sequence)

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, item):
        return self.seq[item]

    def __str__(self):
        return str(self.seq)

    def gc_content(self):
        return GC(self.seq)

    def reverse_complement(self):
        return self.seq.reverse_complement()


class NucleicAcidSequence(BiologicalSequence):
    def complement(self):
        return self.seq.complement()

    def transcribe(self):
        raise NotImplementedError("This method is specific to DNA sequences.")


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        return self.seq.transcribe()


class RNASequence(NucleicAcidSequence):
    def translate(self):
        return self.seq.translate()


class AminoAcidSequence(BiologicalSequence):
    def molecular_weight(self):
        weight_dict = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}
        return sum(weight_dict[aa] for aa in self.seq if aa in weight_dict)


class FastQFilter:
    def __init__(self, input_path, output_path, min_length=0, max_length=float('inf'), min_quality=0, gc_bounds=(0, 100)):
        self.input_path = input_path
        self.output_path = output_path
        self.min_length = min_length
        self.max_length = max_length
        self.min_quality = min_quality
        self.gc_bounds = gc_bounds

    def filter(self):
        with open(self.output_path, 'w') as output_handle:
            for record in SeqIO.parse(self.input_path, "fastq"):
                seq_len = len(record.seq)
                avg_quality = np.mean(record.letter_annotations["phred_quality"])
                gc_content = GC(record.seq)

                if (self.min_length <= seq_len <= self.max_length and
                    self.min_quality <= avg_quality and
                    self.gc_bounds[0] <= gc_content <= self.gc_bounds[1]):
                    SeqIO.write(record, output_handle, "fastq")

