from typing import Dict, Tuple, Union
from Bio import SeqIO
from Bio.SeqUtils import GC

def filter_fastq(file_path: str,
                 gc_bounds: Union[float, Tuple[Union[float, int], Union[float, int]]] = (0, 100),
                 length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
                 quality_threshold: Union[float, int] = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filters sequences based on GC-content, length, and quality threshold.

    Parameters:
    -----------
    file_path (str):
        Path to the FastQ file.
    gc_bounds (Union[float, Tuple[Union[float, int], Union[float, int]]]):
        GC-content filtering bounds.
        - Default is 100 (no filtering).
        - If a single float is provided, it's considered as the upper bound.
        - If an integer is provided, it's considered as the upper bound.
        - If a tuple of two numbers is provided, they are treated as the lower and upper bounds.
    length_bounds (Union[int, Tuple[int, int]]):
        Length filtering bounds.
        - Default is (0, 2^32) (no filtering).
        - If a single integer is provided, it's considered as the upper bound.
        - If a tuple of two integers is provided, they are treated as the lower and upper bounds.
    quality_threshold (Union[float, int]):
        Quality threshold value.
        - Default is 0 (no filtering).
        - Sequences with an average quality below the threshold are discarded.

    Returns:
    --------
    Dict[str, Tuple[str, str]]:
        A filtered dictionary where each value is a tuple containing two strings,
        representing DNA sequences and their quality.
    """
    filtered_seqs = {}

    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            gc_content = GC(record.seq)

            if isinstance(gc_bounds, (float, int)):
                lower_bound = 0
                upper_bound = gc_bounds
            elif isinstance(gc_bounds, tuple) and len(gc_bounds) == 2:
                lower_bound, upper_bound = gc_bounds
            else:
                raise ValueError("gc_bounds should be a single float/int or a tuple of two floats/ints.")

            if lower_bound <= gc_content <= upper_bound:
                seq_length = len(record.seq)

                if isinstance(length_bounds, int):
                    lower_bound = 0
                    upper_bound = length_bounds
                elif isinstance(length_bounds, tuple) and len(length_bounds) == 2:
                    lower_bound, upper_bound = length_bounds
                else:
                    raise ValueError("length_bounds should be a single integer or a tuple of two integers.")

                if lower_bound <= seq_length <= upper_bound:
                    quality_values = [chr(qual) for qual in record.letter_annotations["phred_quality"]]
                    quality_count = sum(record.letter_annotations["phred_quality"])
                    seq_quality = quality_count / seq_length

                    if seq_quality >= quality_threshold:
                        filtered_seqs[record.id] = (str(record.seq), "".join(quality_values))

    return filtered_seqs
