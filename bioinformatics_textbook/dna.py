"""DNA.py

A module for working with DNA sequences through the DNA class
"""

from typing import Optional

from bioinformatics_textbook.ch01.ch01 import compute_hamming_distance


class DNA(str):
    """Representation of a DNA sequence. Contains methods for manipulating a single DNA sequence."""

    def __init__(self, sequence: str) -> None:
        """Initialize the DNA sequence object

        :param sequence: DNA sequence
        :type sequence: str
        """
        super().__init__()
        self.seq = sequence
        self._complementation_table = self._make_dna_complementation_table()

    def generate_d_neighborhood(
        self, num_allowed_mismatches: int, seq: Optional[str] = None
    ) -> list:
        """Find the set of all k-mers whose Hamming distance from a sequence does not exceed a set maximum (d).

        :param num_allowed_mismatches: The maximum allowed Hamming distance (i.e. the maximum number of allowed mismatches).
        :type num_allowed_mismatches: int
        :param seq: Sequence from which a neighborhood is generated. If None, the object sequence is used., defaults to None
        :type seq: Optional[str], optional
        :return: k-mers that are in the d neighborhood of the sequence.
        :rtype: list
        """
        seq = seq if seq is not None else self.seq

        nucleotides = {"A", "C", "G", "T"}

        # if there are no allowed mismatches the seq's neighborhood is itself
        if num_allowed_mismatches == 0:
            return [seq]

        # if the seq is only a single nucleotide, its neighborhood is the set of all nucleotides
        # this  serves as the base condition for the recursion
        if len(seq) == 1:
            return list(nucleotides)

        # making a neighborhood a set instead of a list means duplicates do not have to be explicitly removed
        neighborhood = set()
        suffix_neighbors = self.generate_d_neighborhood(
            num_allowed_mismatches, self._slice_suffix(seq)
        )
        for suffix_neighbor in suffix_neighbors:
            # the suffix neighbors must have a Hamming distance to seq that is greater than or equal to the max allowed Hamming distance
            if (
                compute_hamming_distance(self._slice_suffix(seq), suffix_neighbor)
                < num_allowed_mismatches
            ):
                for nucleotide in nucleotides:
                    neighborhood.add(nucleotide + suffix_neighbor)
            else:
                neighborhood.add(self._slice_first_nucleotide(seq) + suffix_neighbor)

        return list(neighborhood)

    def reverse_complement(self, seq: Optional[str] = None) -> "DNA":
        """Reverse and complement a DNA sequence

        :return: Reverse complemented DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        comp = self.complement(seq)
        rev_comp = self.reverse(comp)

        return DNA(rev_comp)

    def complement(self, seq: Optional[str] = None) -> "DNA":
        """Complement a DNA sequence

        :return: Complemented DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq.translate(self._complementation_table))

    def reverse(self, seq: Optional[str] = None) -> "DNA":
        """Reverse a DNA sequence

        :return: Reversed DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq[::-1])

    def _slice_suffix(self, seq: Optional[str] = None) -> "DNA":
        """Get the suffix of a string

        :param pattern: The string to slice
        :type pattern: str
        :return: The suffix
        :rtype: str
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq[1:])

    def _slice_first_nucleotide(self, seq: Optional[str] = None) -> "DNA":
        """Get the first nucleotide of a DNA string

        :param pattern: The DNA string to slice
        :type pattern: str
        :return: The first nucleotide base
        :rtype: str
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq[:1])

    def _make_dna_complementation_table(self) -> dict:
        """Construct a translation table for complementing DNA

        :return: DNA complementation table
        :rtype: dict
        """
        comp_table = str.maketrans("ATCG", "TAGC")

        return comp_table
