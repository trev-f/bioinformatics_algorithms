"""DNA.py

A module for working with DNA sequences through the DNA class
"""

from typing import Optional

class DNA(str):
    """Representation of a DNA sequence. Contains methods for manipulating a single DNA sequence.
    """

    def __init__(self, sequence: str) -> None:
        """Initialize the DNA sequence object

        :param sequence: DNA sequence
        :type sequence: str
        """
        self.seq = sequence
        self._complementation_table = self._make_dna_complementation_table()
    
    def reverse_complement(self, seq: Optional[str] = None) -> 'DNA':
        """Reverse and complement a DNA sequence

        :return: Reverse complemented DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        comp = self.complement(seq)
        rev_comp = self.reverse(comp)

        return DNA(rev_comp)
    
    def complement(self, seq: Optional[str] = None) -> 'DNA':
        """Complement a DNA sequence

        :return: Complemented DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq.translate(self._complementation_table))
    
    def reverse(self, seq: Optional[str] = None) -> 'DNA':
        """Reverse a DNA sequence

        :return: Reversed DNA sequence
        :rtype: DNA
        """
        seq = seq if seq is not None else self.seq

        return DNA(seq[::-1])
    
    def _make_dna_complementation_table(self) -> dict:
        """Construct a translation table for complementing DNA

        :return: DNA complementation table
        :rtype: dict
        """
        comp_table = str.maketrans("ATCG", "TAGC")

        return comp_table
